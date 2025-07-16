# Purpose:  Pack a pfs file and all files it references into a zip archive (default) or a staging directory. This can e. g.
#           be used for sharing a model including all required files but not the entire directory with result files etc.
#           - Made for .she files, but should work to a certain extent for other pfs files.
#           - Pfs files referenced in the main pfs file will also be processed.
#           - Referenced .m1dx files will also be processed to some extent, also any .mupp, .sqlite and .cs files
#             with the same base name as the .m1dx file will be included.
#           - Referenced .shp files will bring their buddies along
#           - Missing files will be reported, but not treated as errors
#           - References to files outside the directory of the main pfs files are not recommended but are handled
#             by placing the main pfs file in a sub[-sub[-sub[-sub]]] directory of the staging directory.
#             Don't know what happens when you reference files on a different drive - try and report!
#           - .dll files will not be included
# Requires: mikeio, optionally tkinter (for file selection dialog)
# Usage:    A: Run from console with or without arguments. For non-interactive mode specify all required arguments.
#           B: Double click, navigate to pfs file... Requires that .py files are associated with a python interpreter on your system. tkinter installed with python gives you GUI file selection.
#           C: Register as shell command, then right click any pfs file in Windows Explorer and select the command (Windows only)
#              How? Save the following lines as a text file with .reg extension. Set the paths to your python installation and to where you placed this script. Then execute as administrator.
#              (or edit the registry manually to add the below keys)
#                Windows Registry Editor Version 5.00

#                [HKEY_CLASSES_ROOT\*\shell\Pack PFS]
#                @="Pack PFS File"

#                [HKEY_CLASSES_ROOT\*\shell\Pack PFS\Command]
#                @="\"C:\\Program Files\\Python313\\python.exe\" \"C:\\Users\\ME_ME_ME\\scripts\\pfs_pack.py\" \"%1\""
# Warranty: None whatsoever.
# Author: DHI\uha
# Date:   05/2025

from pathlib import Path
from typing import Set, List
import argparse
import mikeio.pfs
import os
import re
import shutil
import sys
import traceback
import xml.etree.ElementTree as ET
import zipfile

tab_file = 0 # printing indentation
IGNORE_BREADCRUMBS = [
  ["MIKESHE_FLOWMODEL", "Results_Post_Processing", "RunStatistics"],
  ["MIKESHE_WaterBalance", "Extraction_Step"],
  ["MIKESHE_WaterBalance", "Postprocessing_Steps"]
]

def is_pfs_file(path: Path) -> bool:
  try:
    _ = mikeio.pfs.read_pfs(path)
    return True
  except Exception as e:
    # Some extensions are known to be pfs files, it is an error if those cannot be parsed!
    if os.path.splitext(path)[1].lower() in [".etv", ".uzs", ".wel", ".wbl", ".sheres", ".mhydro", ".ecolab"]:
      raise
    return False


def extract_file_references(pfs_obj, include_all: bool) -> List[Path]:
  file_paths = []

  def dispatch_node(node, section_dict, breadcrumb):
    if isinstance(node, mikeio.pfs._pfssection.PfsSection):
      dispatch_node(node.to_dict(), section_dict, breadcrumb)
    elif isinstance(node, dict):
      process_section(node, breadcrumb)
    elif isinstance(node, list):
      process_list(node, section_dict, breadcrumb)
    else:
      process_value(node, section_dict, breadcrumb)

  def process_section(section_dict, breadcrumb):
    if breadcrumb in IGNORE_BREADCRUMBS:
      return
    if not include_all:
      if section_dict.get("IsDataUsedInSetup") == 0: # If tag does not exist we do want to continue processing!    
        return
    for key, val in section_dict.items():
      dispatch_node(val, section_dict, breadcrumb + [key])

  def process_list(lst, section_dict, breadcrumb):
    for i, item in enumerate(lst):
      dispatch_node(item, section_dict, breadcrumb + [f"[{i}]"])

  def process_value(val, section_dict, breadcrumb):
    if breadcrumb in IGNORE_BREADCRUMBS:
      return
    if isinstance(val, str):
      match = re.fullmatch(r'\|(.+)\|', val.strip())
      if match:
        filepath = Path(match.group(1))
        if filepath:
          file_paths.append((filepath, breadcrumb.copy()))

  dispatch_node(pfs_obj, [], [])
  return file_paths


def find_m1dx_related_files(m1dx_path: Path) -> List[Path]:
  related = []

  def strip_ns(tag):
    return tag.split('}')[-1] if '}' in tag else tag

  def walk(elem):
    # We are not interested in output files, just input files required to run the model
    if strip_ns(elem.tag) == "ResultSpecifications":
      return

    if strip_ns(elem.tag) == "Path" and elem.text:
      text = elem.text.strip()
      if text and not text.lower().endswith('.html'):
        path = (m1dx_path.parent / Path(text))
        related.append(path)
    for child in elem:
      walk(child)
  try:
    tree = ET.parse(m1dx_path)
    walk(tree.getroot())
  except Exception as e:
    print(f"Error parsing {m1dx_path}: {e}")

  # Add sidecar files with same stem
  # - .sqlite, .mupp: editable setup files
  # - .cs:            plugins
  for ext in [".sqlite", ".mupp", ".cs"]:
    sidecar = m1dx_path.with_suffix(ext)
    if sidecar.exists():
      related.append(sidecar)
  return related


def find_shp_related_files(shp_path: Path) -> List[Path]:
  base = shp_path.with_suffix('')
  extensions = ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']
  return [base.with_suffix(ext) for ext in extensions if base.with_suffix(ext).exists()]


def collect_files(pfs_path: Path, include_all: bool, collected: Set[Path]) -> None:
  global tab_file
  pfs_obj = mikeio.pfs.read_pfs(pfs_path)
  collected.add(pfs_path)
  for ref, breadcrumb in extract_file_references(pfs_obj, include_all):    
    abs_ref = (pfs_path.parent / ref)
    if not abs_ref.exists():
      print(tab_file * "  " + f"Missing file: {abs_ref} ({'/'.join(breadcrumb)})")
      continue
    if abs_ref.suffix.lower() == ".m1dx":
      collected.add(abs_ref)
      for f in find_m1dx_related_files(abs_ref):
        if not f.exists():
          print(tab_file * "  " + f"M1D-Missing file: {f}")
          continue
        collected.add(f)
    elif abs_ref.suffix.lower() == ".shp":
      for shp_related in find_shp_related_files(abs_ref):
        collected.add(shp_related)
    elif abs_ref.suffix.lower() == ".hot":
      collected.add(abs_ref)
      collected.add(Path(os.path.splitext(abs_ref)[0] + ".frf")) # special case - whenever mshe needs a .hot file it also needs the .frf file, not listed in .sheres!
    elif abs_ref.suffix.lower() == ".dll":
      # Only case so far: The python installation dll referenced in the .she file. No good to include that!
      continue
    elif abs_ref.suffix.lower() not in ['.dfs0', '.dfs2', '.dfs3'] and is_pfs_file(abs_ref):
      if not abs_ref in collected:
        collected.add(abs_ref)
        print(tab_file * "  " + f"Sub-PFS-File: {abs_ref}")
      tab_file += 2
      collect_files(abs_ref, include_all, collected)
      tab_file -= 2
    else:
      collected.add(abs_ref)
  # print(collected)


def main(master_pfs: Path, include_all: bool, to_zip: bool, out_path: str, force: bool):
  def pack():
    for file in all_paths:
      try:
        relative_path = file.relative_to(real_root)
        if relative_path == Path("."):
          relative_path = file.name
      except ValueError:
        continue  # skip if file is outside tree
      if file.is_file():
        if to_zip:
          zipf.write(file, relative_path)
        else:
          dest = out_path / relative_path
          dest.parent.mkdir(parents=True, exist_ok=True)
          shutil.copy2(file, dest)
    if to_zip:
      print(f"Zip file created: {out_path}")
    else:
      print(f"Staging directory created: {out_path}")

  master_pfs = Path(master_pfs) # in case a string is passed, otherwise it does an unneccessary copy

  if not is_pfs_file(master_pfs):
    raise ValueError(f'The file provided does not seem to be a valid pfs file:\n\t"{master_pfs}"')

  # Collect all referenced files including the main file
  all_paths = {master_pfs}
  collect_files(master_pfs, include_all, all_paths)

  # Compute common root across all files
  real_root = Path(os.path.commonpath([str(p) for p in all_paths]))

  if to_zip:
    # creating zip file
    if out_path is None:
      out_path = master_pfs.parent / (master_pfs.stem + ".zip")
    
    if os.path.isfile(out_path) and not force:
      answer = input(f"{out_path} exists - overwrite? (y/n): ").strip().lower()
      if answer not in ["y", "yes"]:
        print("Skip writing output file")
        return
    with zipfile.ZipFile(out_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
      pack()
  else:
    if out_path is None:
      out_path = master_pfs.parent / (master_pfs.stem + "_staging")
    if os.path.isdir(out_path) and not force:
      answer = input(f"{out_path} exists - overwrite? (y/n): ").strip().lower()
      if answer not in ["y", "yes"]:
        print("Skip writing output directory")
        return
    pack()


def parse_args():
  parser = argparse.ArgumentParser(description="Pack a PFS file and its dependencies.") 
  parser.add_argument("pfs_path", nargs='?', type=str, help="Path to the main PFS file")
  parser.add_argument("-o", "--output", type=str, help="Path to the output zip file or directory")
  parser.add_argument("-m", "--minimum", action="store_true", help="Include only currently used files")
  parser.add_argument("-a", "--all",     action="store_true", help="Include all files, even if not currently used (default)")
  parser.add_argument("-f", "--force",   action="store_true", help="Overwrite output if it already exists")
  parser.add_argument("-d", "--dir",     action="store_true", help="Create a directory as output instead of a zip file")

  args = parser.parse_args()

  # Prompt for path if not provided
  if not args.pfs_path:
    try:
      import tkinter as tk
      from tkinter import filedialog
      root = tk.Tk()
      root.withdraw()
      path = filedialog.askopenfilename(title="Select main PFS file")
      if not path:
        raise ValueError("No file selected. Exiting.")
      args.pfs_path = path
    except ModuleNotFoundError:
      print("Cannot open file selection dialog (tkinter not available). Please enter the full path to the main PFS file:")
      path = input("Path: ").strip()
      if not path:
        raise ValueError("No path provided. Exiting.")
      args.pfs_path = path

  # Ask for mode if neither flag was set
  if not args.minimum and not args.all:
    while True:
      answer = input("Include all files, even if not currently used? (y/n): ").strip().lower()
      if answer in ["y", "yes"]:
        args.all = True
        break
      elif answer in ["n", "no"]:
        args.minimum = True
        break

  return args


if __name__ == "__main__":
  import sys
  try:
    # print(sys.argv)
    args = parse_args()
    include_all = not args.minimum
    to_zip = not args.dir
    main(Path(args.pfs_path), include_all, to_zip, args.output, args.force)
  except Exception as e:
    traceback.print_exc()
  input("\nPress Enter to exit...")
