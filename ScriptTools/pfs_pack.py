# Purpose:  Copy a pfs files and all files it references into a separate _stage directory. This can be used for sharing a model including all required files
#           but not the entire directory with result files etc.
#           - Made for .she files, but should work to a certain extent for other pfs files.
#           - Pfs files referenced in the main pfs file will also be processed.
#           - Referenced m1dx files will also be processed to some extent, also any .mupp and .sqlte files with the same base name as the m1dx file will be included.
#           - Referenced .shp files will bring their buddies along
#           - Missing files will be reported, but not treated as errors
#           - References to files outside the directory of the main pfs files are not recommended but should be handled correctly by placing the
#             main pfs file in a sub[-sub[-sub[-sub]]] directory of the staging directory.
#             Don't know what happens when you reference files on a different drive - try and report (if you still can after trying)!
#           - .dll files will not be included
# Requires: Python, mikeio, optionally tkinter (for file selection dialog)
# Usage:    A: Run from console with or without arguments. For non-interactive mode specify all required arguments
#           B: Double click, navigate to pfs file... Requires that .py files are associated with a python interpreter on your system. tkinter installed with python gives you GUI file selection.
#           C: Register as shell command, then right click any pfs file in Windows Explorer and select the command (Windows only)
#              How? Save the following lines as a text file with .reg extension. Set the paths to your python installation and to where you placed this script. Then execute as administrator.
#              (or edit the registry manually to add the below keys)
#                Windows Registry Editor Version 5.00

#                [HKEY_CLASSES_ROOT\*\shell\Stage PFS]
#                @="Stage PFS File"

#                [HKEY_CLASSES_ROOT\*\shell\Stage PFS\Command]
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

tab_file = 0

def is_pfs_file(path: Path) -> bool:
  try:
    _ = mikeio.pfs.read_pfs(path)
    return True
  except Exception:
    return False


def extract_file_references(pfs_obj, only_if_used: bool) -> List[Path]:
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
    if only_if_used and section_dict.get("IsDataUsedInSetup") == 0: # If tag does not exist we do want to continue processing!
      return
    for key, val in section_dict.items():
      dispatch_node(val, section_dict, breadcrumb + [key])

  def process_list(lst, section_dict, breadcrumb):
    for i, item in enumerate(lst):
      dispatch_node(item, section_dict, breadcrumb + [f"[{i}]"])

  def process_value(val, section_dict, breadcrumb):
    if isinstance(val, str):
      match = re.fullmatch(r'\|(.+)\|', val.strip())
      if match:
        filepath = Path(match.group(1))
        if filepath:
          file_paths.append((filepath, breadcrumb.copy()))

  dispatch_node(pfs_obj.to_dict(), [], [])
  return file_paths


def find_m1dx_related_files(m1dx_path: Path) -> List[Path]:
  related = []
  stop = [False]

  def strip_ns(tag):
    return tag.split('}')[-1] if '}' in tag else tag

  def walk(elem):
    # We are not interested in output files, just input files required to run the model
    if strip_ns(elem.tag) == "ResultSpecifications":
      return

    if strip_ns(elem.tag) == "Path" and elem.text:
      text = elem.text.strip()
      if text and not text.lower().endswith('.html'):
        path = (m1dx_path.parent / Path(text)).resolve()
        related.append(path)
    for child in elem:
      walk(child)
  try:
    tree = ET.parse(m1dx_path)
    walk(tree.getroot())
  except Exception as e:
    print(f"Error parsing {m1dx_path}: {e}")

  # Add sidecar files (e.g. .sqlite, .mupp) with same stem
  for ext in [".sqlite", ".mupp"]:
    sidecar = m1dx_path.with_suffix(ext)
    if sidecar.exists():
      related.append(sidecar)
  return related


def find_shp_related_files(shp_path: Path) -> List[Path]:
  base = shp_path.with_suffix('')
  extensions = ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']
  return [base.with_suffix(ext) for ext in extensions if base.with_suffix(ext).exists()]


def collect_files(pfs_path: Path, mode: str, collected: Set[Path]) -> None:
  global tab_file
  pfs_obj = mikeio.pfs.read_pfs(pfs_path)
  collected.add(pfs_path.resolve())
  for ref, breadcrumb in extract_file_references(pfs_obj, only_if_used=True):
    abs_ref = (pfs_path.parent / ref).resolve()
    if not abs_ref.exists():
      print(tab_file * "  " + f"Missing file: {abs_ref} ({'/'.join(breadcrumb)})")
      continue
    if abs_ref.suffix.lower() == ".m1dx":
      collected.add(abs_ref)
      for f in find_m1dx_related_files(abs_ref):
        if not f.exists():
          print(tab_file * "  " + f"M1D-Missing file: {f}")
          continue
        collected.add(f.resolve())
    elif abs_ref.suffix.lower() == ".shp":
      for shp_related in find_shp_related_files(abs_ref):
        collected.add(shp_related.resolve())
    elif abs_ref.suffix.lower() == ".dll":
      # Only case so far: The python installation dll referenced in the she file. No good to include that!
      continue
    elif abs_ref.suffix.lower() not in ['.dfs0', '.dfs2', '.dfs3'] and is_pfs_file(abs_ref):
      if not abs_ref in collected:
        collected.add(abs_ref)
        print(tab_file * "  " + f"Sub-PFS-File: {abs_ref}")
      tab_file += 2
      collect_files(abs_ref, mode, collected)
      tab_file -= 2
    else:
      collected.add(abs_ref)
  # print(collected)


def create_zip(master_pfs: Path, mode: str, out_path: str, force):
  master_pfs = master_pfs.resolve()

  if not is_pfs_file(master_pfs):
    raise ValueError(f'The file provided does not seem to be a valid pfs file:\n\t"{master_pfs}"')

  # Collect all referenced files including the main file
  all_paths = {master_pfs}
  collect_files(master_pfs, mode, all_paths)

  # Compute common root across all files
  real_root = Path(os.path.commonpath([str(p) for p in all_paths]))
  
  # creating zip file
  if out_path is None:
    zip_file = master_pfs.parent / (master_pfs.stem + ".zip")
  
  if os.path.isfile(zip_file) and not force:
    answer = input(f"{zip_file} exists - overwrite? (y/n): ").strip().lower()
    if answer not in ["y", "yes"]:
      print("Skip writing output file")
      return
  with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
    for file in all_paths:
      try:
        relative_path = file.relative_to(real_root)
        if relative_path == Path("."):
          relative_path = file.name
      except ValueError:
        continue  # skip if file is outside tree
      if file.is_file():
        # shutil.copy2(file, dest)
        zipf.write(file, relative_path)
  print(f"Zip file created: {zip_file}")


def parse_args():
  parser = argparse.ArgumentParser(description="Stage a PFS file and its dependencies.") 
  parser.add_argument("pfs_path", nargs='?', type=str, help="Path to the main PFS file")
  parser.add_argument("-o", "--output", type=str, help="Path to the output zip file")
  parser.add_argument("-m", "--minimum", action="store_true", help="Inlcude only currently used files")
  parser.add_argument("-a", "--all", action="store_true", help="Include all files, even if not currently used")
  parser.add_argument("-f", "--force", action="store_true", help="Overwrite output if it already exists")

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
      print("tkinter not available. Please enter the full path to the main PFS file:")
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
    print(sys.argv)
    args = parse_args()
    mode = "minimum" if args.minimum else "all"
    create_zip(Path(args.pfs_path), mode, args.output, args.force)
  except Exception as e:
    traceback.print_exc()
  input("\nPress Enter to exit...")
