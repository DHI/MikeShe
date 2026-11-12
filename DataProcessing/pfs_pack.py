# Purpose:  Pack a pfs file and the files it references into a zip archive (default) or a staging directory. This can e. g.
#           be used for sharing a model including all required files but not the entire directory with result files etc.
#           - Made for .she files, but should work to a certain extent for other pfs files.
#           - By default only the files in use are packed (-m): references in a section with IsDataUsedInSetup = 0
#             are left out. -a packs all referenced files.
#           - Pfs files referenced in the main pfs file will also be processed.
#           - Referenced .m1dx files will also be processed to some extent, also any .mupp, .sqlite and .cs files
#             with the same base name as the .m1dx file will be included.
#           - Referenced .shp files will bring their buddies along
#           - Missing files will be reported, but not treated as errors
#           - References to files outside the directory of the main pfs files are not recommended but are handled
#             by placing the main pfs file in a sub[-sub[-sub[-sub]]] directory of the zip file or staging
#             directory (-d). With -p the main pfs file lands directly in the output directory instead, and such
#             references land above it.
#             Files on a different drive than the main pfs file cannot be packed: the script lists them and
#             stops, unless you confirm (or pass -i) to pack without them.
#           - .dll files will not be included
#           - A staging directory is filled like robocopy does it: a file that is already there with the same size
#             and time is skipped, a differing file is updated from the source after confirmation (or with -f), other
#             files in the directory are kept, and a re-run of the same command continues an interrupted copy.
#             On Windows robocopy does the copying (retries and resumes after a network error).
# Requires: mikeio, optionally tkinter (for file selection dialog)
# Usage:    A: Run from console with or without arguments. For non-interactive mode specify all required arguments
#              (pfs file, -m or -a, -o with -p, -f if the target already exists, and -i if files on another
#              drive may be left out). Without -m or -a the script waits for Enter at the end, so that the
#              console window of usage B and C stays open: a batch file has to pass one of them.
#              Exit code: 0 on success, 1 on failure, 2 on a command line error.
#           B: Windows Explorer: double click the script, or right click it and choose "Open with" and a python interpreter.
#              Then navigate to the pfs file. The double click requires that .py files are associated with a python
#              interpreter on your system. tkinter installed with python gives you GUI file selection.
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

from datetime import datetime
from pathlib import Path
from typing import Set, List, Tuple
import argparse
import mikeio.pfs
import os
import re
import shutil
import subprocess
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


class UserAbort(Exception):
  """The user gave no input where one was needed (empty prompt, cancelled dialog). Reported without a traceback."""


class UserError(Exception):
  """A problem with the input that the user can fix, with a message written for the user. Reported without a
  traceback, which stays reserved for defects of the script."""


def canonical(path) -> Path:
  """The one form every path takes before it enters the packed set: absolute, no "..", native separators.
  The same file referenced from different directories is then one set member, and the common root is well
  defined. normpath rather than resolve(): resolve() would follow junctions and symbolic links, moving the
  file out of the model tree the user sees."""
  return Path(os.path.normpath(Path(path).absolute()))


def try_read_pfs(path: Path):
  """The parsed pfs document, or None if the file is not a pfs file. Parsing is the expensive part of this
  script, so the result is handed on to collect_files() instead of being parsed a second time there."""
  try:
    return mikeio.pfs.read_pfs(path)
  except Exception as e:
    # Some extensions are known to be pfs files, it is an error if those cannot be parsed!
    if os.path.splitext(path)[1].lower() in [".she", ".etv", ".uzs", ".wel", ".wbl", ".sheres", ".mhydro", ".ecolab"]:
      raise UserError(f'Cannot parse pfs file "{path}": {e}') from e
    return None


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
        related.append(canonical(m1dx_path.parent / text))
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
  """The referenced .shp (which exists, the caller checked) plus the companions present next to it.
  The companion extensions follow the case of the referenced extension, so that RIV.SHP finds RIV.DBF on a
  case-sensitive file system. The stem is kept as is: on such a system A.shp and a.shp are different files."""
  base, shp_ext = os.path.splitext(str(shp_path)) # splitext keeps a dotted stem such as "riv.v2" intact
  companion_extensions = ['.shx', '.dbf', '.prj', '.cpg', '.qpj']
  if shp_ext.isupper():
    companion_extensions = [ext.upper() for ext in companion_extensions]
  related = [shp_path]
  related.extend(Path(base + ext) for ext in companion_extensions if os.path.isfile(base + ext))
  return related


def collect_files(pfs_path: Path, pfs_obj, include_all: bool, collected: Set[Path], walked: Set[Path]) -> None:
  """Adds pfs_path (already parsed into pfs_obj) and everything it references, recursively, to `collected`
  (the files to pack). `walked` holds the pfs files whose references have already been extracted."""
  global tab_file
  collected.add(pfs_path)
  walked.add(pfs_path)
  for ref, breadcrumb in extract_file_references(pfs_obj, include_all):    
    abs_ref = canonical(pfs_path.parent / ref)
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
      frf = Path(os.path.splitext(abs_ref)[0] + ".frf") # special case - whenever mshe needs a .hot file it also needs the .frf file, not listed in .sheres!
      if frf.is_file():
        collected.add(frf)
      else:
        print(tab_file * "  " + f"Missing file: {frf} (needed together with {abs_ref.name})")
    elif abs_ref.suffix.lower() == ".dll":
      # Only case so far: The python installation dll referenced in the .she file. No good to include that!
      continue
    elif abs_ref in walked:
      # A file that is referenced twice, or refers back to its parent, is walked (and parsed) once.
      continue
    elif abs_ref.suffix.lower() not in ['.dfs0', '.dfs2', '.dfs3'] and (sub_pfs := try_read_pfs(abs_ref)) is not None:
      print(tab_file * "  " + f"Sub-PFS-File: {abs_ref}")
      tab_file += 2
      collect_files(abs_ref, sub_pfs, include_all, collected, walked)
      tab_file -= 2
    else:
      collected.add(abs_ref)
  # print(collected)


def confirm_overwrite(question: str, force: bool, flag_hint: str = "-f to force update") -> bool:
  """True if the flag that answers the question with yes was given (`force`), or the user answered yes at the
  prompt. `flag_hint` names that flag for the case that nobody can answer."""
  if force:
    return True
  try:
    answer = input(f"{question} (y/n): ").strip().lower()
  except (EOFError, OSError): # no usable stdin (CI, scheduled task, "< NUL"): stop with the cause instead of a traceback
    raise UserAbort(f"\nNo stdin available to answer the question, run again with {flag_hint}.")
  return answer in ["y", "yes"]


def describe(path: Path) -> str:
  """Size and modification time of a file, for the overview of what a copy would replace."""
  stat = path.stat()
  return f"{stat.st_size:,} bytes, {datetime.fromtimestamp(stat.st_mtime):%Y-%m-%d %H:%M:%S}"


def is_identical(src: Path, dest: Path) -> bool:
  """Same size and modification time, which both copy routines preserve: the file has been copied completely.
  A partial file from an interrupted copy never matches, its time is the time of the interruption. The two
  second tolerance is the FAT time stamp resolution, the same tolerance robocopy applies with /FFT."""
  src_stat, dest_stat = src.stat(), dest.stat()
  return src_stat.st_size == dest_stat.st_size and abs(src_stat.st_mtime - dest_stat.st_mtime) < 2


def copy_files(jobs: List[Tuple[Path, Path]]) -> bool:
  """Copies every (source, destination) pair with its time stamps. True on success. On Windows robocopy does it,
  one call per (source directory, destination directory) pair with the file names as arguments: /Z resumes a copy
  that robocopy itself interrupted, /R and /W ride out a transient network error, and exit codes below 8 mean
  success (1: files copied). Elsewhere a plain copy, which a re-run of the script continues file by file."""
  if os.name == "nt" and shutil.which("robocopy"):
    groups = {}
    for src, dest in jobs:
      groups.setdefault((src.parent, dest.parent), []).append(src.name)
    for (src_dir, dest_dir), names in groups.items():
      command = ["robocopy", str(src_dir), str(dest_dir), *names, "/Z", "/R:20", "/W:30", "/NJH", "/NJS", "/NDL"]
      sys.stdout.flush() # keeps the plan before robocopy's lines when stdout is a pipe
      exit_code = subprocess.run(command).returncode
      if exit_code >= 8:
        print(f"\nrobocopy failed with exit code {exit_code} while copying into {dest_dir}")
        print("Re-run the same command to continue where it stopped (-f skips the confirmation).")
        return False
    return True
  for src, dest in jobs:
    print(f"{src} ({describe(src)})")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
  return True


def pack_entries(all_paths: Set[Path], real_root: Path, place: bool) -> List[Tuple[Path, Path]]:
  """(file, path relative to real_root) for every file to pack, in path order. With `place` the relative path
  may start with "..". A path that is not a file is reported and left out."""
  entries = []
  for file in sorted(all_paths):
    relative_path = Path(os.path.relpath(file, real_root)) if place else file.relative_to(real_root)
    if relative_path == Path("."):
      relative_path = Path(file.name)
    if not file.is_file():
      print(f"Skipped, not a file: {file}")
      continue
    entries.append((file, relative_path))
  return entries


def fill_staging_dir(entries: List[Tuple[Path, Path]], out_path: Path, force: bool, notes: List[str]) -> bool:
  """Brings the staging directory out_path up to date with the entries of pack_entries(). Differing files that
  the user wants to keep are added to `notes`. False when the copy failed."""
  new = []       # (source, destination) of the files that are not there yet
  differing = [] # (source, destination) of the files that are there but differ; updating them needs -f or a yes
  n_complete = 0
  for file, relative_path in entries:
    dest = canonical(out_path / relative_path)
    if relative_path.parts[0] == "..":
      print(f"Lands above the output directory: {dest}")
    if not dest.exists():
      new.append((file, dest))
    elif is_identical(file, dest):
      n_complete += 1
    else:
      differing.append((file, dest))

  # The plan, then the copy. Nothing is deleted, ever.
  print(f"  {n_complete} file(s) complete, skipped")
  print(f"  {len(new)} file(s) new, will be copied")
  jobs = new
  if differing:
    print(f"  {len(differing)} file(s) differ, will be updated:")
    for file, dest in differing:
      print(f"    {dest}")
      print(f"      here:   {describe(dest)}")
      print(f"      source: {describe(file)}")
    if confirm_overwrite("Update the differing file(s)?", force):
      jobs = new + differing
    else:
      notes.append(f"{len(differing)} file(s) kept as they are")
  return not jobs or copy_files(jobs)


def main(master_pfs: Path, include_all: bool, to_zip: bool, place: bool, out_path: str, force: bool, incomplete: bool) -> bool:
  """Returns True when the package was written, False when the user declined to overwrite the zip file or to
  pack without the files on another drive, or the copy failed."""
  def confirm_incomplete() -> bool:
    """True if nothing has to be left out, or the user accepts (-i or a yes) a package without the files on
    another drive."""
    if not other_drive:
      return True
    print("Referenced file(s) on a different drive than the main pfs file cannot be packed:")
    for file in other_drive:
      print(f"  {file}")
    if confirm_overwrite("Pack without these file(s)? The package will be incomplete.", incomplete, "-i to pack without these files"):
      return True
    print("Nothing packed")
    return False

  master_pfs = canonical(master_pfs)

  if not master_pfs.is_file():
    raise UserError(f'File not found:\n\t"{master_pfs}"')
  master_obj = try_read_pfs(master_pfs)
  if master_obj is None:
    raise UserError(f'The file provided does not seem to be a valid pfs file:\n\t"{master_pfs}"')

  # Collect all referenced files including the main file
  if not include_all:
    print("Packing only the files in use. Use -a to include all referenced files.")
  all_paths = {master_pfs}
  collect_files(master_pfs, master_obj, include_all, all_paths, walked=set())

  # A path on a different drive (or UNC share) than the main pfs file has no place in a relative tree. It
  # leaves the set here; a file among them makes the package incomplete, which confirm_incomplete() asks about
  # once the output path is known to be valid.
  master_drive = os.path.splitdrive(master_pfs)[0].lower()
  off_drive = {p for p in all_paths if os.path.splitdrive(p)[0].lower() != master_drive}
  all_paths -= off_drive
  other_drive = sorted(p for p in off_drive if p.is_file())
  for path in sorted(off_drive.difference(other_drive)):
    print(f"Skipped, not a file: {path}")

  if place:
    # The master lands directly in the output directory, every other file keeps its position relative to it
    real_root = master_pfs.parent
  else:
    real_root = Path(os.path.commonpath([str(p) for p in all_paths]))

  # What keeps the package from being complete, for the closing lines
  notes = [f"{len(other_drive)} file(s) on another drive left out"] if other_drive else []

  if to_zip:
    # creating zip file
    if out_path is None:
      out_path = master_pfs.parent / (master_pfs.stem + ".zip")
    out_path = Path(out_path)
    if out_path.is_dir():
      raise UserError(f'The output path is a directory. Use -d to create a staging directory:\n\t"{out_path}"')
    if not out_path.suffix:
      out_path = out_path.with_suffix(".zip")
    if not confirm_incomplete():
      return False
    if out_path.exists() and not confirm_overwrite(f"{out_path} exists - overwrite?", force):
      print("Skip writing output file")
      return False
    with zipfile.ZipFile(out_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
      for file, relative_path in pack_entries(all_paths, real_root, place):
        zipf.write(file, relative_path)
    print(f"Zip file created: {out_path}")
  else:
    if out_path is None:
      out_path = master_pfs.parent / (master_pfs.stem + "_staging")
    out_path = Path(out_path)
    if out_path.is_file():
      raise UserError(f'The output path is a file. Omit -d/-p to create a zip file:\n\t"{out_path}"')
    if not confirm_incomplete():
      return False
    if out_path.exists():
      print(f"Writing into the existing directory: {out_path}")
      print("Files are compared by size and modification time, like robocopy does. Nothing in the directory is deleted.")
    if not fill_staging_dir(pack_entries(all_paths, real_root, place), out_path, force, notes):
      return False
    print(f"Staging directory written: {out_path}")

  for note in notes:
    print(f"  {note}")
  return True


def parse_args():
  parser = argparse.ArgumentParser(description="Pack a PFS file and its dependencies.") 
  parser.add_argument("pfs_path", nargs='?', type=str, help="Path to the main PFS file")
  parser.add_argument("-o", "--output", type=str, help="Path to the output zip file or directory")
  mode = parser.add_mutually_exclusive_group()
  mode.add_argument("-m", "--minimum", action="store_true", help="Include only currently used files (default)")
  mode.add_argument("-a", "--all",     action="store_true", help="Include all files, even if not currently used")
  parser.add_argument("-f", "--force",   action="store_true", help="Replace an existing zip file, or update differing files in the output directory, without asking")
  parser.add_argument("-i", "--incomplete", action="store_true", help="Pack without asking even if referenced files are on a different drive than the main pfs file; they are left out")
  layout = parser.add_mutually_exclusive_group()
  layout.add_argument("-d", "--dir",   action="store_true", help="Create a directory as output instead of a zip file, with the same layout as the zip file")
  layout.add_argument("-p", "--place", action="store_true", help="Like -d, but the main pfs file lands directly in the output directory (-o required); "
                                                                 "referenced files keep their position relative to it, ..\\ references land above the directory")

  args = parser.parse_args()
  if args.place and not args.output:
    parser.error("-p/--place requires -o/--output")

  # Interactive use (double click, Explorer command, console without arguments) prompts for a missing path
  # and pauses at the end, so that the console window stays open. The Explorer command passes the path, so
  # the absence of -m/-a is what tells it from a complete command line, which does not pause.
  args.interactive = not args.pfs_path or not (args.minimum or args.all)

  # Prompt for path if not provided
  if not args.pfs_path:
    try:
      import tkinter as tk
      from tkinter import filedialog
      dialog = tk.Tk()
      dialog.withdraw()
    except Exception as e: # tkinter not installed, or no display
      print(f"Cannot open file selection dialog ({e}). Please enter the full path to the main PFS file:")
      path = input("Path: ").strip()
      if not path:
        raise UserAbort("No path provided.")
    else:
      path = filedialog.askopenfilename(title="Select main PFS file", filetypes=[("MIKE SHE setup files", "*.she"), ("MIKE SHE water balance files", "*.wbl"), ("All files", "*.*")])
      dialog.destroy()
      if not path:
        raise UserAbort("No file selected.")
    args.pfs_path = path

  return args


if __name__ == "__main__":
  interactive = True
  exit_code = 1
  try:
    args = parse_args()
    interactive = args.interactive
    to_zip = not (args.dir or args.place)
    if main(Path(args.pfs_path), args.all, to_zip, args.place, args.output, args.force, args.incomplete):
      exit_code = 0
  except UserAbort as e: # the user's own decision, not an error: no traceback, but the pause below still applies
    print(f"{e} Exiting.")
  except UserError as e:
    print(f"Error: {e}")
  except Exception as e:
    traceback.print_exc()
  if interactive:
    try:
      input("\nPress Enter to exit...")
    except EOFError:
      pass
  sys.exit(exit_code)
