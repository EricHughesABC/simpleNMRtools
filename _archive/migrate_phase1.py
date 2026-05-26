"""
simpleNMRtools — Phase 1 migration script
==========================================
Reorganises the flat directory structure into:
  core/      NMR processing modules
  utils/     helper modules
  config/    NMR constants (globals, column headers)
  notebooks/ Jupyter notebooks
  archive/   backup directories and old files

SAFE TO RUN:  originals are kept in place until Phase 2 (import updates).
IDEMPOTENT:   can be re-run without side effects.

Prerequisites
-------------
  git add -A && git commit -m "pre-migration snapshot"

Run from the project root:
  python migrate_phase1.py
"""

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).parent.resolve()


# ── helpers ───────────────────────────────────────────────────────────────────

def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def copy(src: Path, dst: Path) -> None:
    if not src.exists():
        print(f"  SKIP  (not found) {src.relative_to(ROOT)}")
        return
    if dst.exists():
        print(f"  SKIP  (already exists) {dst.relative_to(ROOT)}")
        return
    shutil.copy2(src, dst)
    print(f"  COPY  {src.relative_to(ROOT)}  →  {dst.relative_to(ROOT)}")


def move_to_archive(src: Path, archive: Path) -> None:
    dst = archive / src.name
    if not src.exists():
        print(f"  SKIP  (not found) {src.relative_to(ROOT)}")
        return
    if dst.exists():
        print(f"  SKIP  (already archived) {src.relative_to(ROOT)}")
        return
    shutil.move(str(src), str(dst))
    print(f"  ARCHIVE  {src.relative_to(ROOT)}  →  archive/{src.name}")


def write_if_missing(path: Path, content: str = "") -> None:
    if path.exists():
        print(f"  SKIP  (already exists) {path.relative_to(ROOT)}")
        return
    path.write_text(content, encoding="utf-8")
    print(f"  CREATE  {path.relative_to(ROOT)}")


# ── Step 1: create new package directories ───────────────────────────────────

def create_directories() -> None:
    print("\n── Step 1: Create directories ──")
    for d in ["core", "utils", "config", "notebooks", "archive"]:
        mkdir(ROOT / d)
        print(f"  DIR  {d}/")


# ── Step 2: archive backup directories and old files ─────────────────────────

ARCHIVE_DIRS = [
    "bckup_files",
    "deprecated",
    "MNOVAscripts_1",
    "MNOVAscripts_2",
    "MNOVAscripts_3",
    "MNOVAscripts_17",
    "MNOVA_scripts_bckup_19_29",
    "MNOVAscripts_3Nov2025_bckup",
]

ARCHIVE_FILES = [
    "simulatedAnnealing_v5.py",
    "simpleNMRtest_app_bckup_21Nov2025.py",
]


def archive_old_files() -> None:
    print("\n── Step 2: Archive old files and directories ──")
    archive = ROOT / "archive"
    for name in ARCHIVE_DIRS:
        move_to_archive(ROOT / name, archive)
    for name in ARCHIVE_FILES:
        move_to_archive(ROOT / name, archive)


# ── Step 3: move notebooks ───────────────────────────────────────────────────

def move_notebooks() -> None:
    print("\n── Step 3: Move Jupyter notebooks ──")
    nb_dir = ROOT / "notebooks"
    for nb in ROOT.glob("*.ipynb"):
        dst = nb_dir / nb.name
        if dst.exists():
            print(f"  SKIP  (already exists) {nb.name}")
        else:
            shutil.move(str(nb), str(dst))
            print(f"  MOVE  {nb.name}  →  notebooks/{nb.name}")


# ── Step 4: copy core NMR modules ────────────────────────────────────────────

CORE_FILES = {
    "nmrsolution.py":            "nmrsolution.py",
    "expectedmolecule.py":       "expectedmolecule.py",
    "html_from_assignments.py":  "html_from_assignments.py",
    "simulatedAnnealing_v5a.py": "simulated_annealing.py",
}


def copy_core_modules() -> None:
    print("\n── Step 4: Copy core NMR modules to core/ ──")
    dst_dir = ROOT / "core"
    write_if_missing(dst_dir / "__init__.py",
        '"""NMR processing core."""\n')
    for src_name, dst_name in CORE_FILES.items():
        copy(ROOT / src_name, dst_dir / dst_name)


# ── Step 5: copy utility modules ─────────────────────────────────────────────

UTILS_FILES = {
    "javaUtils.py":        "java_utils.py",
    "jsonUtils.py":        "json_utils.py",
    "mnovautils.py":       "mnova_utils.py",
    "registrationUtils.py": "registration_utils.py",
}


def copy_utils_modules() -> None:
    print("\n── Step 5: Copy utility modules to utils/ ──")
    dst_dir = ROOT / "utils"
    write_if_missing(dst_dir / "__init__.py",
        '"""Utility helpers."""\n')
    for src_name, dst_name in UTILS_FILES.items():
        copy(ROOT / src_name, dst_dir / dst_name)


# ── Step 6: copy config/constants modules ────────────────────────────────────

CONFIG_FILES = {
    "globals.py":      "globals.py",
    "excelheaders.py": "excel_headers.py",
}


def copy_config_modules() -> None:
    print("\n── Step 6: Copy NMR constants to config/ ──")
    dst_dir = ROOT / "config"
    write_if_missing(dst_dir / "__init__.py",
        '"""NMR configuration constants."""\n')
    for src_name, dst_name in CONFIG_FILES.items():
        copy(ROOT / src_name, dst_dir / dst_name)


# ── Step 7: fix app/ package ─────────────────────────────────────────────────

def fix_app_package() -> None:
    print("\n── Step 7: Fix app/ package ──")
    app_dir = ROOT / "app"

    init_py = app_dir / "__init__.py"
    init_wrongname = app_dir / "init.py"

    if init_py.read_text(encoding="utf-8").strip() == "" and init_wrongname.exists():
        # Copy content of init.py into __init__.py, rename init.py to init.py.bak
        content = init_wrongname.read_text(encoding="utf-8")
        init_py.write_text(content, encoding="utf-8")
        init_wrongname.rename(app_dir / "init.py.bak")
        print("  FIX   app/init.py  →  app/__init__.py  (content moved)")
    else:
        print("  SKIP  app/__init__.py (already has content)")

    # Fix models/__init__.py similarly
    models_dir = app_dir / "models"
    if models_dir.exists():
        models_init = models_dir / "__init__.py"
        models_wrongname = models_dir / "init.py"
        if models_init.exists() and models_wrongname.exists():
            existing = models_init.read_text(encoding="utf-8").strip()
            if existing == "":
                wrong_content = models_wrongname.read_text(encoding="utf-8")
                models_init.write_text(wrong_content, encoding="utf-8")
                models_wrongname.rename(models_dir / "init.py.bak")
                print("  FIX   app/models/init.py  →  app/models/__init__.py")
            else:
                print("  SKIP  app/models/__init__.py (already has content)")


# ── Step 8: move stray template backup ───────────────────────────────────────

def tidy_templates() -> None:
    print("\n── Step 8: Tidy templates/ ──")
    templates = ROOT / "templates"
    backups = templates / "backups"
    mkdir(backups)

    stray = templates / "d3molplotmnova_template_bckup_31May2025.html"
    if stray.exists():
        shutil.move(str(stray), str(backups / stray.name))
        print(f"  MOVE  templates/{stray.name}  →  templates/backups/")
    else:
        print("  SKIP  no stray template backups found")


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    print(f"simpleNMRtools Phase 1 migration")
    print(f"Root: {ROOT}\n")

    # Safety check
    if not (ROOT / "simpleNMRtest_app.py").exists():
        print("ERROR: simpleNMRtest_app.py not found.")
        print("       Run this script from the project root directory.")
        sys.exit(1)

    create_directories()
    archive_old_files()
    move_notebooks()
    copy_core_modules()
    copy_utils_modules()
    copy_config_modules()
    fix_app_package()
    tidy_templates()

    print("\n── Migration complete ──")
    print("""
Next steps
----------
1. Run the tests / start the server — original flat files are still in place,
   so simpleNMRtest_app.py should work unchanged.

2. Phase 2: update import paths in the new core/, utils/, config/ files
   (Claude will produce updated files in the next session).

3. Phase 3: extract routes from simpleNMRtest_app.py into app/routes.py,
   complete the app factory, update the PythonAnywhere WSGI file.
""")


if __name__ == "__main__":
    main()
