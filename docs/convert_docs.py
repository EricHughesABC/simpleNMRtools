#!/usr/bin/env python3
"""
convert_docs.py — batch-convert docs/source/<chapter>/<chapter>.docx to .rst

Expects the layout:
    docs/source/<chapter>/<chapter>.docx

For each one found, it:
  1. Makes a throwaway TEMP COPY of the docx — the original master file in
     docs/source/<chapter>/ is never modified.
  2. On that temp copy: flattens Word's SEQ auto-number fields (used for
     Figure/Table/Equation captions) into plain text, since pandoc silently
     drops them otherwise, producing captions like "Figure " with no number.
  3. On that temp copy: (unless --no-heading-fix) promotes oversized
     "Normal"-styled paragraphs to proper Heading styles, using font size as
     a heuristic, and prints exactly what it changed so you can sanity-check
     it — since pandoc only turns real Heading-styled paragraphs into RST
     section titles.
  4. Runs pandoc against the temp copy, writing <chapter>.rst and extracted
     media into the chapter folder itself (matching your existing
     run_pandoc.bat convention of running pandoc from inside the folder).

Usage:
    python3 convert_docs.py                     # looks under ./docs/source
    python3 convert_docs.py path/to/docs/source  # explicit source dir
    python3 convert_docs.py --dry-run            # preview only, no changes
    python3 convert_docs.py --no-heading-fix     # skip the heading heuristic

Requires: pandoc on PATH, and `pip3 install python-docx`.
"""
import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

try:
    import docx
except ImportError:
    sys.exit("Missing dependency. Install with: pip3 install python-docx")


# ---------------------------------------------------------------------------
# Step 1: flatten SEQ fields (Figure/Table/Equation auto-numbers) to text
# ---------------------------------------------------------------------------

def _flatten_seq_fields_xml(document_xml: str) -> tuple[str, int]:
    """Replace <w:fldSimple ... SEQ ...>...</w:fldSimple> blocks with plain
    text runs containing the field's last-cached value (the actual number
    Word had computed), which pandoc otherwise drops entirely."""

    def _flatten(m: "re.Match[str]") -> str:
        block = m.group(0)
        texts = re.findall(r"<w:t[^>]*>([^<]*)</w:t>", block)
        cached_text = "".join(texts)
        return f'<w:r><w:t xml:space="preserve">{cached_text}</w:t></w:r>'

    return re.subn(
        r'<w:fldSimple\b[^>]*\bw:instr="[^"]*SEQ[^"]*"[^>]*>.*?</w:fldSimple>',
        _flatten,
        document_xml,
        flags=re.DOTALL,
    )


def flatten_seq_fields(src: Path, dst: Path) -> int:
    """Read docx at src, write a copy at dst with SEQ fields flattened.
    Returns the number of fields flattened."""
    count = 0
    with zipfile.ZipFile(src, "r") as zin, zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/document.xml":
                xml = data.decode("utf-8")
                xml, count = _flatten_seq_fields_xml(xml)
                data = xml.encode("utf-8")
            zout.writestr(item, data)
    return count


# ---------------------------------------------------------------------------
# Step 2: heuristic heading promotion
# ---------------------------------------------------------------------------

def promote_headings(path: Path) -> list[str]:
    """Find 'Normal'-styled paragraphs with a uniform run font size larger
    than typical body text, and promote them to Heading 1/2/3 based on
    relative size ranking (largest size in the doc -> Heading 1, and so on).
    Modifies the docx at `path` in place. Returns a list of human-readable
    descriptions of what changed, for review."""
    d = docx.Document(str(path))
    changes: list[str] = []

    candidates: list[tuple[object, float]] = []
    for p in d.paragraphs:
        if not p.text.strip() or p.style.name != "Normal":
            continue
        sizes = {r.font.size.pt for r in p.runs if r.font.size}
        if len(sizes) == 1:
            size = sizes.pop()
            if size >= 14:  # heuristic: body text is normally 10-12pt
                candidates.append((p, size))

    if not candidates:
        return changes

    distinct_sizes = sorted({s for _, s in candidates}, reverse=True)[:3]  # cap at Heading 3
    size_to_level = {s: i + 1 for i, s in enumerate(distinct_sizes)}

    for p, size in candidates:
        level = size_to_level.get(size)
        if level is None:
            continue
        style_name = f"Heading {level}"
        p.style = d.styles[style_name]
        changes.append(f'    "{p.text.strip()[:60]}" -> {style_name} ({size:g}pt)')

    if changes:
        d.save(str(path))
    return changes


# ---------------------------------------------------------------------------
# Per-chapter pipeline
# ---------------------------------------------------------------------------

def process_chapter(docx_path: Path, dry_run: bool, fix_headings: bool) -> bool:
    chapter_dir = docx_path.parent
    stem = docx_path.stem
    print(f"\n=== {chapter_dir.name}/{docx_path.name} ===")

    if dry_run:
        print("  [dry run] would flatten SEQ fields"
              + (" and promote headings" if fix_headings else "") + ", then run pandoc")
        return True

    with tempfile.TemporaryDirectory() as tmp:
        tmp_docx = Path(tmp) / docx_path.name

        n_fields = flatten_seq_fields(docx_path, tmp_docx)
        if n_fields:
            print(f"  flattened {n_fields} auto-number field(s) (e.g. Figure/Table captions)")

        if fix_headings:
            changes = promote_headings(tmp_docx)
            if changes:
                print("  promoted headings:")
                for c in changes:
                    print(c)
            else:
                print("  no heading changes (already styled, or nothing matched the heuristic)")

        cmd = ["pandoc", "--extract-media", ".", "-o", f"{stem}.rst", str(tmp_docx)]
        result = subprocess.run(cmd, cwd=chapter_dir, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  PANDOC FAILED:\n{result.stderr}")
            return False

        print(f"  wrote {chapter_dir.name}/{stem}.rst")
        return True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source_dir", nargs="?", default="docs/source",
                     help="path containing <chapter>/<chapter>.docx folders (default: docs/source)")
    ap.add_argument("--dry-run", action="store_true", help="show what would happen, make no changes")
    ap.add_argument("--no-heading-fix", action="store_true", help="skip the heading-promotion heuristic")
    args = ap.parse_args()

    if shutil.which("pandoc") is None:
        sys.exit("pandoc not found on PATH. Install it first (e.g. `brew install pandoc`).")

    source_dir = Path(args.source_dir)
    if not source_dir.is_dir():
        sys.exit(f"Not a directory: {source_dir}")

    all_docx = sorted(source_dir.glob("*/*.docx"))
    matching = [f for f in all_docx if f.stem == f.parent.name]
    stray = [f for f in all_docx if f.stem != f.parent.name]

    if stray:
        print("Skipping files that don't match the <chapter>/<chapter>.docx convention:")
        for f in stray:
            print(f"  {f}")

    if not matching:
        print(f"No matching <chapter>/<chapter>.docx files found under {source_dir}")
        return

    print(f"Found {len(matching)} chapter docx file(s) under {source_dir}")

    failures = []
    for f in matching:
        ok = process_chapter(f, args.dry_run, fix_headings=not args.no_heading_fix)
        if not ok:
            failures.append(f)

    print()
    if failures:
        print(f"Done with {len(failures)} failure(s):")
        for f in failures:
            print(f"  {f}")
        sys.exit(1)
    else:
        print("Done — all chapters converted.")


if __name__ == "__main__":
    main()
