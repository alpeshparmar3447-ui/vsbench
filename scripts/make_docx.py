"""Convert the Markdown manuscript and its tables to Word documents.

Two files, matching how journals take a submission: the manuscript with its figure legends,
and the tables separately in landscape because the widest result table has sixteen columns.
Styling comes from the reference documents in ``paper/style``, which set a serif face at 12 pt,
1.5 line spacing and black headings rather than pandoc's default sans-serif and blue.

Usage: python scripts/make_docx.py [OUT_DIR]      (default: paper/)
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STYLE = ROOT / "paper" / "style"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "paper"


def need(path: Path, what: str) -> Path:
    if not path.exists():
        raise SystemExit(f"missing {what}: {path}")
    return path


def main() -> int:
    if shutil.which("pandoc") is None:
        raise SystemExit("pandoc is required: https://pandoc.org/installing.html")
    OUT.mkdir(parents=True, exist_ok=True)

    manuscript = need(ROOT / "paper" / "manuscript_draft.md", "manuscript")
    ref_ms = need(STYLE / "reference_manuscript.docx", "manuscript style reference")
    subprocess.run(
        ["pandoc", str(manuscript), "--from=markdown+superscript+subscript",
         f"--reference-doc={ref_ms}", "-o", str(OUT / "manuscript_draft.docx")],
        check=True)
    print(f"wrote {OUT / 'manuscript_draft.docx'}")

    # Tables live in the run's own output, so they are current by construction. Which run,
    # though, has to be decided rather than discovered: more than one arm can carry a report,
    # and taking whichever the filesystem lists first would ship one arm's numbers unannounced.
    arms = sorted(p.parent for p in (ROOT / "results").glob("*/tables.md"))
    named = os.environ.get("VSBENCH_RESULTS")
    if named:
        arms = [a for a in arms if a.name == named or str(a) == named]
        if not arms:
            raise SystemExit(f"no report for results arm {named!r}")
    if not arms:
        print("no results/*/tables.md found; skipped the tables document")
        return 0
    if len(arms) > 1:
        raise SystemExit(f"several arms carry a report ({', '.join(a.name for a in arms)}); "
                         f"name one with VSBENCH_RESULTS=<arm>")
    tables_md = arms[0] / "tables.md"
    ref_tbl = need(STYLE / "reference_tables.docx", "tables style reference")
    title = ("# Tables\n\nSupporting tables for *Docking does not outrank a ligand-similarity "
             "baseline on property-matched experimentally confirmed inactives*.\n\n")
    src = OUT / "_tables_src.md"
    src.write_text(title + tables_md.read_text())
    try:
        subprocess.run(
            ["pandoc", str(src), "--from=markdown", f"--reference-doc={ref_tbl}",
             "-o", str(OUT / "manuscript_tables.docx")], check=True)
    finally:
        src.unlink(missing_ok=True)
    print(f"wrote {OUT / 'manuscript_tables.docx'}  (from {tables_md})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
