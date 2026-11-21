"""Build one self-contained submission document: manuscript, figures, tables.

A reviewer cannot check a result they cannot see. The plain manuscript carries figure legends
and table callouts but no artwork and no numbers, so this assembles the complete evidence
package in reading order: the manuscript, then every figure placed under its own legend, then
every table in landscape.

Figures are embedded from the PNGs the analysis wrote, so the document cannot show artwork that
disagrees with the result files.

The main-text tables go into the submission document; the supplementary tables become
Additional file 1, as a separate document, which is how the journal expects them.

Usage: python scripts/make_submission_docx.py [OUT_DIR]      (default: paper/)
       VSBENCH_RESULTS=<arm> python scripts/make_submission_docx.py   (pick a result arm)
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

#: A4 in twentieths of a point, long edge horizontal.
LANDSCAPE = ('<w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/>'
             '<w:pgMar w:top="1134" w:right="1134" w:bottom="1134" w:left="1134"/>')

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "paper"
STYLE = ROOT / "paper" / "style"

#: legend heading -> figure file, in the order the manuscript introduces them
FIGURES = [
    ("Graphical abstract", "graphical_abstract.png"),
    ("Figure 1", "redocking.png"),
    ("Figure 2", "auc_summary.png"),
    ("Figure 3", "size_control.png"),
    ("Figure 4", "negative_sets_consensus.png"),
    ("Figure 5", "prioritization.png"),
    ("Figure 6", "template_difficulty.png"),
    ("Figure 7", "site_recovery.png"),
    ("Figure 8", "grid_strategies.png"),
    ("Figure 9", "performance.png"),
]


def _results_dir() -> Path:
    """The run whose numbers go into the document, resolved explicitly.

    Several arms can sit side by side under ``results/``, so taking whichever the filesystem
    lists first would put one arm's figures beside another's tables with nothing to show it
    had happened. One arm is used without ceremony; more than one has to be named.
    """
    named = os.environ.get("VSBENCH_RESULTS")
    if named:
        d = Path(named) if Path(named).is_absolute() else ROOT / "results" / named
        if not (d / "tables.md").exists():
            raise SystemExit(f"no report in {d}: run `vsbench report` for that arm first")
        return d
    dirs = sorted(p.parent for p in (ROOT / "results").glob("*/tables.md"))
    if not dirs:
        raise SystemExit("run the report first: no results/*/tables.md found")
    if len(dirs) > 1:
        raise SystemExit(f"several arms carry a report ({', '.join(p.name for p in dirs)}); "
                         f"name one with VSBENCH_RESULTS=<arm>")
    return dirs[0]


def _render(src_md: Path, out: Path, ref: Path, resource: Path | None = None) -> None:
    cmd = ["pandoc", str(src_md), "--from=markdown+superscript+subscript",
           f"--reference-doc={ref}"]
    if resource:
        cmd += ["--resource-path", str(resource)]
    subprocess.run(cmd + ["-o", str(out)], check=True)
    src_md.unlink(missing_ok=True)


def _landscape_last_section(docx: Path) -> None:
    """Turn the document's final section landscape.

    Word stores page orientation per section, and a section's properties live at its *end*: the
    sectPr inside a paragraph closes the section containing that paragraph, while the one at the
    end of the body governs everything after the last such paragraph. The tables are therefore
    governed by the body-level sectPr, which the reference document leaves without a pgSz, so it
    has to be given one here. Child order inside sectPr is schema-enforced, hence the insertion
    after footnotePr rather than at the start.
    """
    with zipfile.ZipFile(docx) as z:
        items = {i.filename: z.read(i.filename) for i in z.infolist() if not i.is_dir()}
    xml = items["word/document.xml"].decode("utf-8")

    starts = [m.start() for m in re.finditer(r"<w:sectPr\b", xml)]
    if not starts:
        raise SystemExit("no section properties in the built document")
    head, tail = xml[:starts[-1]], xml[starts[-1]:]

    # Any page size already there is replaced, not appended to.
    tail = re.sub(r"<w:pgSz\b[^/>]*/>|<w:pgMar\b[^/>]*/>", "", tail, count=2)
    # pgSz and pgMar follow footnotePr, endnotePr and type, and precede everything else.
    anchors = [m.end() for m in re.finditer(
        r"</w:footnotePr>|</w:endnotePr>|<w:type\b[^/>]*/>", tail)]
    at = max(anchors) if anchors else re.match(r"<w:sectPr\b[^>]*>", tail).end()
    items["word/document.xml"] = (head + tail[:at] + LANDSCAPE + tail[at:]).encode("utf-8")

    with zipfile.ZipFile(docx, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in items.items():
            z.writestr(name, data)


def main() -> int:
    if shutil.which("pandoc") is None:
        raise SystemExit("pandoc is required")
    manuscript = (ROOT / "paper" / "manuscript_draft.md").read_text()
    res = _results_dir()
    figdir = res / "figures"
    main_tables, supp_tables = res / "tables_main.md", res / "tables_supplementary.md"
    if not figdir.is_dir():
        raise SystemExit(f"no figures in {res}: run the analysis first")
    for f in (main_tables, supp_tables):
        if not f.exists():
            raise SystemExit(f"{f.name} is missing from {res}; re-run `vsbench report` "
                             f"so the main and supplementary tables are written separately")

    # place each figure image directly beneath its legend
    body = manuscript
    missing = []
    for label, fname in FIGURES:
        img = figdir / fname
        if not img.exists():
            missing.append(fname)
            continue
        pattern = re.compile(rf'(^\*\*{re.escape(label)}[.,][^\n]*\n)', re.M)
        if not pattern.search(body):
            missing.append(f"{label} (no legend found)")
            continue
        # artwork above its caption, which is the usual convention
        body = pattern.sub(rf'![]({img.as_posix()})\n\n\1', body, count=1)
    if missing:
        print("warning: not embedded ->", ", ".join(missing))

    # Append every table after the figures, in a landscape section: the widest result table
    # has sixteen columns and is unreadable portrait. The section break is raw OOXML because
    # page orientation is a Word section property with no Markdown equivalent.
    landscape = (
        '```{=openxml}\n'
        '<w:p><w:pPr><w:sectPr>'
        '<w:pgSz w:w="11906" w:h="16838"/>'
        '<w:pgMar w:top="1134" w:right="1134" w:bottom="1134" w:left="1134"/>'
        '</w:sectPr></w:pPr></w:p>\n'
        '```\n')
    body += ("\n\n---\n\n" + landscape + "\n# Tables\n\n"
             "Every main-text table, in the order the Results cite them. The supplementary "
             "tables are in Additional file 1.\n\n"
             + main_tables.read_text())

    OUT.mkdir(parents=True, exist_ok=True)
    ref = STYLE / "reference_manuscript.docx"
    out = OUT / "manuscript_submission.docx"
    src = OUT / "_submission_src.md"
    src.write_text(body)
    _render(src, out, ref, figdir)
    _landscape_last_section(out)

    # Additional file 1, as its own document: the journal takes supplementary material
    # separately, and a reviewer should not have to page past it to reach the references.
    add_out = OUT / "manuscript_additional_file_1.docx"
    add_src = OUT / "_additional_src.md"
    add_src.write_text("# Additional file 1: supplementary tables\n\n"
                       "Supporting tables for *Docking does not outrank a ligand-similarity "
                       "baseline on property-matched experimentally confirmed inactives: a "
                       "two-engine benchmark on LIT-PCBA*.\n\n" + supp_tables.read_text())
    _render(add_src, add_out, STYLE / "reference_tables.docx")
    _landscape_last_section(add_out)

    def n_tables(p: Path) -> int:
        return len(re.findall(r"^\*\*Table ", p.read_text(), re.M))

    n_main, n_supp = n_tables(main_tables), n_tables(supp_tables)
    print(f"wrote {out} ({len(FIGURES) - len(missing)} figures embedded, "
          f"{n_main} main tables, landscape appendix)")
    print(f"wrote {add_out} ({n_supp} supplementary tables)")
    print(f"both built from results/{res.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
