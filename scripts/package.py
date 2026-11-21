"""Copy the manuscript, protocol and a clean code archive to a delivery folder.
Usage: python scripts/package.py DEST_DIR"""
import shutil, subprocess, sys, zipfile
from pathlib import Path

root = Path(__file__).resolve().parents[1]
dest = Path(sys.argv[1]); dest.mkdir(parents=True, exist_ok=True)
# Rebuild the Word versions from the Markdown so the two cannot drift apart. The submission
# document is rebuilt here too: it embeds the figures and tables from the result files, so a
# stale copy is the one most likely to disagree with what the analysis currently says.
for script in ("make_docx.py", "make_submission_docx.py"):
    try:
        subprocess.run([sys.executable, str(root / "scripts" / script)],
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        print(f"note: {script} failed, shipping whatever exists ({exc.stderr.strip()[:200]})")
    except Exception as exc:                   # pandoc absent: ship the Markdown alone
        print(f"note: {script} skipped ({exc})")

copied = []
for name in ("paper/manuscript_draft.md", "paper/manuscript_draft.docx",
             "paper/manuscript_submission.docx",
             "paper/manuscript_additional_file_1.docx", "paper/manuscript_tables.docx",
             "PROTOCOL.md", "README.md",
             "LICENSE", "CITATION.cff", "requirements-lock.txt"):
    f = root / name
    if f.exists():
        shutil.copy2(f, dest / f.name)
        copied.append(f.name)
    else:
        print(f"note: {name} not found, not delivered")
for extra in sorted((root / "results").glob("*/tables*.md")):
    shutil.copy2(extra, dest / f"{extra.stem}_{extra.parent.name}.md")
# figures the manuscript's legends refer to, one folder per run
for figdir in sorted((root / "results").glob("*/figures")):
    out = dest / f"figures_{figdir.parent.name}"
    out.mkdir(exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        for img in sorted(figdir.glob(f"*.{ext}")):
            shutil.copy2(img, out / img.name)
tables = dest / "template_tables"
tables.mkdir(exist_ok=True)
for name in ("metal_audit", "structure_qc", "sites", "site_controls"):
    f = root / "results" / f"{name}.csv"
    if f.exists():
        shutil.copy2(f, tables / f.name)
skip = {"__pycache__", ".pytest_cache"}
with zipfile.ZipFile(dest / "vsbench_code.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for top in ("src", "tests", "configs", "scripts"):
        for p in sorted((root / top).rglob("*")):
            if p.is_file() and not (set(p.parts) & skip) and p.suffix != ".pyc":
                z.write(p, p.relative_to(root))
    for f in ("pytest.ini", "README.md", "PROTOCOL.md", "LICENSE",
              "CITATION.cff", "pyproject.toml", "requirements-lock.txt"):
        if (root / f).exists():
            z.write(root / f, f)
# Report what this run copied, not what the folder happens to contain: a listing would show a
# stale file from an earlier delivery as though it had just been written.
print(f"delivered to {dest}")
for n in copied:
    print("  ", n)
print("   figures, template tables and vsbench_code.zip")
