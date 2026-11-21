#!/usr/bin/env bash
# Launch a benchmark run detached from this shell, at low priority, resumable.
# Usage: scripts/launch.sh configs/standard.json
set -euo pipefail
cd "$(dirname "$0")/.."
cfg="${1:?usage: launch.sh CONFIG.json}"
name="$(basename "$cfg" .json)"
mkdir -p "results/$name"

# Prefer the installed console script; fall back to a project virtualenv, then to the
# source tree, so this works in a pip install, a checkout with .venv, and a bare checkout.
if command -v vsbench > /dev/null 2>&1; then
    cmd=(vsbench)
elif [ -x .venv/bin/vsbench ]; then
    cmd=(.venv/bin/vsbench)
elif [ -x .venv/bin/python ]; then
    cmd=(env PYTHONPATH=src .venv/bin/python -m vsbench.cli)
else
    cmd=(env PYTHONPATH=src python3 -m vsbench.cli)
fi

# Fail now rather than 30 s into a detached run that only logs to a file.
"${cmd[@]}" doctor > "results/$name/doctor.log" 2>&1 || {
    echo "environment check failed; see results/$name/doctor.log" >&2; exit 1; }

# setsid: own session, so stopping an unrelated job cannot signal this one (and vice versa).
# nice: the laptop stays usable. Re-running the same command resumes where it stopped.
setsid nice -n 10 "${cmd[@]}" run --config "$cfg" \
    > "results/$name/launch.log" 2>&1 < /dev/null &
echo "started pid $! ; follow with: tail -f results/$name/run.log"
