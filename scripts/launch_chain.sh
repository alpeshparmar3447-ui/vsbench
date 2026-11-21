#!/usr/bin/env bash
# Run the whole q1 pipeline in order, detached and at low priority. Each stage resumes if re-run.
# Order matters: `perf` runs right after the screen so its timings are taken on an otherwise idle GPU.
# Usage: scripts/launch_chain.sh configs/q1.json
set -uo pipefail
cd "$(dirname "$0")/.."
cfg="${1:?usage: launch_chain.sh CONFIG.json}"
name="$(basename "$cfg" .json)"
mkdir -p "results/$name"
log="results/$name/chain.log"

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

# Re-exec once into its own session. A chain started from an interactive shell would otherwise
# share that shell's session, where stopping an unrelated job can signal this one as well.
if [ "${VSBENCH_CHAIN_DETACHED:-}" != "1" ]; then
    VSBENCH_CHAIN_DETACHED=1 setsid nice -n 10 "$0" "$@" < /dev/null > /dev/null 2>&1 &
    echo "started chain pid $! ; follow with: tail -f $log"
    exit 0
fi

stage() {
    echo "$(date '+%F %T')  START $*" >> "$log"
    "${cmd[@]}" "$@" >> "results/$name/chain_$1.out" 2>&1
    echo "$(date '+%F %T')  END   $1 exit $?" >> "$log"
}

stage run --config "$cfg"
stage perf --config "$cfg"
stage ablation --config "$cfg"
stage grid --config "$cfg"
stage report --config "$cfg"
echo "$(date '+%F %T')  CHAIN DONE" >> "$log"
