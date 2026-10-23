#!/usr/bin/env bash
# collect_pipeline.sh -- gather finished MD_RESULTS folders into one place.
#
#   collect_pipeline.sh [SOURCE ...] [--all] [--dest DIR] [--dry-run]
#
# Each SOURCE is a project's MD_RESULTS/ directory. Each is copied whole, byte
# for byte, to <dest>/<PREFIX>_MD_RESULTS/. Default dest:
#   $MD_COLLECT_DEST, else ~/MD_NDP_DATA
#
#   --all                   find every MD_RESULTS under the roots in sources.yaml
#   --dry-run               say what would happen; write nothing at all
#   --force                 re-copy sources that have not changed
#   --verify-hashes         re-read the share and check it against the manifests
#   --no-production-files   omit the production_files/ symlinks
#
# This is a SEPARATE pipeline. It builds no records and analyses nothing: it
# reads finished output and puts it where the lab can see it. NDP_PIPELINE and
# MD_ANALYSES are untouched by it, and it never writes into a source.
set -uo pipefail

KIT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

PY="${SUBSTRATE_PYTHON:-python3}"
if ! "$PY" -c "import json, hashlib" 2>/dev/null; then
    echo "  $PY cannot run this. Name one that can:" >&2
    echo "      export SUBSTRATE_PYTHON=/path/to/python3" >&2
    exit 2
fi

# Refuse to run shared code that someone other than its owner could have
# changed. Same guard, same reason, as the other two kits.
if [[ -f "$KIT/trust_check.sh" ]]; then
    # shellcheck source=/dev/null
    source "$KIT/trust_check.sh"
    trust_check "$KIT" || exit 2
fi

exec "$PY" "$KIT/collect.py" "$@"
