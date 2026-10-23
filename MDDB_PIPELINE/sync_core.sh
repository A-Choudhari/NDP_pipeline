#!/usr/bin/env bash
# sync_core.sh -- keep this pipeline's core/ in step with the NDP pipeline's.
#
#   ./sync_core.sh            show what differs
#   ./sync_core.sh --apply    copy the NDP pipeline's modules over ours
#
# `mwf.py` and `workdir.py` belong to this pipeline only and are never touched.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SRC="$(cd "$HERE/../NDP_PIPELINE/core" && pwd -P)" || { echo "no NDP_PIPELINE beside this one" >&2; exit 2; }
OURS_ONLY=(mwf.py workdir.py)

apply=0; [[ "${1:-}" == "--apply" ]] && apply=1
changed=0
for f in "$SRC"/*.py; do
    name="$(basename "$f")"
    skip=0
    for own in "${OURS_ONLY[@]}"; do [[ "$name" == "$own" ]] && skip=1; done
    [[ $skip == 1 ]] && continue
    if ! cmp -s "$f" "$HERE/core/$name"; then
        changed=$((changed + 1))
        echo "--- $name"
        diff -u "$HERE/core/$name" "$f" | head -40
        [[ $apply == 1 ]] && cp "$f" "$HERE/core/$name" && echo "    copied"
    fi
done
if [[ $changed == 0 ]]; then
    echo "core/ is in step with $SRC"
elif [[ $apply == 0 ]]; then
    echo
    echo "$changed file(s) differ. Copy them with:  $(basename "$0") --apply"
fi
