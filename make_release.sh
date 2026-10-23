#!/usr/bin/env bash
# make_release.sh -- build the archive to share with other labs.
#
#   ./make_release.sh [OUTDIR]        # default: the folder above the kit
#
# The archive is this kit with its own lab.yaml replaced by the empty
# lab.template.yaml, so a lab that unzips it starts from nothing of yours.
# Your lab.yaml stays as it is here. Caches and editor folders are left out.
# The version is the first `## vX.Y (date)` heading of CHANGELOG.md.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
OUT="$(cd "${1:-$HERE/..}" && pwd -P)"
NAME=MD_NDP_PIPELINE
HEAD="$(grep -m1 -E '^## v[0-9]' "$HERE/CHANGELOG.md")"
VERSION="$(sed -E 's/^## (v[0-9.]+).*/\1/' <<< "$HEAD")"
DATE="$(sed -E 's/.*\(([0-9-]+)\).*/\1/' <<< "$HEAD")"
ZIP="$OUT/${NAME}_${VERSION}_${DATE}.zip"

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
rsync -a --exclude __pycache__ --exclude .pytest_cache --exclude .claude \
      --exclude '*.pyc' --exclude .DS_Store "$HERE/" "$STAGE/$NAME/"
cp "$HERE/lab.template.yaml" "$STAGE/$NAME/lab.yaml"
rm -f "$ZIP"
(cd "$STAGE" && zip -qr -X "$ZIP" "$NAME")
echo "$ZIP  ($(unzip -l "$ZIP" | tail -1 | awk '{print $2}') files)"
