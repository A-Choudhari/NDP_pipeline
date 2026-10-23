#!/usr/bin/env bash
# ndp_pipeline.sh -- one folder in, one JSON record per simulated system out.
#
#   ndp_pipeline.sh /path/to/<folder> [--yes] [options]
#
# Steps, in order:
#   1. set up          settings.yaml (and move a folder from the old layout)
#   2. index           read every trajectory header once
#   3. datasets        propose what to describe, into settings.yaml
#   4. PDB metadata    fill EMPTY descriptions and structure citations
#   5. glossary        one stable entry per MD, from configs, logs and headers
#   6. records         one JSON per system + catalog.json + the NDP sheet
#
# Writes ONLY into <folder>/MD_RESULTS/ and <folder>/.md_pipeline/, and reads
# nothing outside the folder you name. Your data is never written to.
#
#   --yes                   do not ask
#   --pdb-offline           use the cached RCSB metadata only
#   --skip-pdb-metadata     leave settings.yaml untouched in step 4
#   --deep                  scan .xtc/ASCII .mdcrd for frame counts (slow)
#   --production-view       also link the published production files (and only
#                           those) into MD_RESULTS/production_files/
#   --rebase-root           this folder was moved, renamed, or its settings were
#                           copied from somewhere else: point `folder: root:` at
#                           where the settings file actually is, then build
#   --include-nested        build anyway when this folder CONTAINS folders that
#                           have already been built (see below)

set -uo pipefail

KIT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

FOLDER=""; YES=""; DEEP=""; PDB=1; VIEW=0; REBASE=0; NESTED=0; PDB_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --yes|-y)            YES="--yes" ;;
        --deep)              DEEP="--deep" ;;
        --skip-pdb-metadata) PDB=0 ;;
        --pdb-offline)       PDB_ARGS+=("--offline") ;;
        --production-view)   VIEW=1 ;;
        --rebase-root)       REBASE=1 ;;
        --include-nested)    NESTED=1 ;;
        -h|--help)           sed -n '2,20p' "$0"; exit 0 ;;
        --*)                 echo "unknown option: $1" >&2; exit 2 ;;
        *)                   FOLDER="$1" ;;
    esac
    shift
done
[[ -n "$FOLDER" ]] || { sed -n '2,20p' "$0" >&2; exit 2; }
FOLDER="$(cd "$FOLDER" 2>/dev/null && pwd -P)" || { echo "no such folder" >&2; exit 2; }

PY="${SUBSTRATE_PYTHON:-python3}"
if ! "$PY" -c "import pandas, pyarrow, yaml" 2>/dev/null; then
    echo
    echo "  $PY cannot import pandas, pyarrow and pyyaml."
    echo "  Name one that can:"
    echo "      export SUBSTRATE_PYTHON=/path/to/python3"
    echo "      $(basename "$0") $FOLDER"
    echo
    exit 2
fi

# ONE FOLDER, ONE BUILD. Pointed at a PARENT of a folder that has already
# been built, this indexes the child's trajectories again and writes a second
# record for the same simulations. Nothing about the result looks wrong: the
# lab-wide build at the root of this tree holds 79 records that duplicate the
# 13 real ones, and it took a survey to notice. So a folder that CONTAINS
# other builds is refused, because the first instinct of somebody new is to
# point the tool at the top of their directory.
#
# The test is evidence, not a rule of thumb: every one of the 13 project
# folders has zero MD_RESULTS below it, and the root has 13.
if [[ $NESTED -eq 0 ]]; then
    NESTED_FOUND="$(find "$FOLDER" -mindepth 2 -maxdepth 6 -type d \
                        -name MD_RESULTS 2>/dev/null | head -5)"
    if [[ -n "$NESTED_FOUND" ]]; then
        N_TOTAL="$(find "$FOLDER" -mindepth 2 -maxdepth 6 -type d \
                       -name MD_RESULTS 2>/dev/null | wc -l)"
        echo
        echo "  Refused: this folder contains $N_TOTAL folder(s) that have"
        echo "  already been built:"
        echo
        # shellcheck disable=SC2001
        echo "$NESTED_FOUND" | sed "s|^|      |"
        [[ "$N_TOTAL" -gt 5 ]] && echo "      ... and $((N_TOTAL - 5)) more"
        echo
        echo "  Building here would read their trajectories a second time and"
        echo "  write a second record for each of the same simulations, with"
        echo "  nothing to say so. Point the pipeline at ONE folder holding one"
        echo "  set of MDs -- usually one of the folders listed above."
        echo
        echo "  If a lab-wide index over all of them is what you want, and the"
        echo "  duplication is understood, pass --include-nested."
        echo
        exit 2
    fi
fi

# THE ROOT IS CHECKED BEFORE ANY STEP RUNS. `folder: root:` used to be trusted
# forever, so a copied settings.yaml or a renamed folder made the build read one
# tree and write into another, silently. --rebase-root is the deliberate act
# that says "this folder is where it now is".
if [[ $REBASE -eq 1 ]]; then
    "$PY" - "$FOLDER" "$KIT" <<'PYREBASE'
import sys, os
sys.path.insert(0, sys.argv[2])
from core import settings
was, now, changed = settings.rebase_root(sys.argv[1])
if changed:
    print(f"\n  root rebased\n      was : {was}\n      now : {now}")
    print("  The previous root may mean the records here describe another "
          "folder's data. Check them.\n")
else:
    print(f"\n  root already points here: {now}\n")
PYREBASE
    [[ $? -eq 0 ]] || exit 2
fi

# Refuse to run shared code that someone other than its owner could have
# changed. See trust_check.sh for what is checked and why.
# shellcheck source=/dev/null
if [[ -f "$KIT/trust_check.sh" ]]; then
    source "$KIT/trust_check.sh"
    trust_check "$KIT" || exit 2
fi

LOGDIR="$FOLDER/.md_pipeline/logs"
mkdir -p "$LOGDIR" 2>/dev/null
STATUS="$LOGDIR/status.json"
START="$(date +%s)"
STEP=""

step() { STEP="$1"; printf '\n==== [%s] %s\n' "$1" "$2"; }
die()  {
    printf '\n  STOPPED at %s -- see the messages above.\n\n' "$1"
    "$PY" - "$STATUS" "$FOLDER" "$STEP" "$1" "$START" <<'PYS' 2>/dev/null || true
import json, os, sys, time
status, folder, step, why, start = sys.argv[1:6]
os.makedirs(os.path.dirname(status), exist_ok=True)
json.dump({"status": "finished", "folder": folder, "exit": 1,
           "failed_step": step, "note": why,
           "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "elapsed_s": int(time.time()) - int(start)}, open(status, "w"), indent=2)
PYS
    exit 1
}

echo
echo "  folder : $FOLDER"
echo "  python : $(command -v "$PY")"
echo "  kit    : $KIT"

step 1/6 "set up"
"$PY" - "$FOLDER" "$KIT" <<'PYM' || die "migrate"
import sys
sys.path.insert(0, sys.argv[2])
from core import migrate
migrate.migrate(sys.argv[1])
PYM
if [[ -f "$FOLDER/MD_RESULTS/settings.yaml" ]]; then
    echo "  existing MD_RESULTS/settings.yaml kept (edit it to change anything)"
else
    "$PY" "$KIT/init.py" "$FOLDER" $YES || die "init.py"
fi

step 2/6 "index (the first run on a large folder can take a while)"
"$PY" "$KIT/index_build.py" "$FOLDER" $DEEP || die "index_build.py"

step 3/6 "datasets"
if "$PY" -c "import sys,yaml;d=yaml.safe_load(open(sys.argv[1]));sys.exit(0 if (d or {}).get('datasets') else 1)" \
        "$FOLDER/MD_RESULTS/settings.yaml" 2>/dev/null; then
    echo "  datasets already declared in settings.yaml (your edits are safe)"
else
    "$PY" "$KIT/datasets_init.py" "$FOLDER" --yes
    rc=$?
    case "$rc" in
        0) ;;
        1) echo
           echo "  No production runs to describe in this folder."
           echo "  The index and $FOLDER/MD_RESULTS/README.md are complete on their own."
           echo
           exit 0 ;;
        *) die "datasets_init.py (exit $rc)" ;;
    esac
fi

step 4/6 "PDB metadata (preserving existing values)"
if [[ "$PDB" == 1 ]]; then
    "$PY" "$KIT/enrich_pdb.py" "$FOLDER" --apply "${PDB_ARGS[@]+"${PDB_ARGS[@]}"}" \
        || die "enrich_pdb.py"
else
    echo "  skipped (--skip-pdb-metadata)"
fi

step 5/6 "MD glossary"
"$PY" "$KIT/glossary_build.py" "$FOLDER" || die "glossary_build.py"

step 6/6 "records"
"$PY" "$KIT/ndp_build.py" "$FOLDER" || die "ndp_build.py"

if [[ $VIEW -eq 1 ]]; then
    "$PY" "$KIT/production_view.py" "$FOLDER" || die "production_view.py"
fi

"$PY" - "$STATUS" "$FOLDER" "$START" <<'PYS' 2>/dev/null || true
import json, os, sys, time
status, folder, start = sys.argv[1:4]
os.makedirs(os.path.dirname(status), exist_ok=True)
json.dump({"status": "finished", "folder": folder, "exit": 0,
           "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "elapsed_s": int(time.time()) - int(start)}, open(status, "w"), indent=2)
PYS

cat <<EOF

  Done in $(( $(date +%s) - START ))s. Everything is in $FOLDER/MD_RESULTS/ :

      README.md              start here: what was found, what to check
      settings.yaml          the one file to edit
      glossary.md            every MD  (glossary.parquet for scripts)
      catalog.json           every record, for the lab catalogue
      NDP_registration.md    what to paste into the NDP form
      *.json                 one record per simulated system  (production only:
                             no minimisation, heating or NPT equilibration)

  The machinery is in $FOLDER/.md_pipeline/ and never needs opening.

  Next: fill in what only you know in settings.yaml -- the identifier, the
  topic, the title, the paper, the construct, the glycans, the ions -- and
  re-run. Your text is never overwritten.
EOF
