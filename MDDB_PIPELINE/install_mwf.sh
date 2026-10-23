#!/usr/bin/env bash
# install_mwf.sh -- install the MDDB workflow (mwf) in your own space, no sudo.
#
#   install_mwf.sh                 # installs into ~/Software
#   install_mwf.sh /some/dir       # or wherever you prefer
#
# Creates:
#   <dir>/mwf        a clone of https://github.com/mmb-irb/MDDB-workflow
#   conda env        mwf_env  (ambertools, GROMACS, VMD, rdkit, ...)
#
# Needs conda and internet access. The environment is several GB and the first
# solve can take 20-60 minutes.

set -euo pipefail

DEST="${1:-$HOME/Software}"
mkdir -p "$DEST"
DEST="$(cd "$DEST" && pwd -P)"

# Never inside the lab data tree (lab.yaml `storage: data_roots:`): a new
# top-level directory there is picked up by the SUBSTRATE index as
# unexplained content.
# If Python cannot read lab.yaml (no PyYAML, a broken file), the roots are
# taken from its `data_roots: [...]` line by hand; failing open would let an
# install land in the data tree.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
LAB_FILE="${MD_PIPELINE_LAB:-$HERE/../lab.yaml}"
if ! ROOTS="$(cd "$HERE" && "${SUBSTRATE_PYTHON:-python3}" -c \
        'from core import lab; print("\n".join(lab.data_roots()))' 2>/dev/null)"; then
    ROOTS="$(sed -n 's/^[[:space:]]*data_roots:[[:space:]]*\[\(.*\)\].*/\1/p' "$LAB_FILE" 2>/dev/null \
             | tr ',' '\n' | sed "s/^[[:space:]'\"]*//; s/[[:space:]'\"]*\$//")"
fi
while IFS= read -r root; do
    [[ -n "${root%/}" ]] || continue          # `/` alone would refuse everything
    case "$DEST/" in
        "${root%/}"/*)
            echo "Please install outside $root (e.g. your home: ~/Software)." >&2
            exit 2 ;;
    esac
done <<< "$ROOTS"

# --- find conda ---------------------------------------------------------------
CONDA_SH=""
for c in "${CONDA_EXE:-}" "$(command -v conda 2>/dev/null || true)" \
         "$HOME/Software/anaconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
         "$HOME/miniconda3/bin/conda" "$HOME/miniforge3/bin/conda"; do
    [[ -n "$c" && -x "$c" ]] || continue
    CONDA_SH="$(dirname "$(dirname "$c")")/etc/profile.d/conda.sh"
    [[ -f "$CONDA_SH" ]] && break
    CONDA_SH=""
done
if [[ -z "$CONDA_SH" ]]; then
    echo "conda not found. Install Miniforge in your home directory first:" >&2
    echo "  https://github.com/conda-forge/miniforge" >&2
    exit 2
fi
# shellcheck disable=SC1090
source "$CONDA_SH"

# --- clone --------------------------------------------------------------------
if [[ -d "$DEST/mwf/.git" ]]; then
    echo "repository already present: $DEST/mwf"
else
    git clone --depth 1 https://github.com/mmb-irb/MDDB-workflow "$DEST/mwf"
fi

# --- environment ----------------------------------------------------------------
# Captured first: `grep -q` quitting early can SIGPIPE the producer, and under
# pipefail that turns "found" into "not found".
ENVS="$(conda env list | awk '{print $1}')"
if grep -qx mwf_env <<<"$ENVS"; then
    echo "conda environment mwf_env already exists"
else
    echo "creating conda environment mwf_env (this is the slow part) ..."
    conda env create --file "$DEST/mwf/envs/environment.yml" -y
fi

# Activation runs every package's own activation script, and GROMACS's GMXRC
# references variables that are never set (`shell`, `GMXLDLIB`). Under `set -u`
# that is fatal and kills this script outright, so relax it for activation only.
set +u
conda activate mwf_env
set -u
python -m pip install --quiet -e "$DEST/mwf"

echo
if mwf --help >/dev/null 2>&1; then
    echo "mwf installed and working."
    echo
    echo "  activate:  source $CONDA_SH && conda activate mwf_env"
    echo "  validate:  $(dirname "$0")/validate_with_mwf.sh /path/to/<folder>"
else
    echo "mwf was installed but does not start. Run 'mwf --help' in mwf_env to see why."
    exit 1
fi
