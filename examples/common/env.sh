# Sourced by every examples/*/run.sh.
set -euo pipefail
# Pin threads: unpinned BLAS thrashes on many-core machines, and GROMACS refuses
# to start when OMP_NUM_THREADS and -ntomp disagree.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MMFF_THREADS="${MMFF_THREADS:-$OMP_NUM_THREADS}" POLYLEGO_MMFF_THREADS="${POLYLEGO_MMFF_THREADS:-$OMP_NUM_THREADS}"
export PYTHONHASHSEED="${PYTHONHASHSEED:-0}"
PY="${PYTHON:-python}"
GMX="${POLYLEGO_GMX:-$(command -v gmx || command -v gmx_mpi || true)}"
if [ -z "$GMX" ]; then
    echo "GROMACS not found: put gmx (or gmx_mpi) on PATH, or set POLYLEGO_GMX." >&2
    exit 2
fi
export POLYLEGO_GMX="$GMX"
# GPU for the GROMACS runs: auto (use a compatible GPU if there is one, else CPU),
# none (CPU only) or a GPU id, e.g.  GPU=0 bash run.sh
GPU="${GPU:-${POLYLEGO_GPU:-auto}}"
EXAMPLES="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="$EXAMPLES/common${PYTHONPATH:+:$PYTHONPATH}"
step() { echo; echo "=== $* ==="; }
