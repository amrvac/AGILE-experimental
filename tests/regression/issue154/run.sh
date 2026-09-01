#!/bin/sh
# Build and run tests/regression/issue154/repro154.f90 with each toolchain.
# Each run reports 4 variants: 1 = the buggy collapse(3)+call form (runtime
# nw), 2 & 3 = the documented workarounds, 4 = variant 1 with nw a compile-
# time parameter.  Expected:
#   gfortran -fopenacc      all 4 PASS
#   nvfortran, no -acc      all 4 PASS
#   nvfortran -acc=gpu      variant 1 FAIL, variants 2-4 PASS  (issue #154)
#                           (flip radius_dependent_flow to also fail variant 2)
set -e
cd "$(dirname "$0")"

GF=${GF:-gfortran}
NVF=${NVF:-nvfortran}      # e.g. /opt/nvidia/hpc_sdk/Linux_x86_64/24.11/compilers/bin/nvfortran

echo "== gfortran -fopenacc (host) =="
$GF -O2 -fopenacc repro154.f90 -o repro154.gnu   && ./repro154.gnu   || true

echo "== nvfortran, no -acc (host) =="
$NVF -O3 -fast repro154.f90 -o repro154.nvhost   && ./repro154.nvhost || true

echo "== nvfortran -acc=gpu (AGILE arch/nvidia.mk flags) =="
$NVF -cpp -Mfree -acc=gpu -Mvect=levels:5 -Minline -g -O3 -fast \
     repro154.f90 -o repro154.nvgpu && ./repro154.nvgpu || true
