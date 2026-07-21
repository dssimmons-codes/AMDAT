#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

: "${CONDA_PREFIX:?Run this test inside the AMDAT conda/micromamba environment}"

cxx="${CXX:-x86_64-conda-linux-gnu-g++}"
cc="${CC:-x86_64-conda-linux-gnu-gcc}"
sanitizers=(-fsanitize=address,leak,undefined -fno-omit-frame-pointer)
cxxflags=(-std=c++17 -O1 -g "${sanitizers[@]}" -fopenmp)
sanitizer_string="${sanitizers[*]}"

make clean
make -j"${JOBS:-2}" \
  CC="$cc" \
  CXX="$cxx $sanitizer_string" \
  WARN="$sanitizer_string" \
  MODE=debug \
  OMP=1 \
  FFTW_ROOT="$CONDA_PREFIX"

mapfile -t objects < <(
  find build -type f -name '*.o' ! -path 'build/amdat.o' -print | sort
)

"$cxx" \
  -Isrc \
  "${cxxflags[@]}" \
  tests/padded_analysis_lifecycle.cpp \
  "${objects[@]}" \
  "${sanitizers[@]}" \
  -Lthird_party/voro++-0.4.6/src \
  -L"$CONDA_PREFIX/lib" \
  -lvoro++ \
  -lfftw3 \
  -fopenmp \
  -lm \
  -o build/padded_analysis_lifecycle_test

ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1 \
  build/padded_analysis_lifecycle_test
