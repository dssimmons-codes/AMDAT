# Making AMDAT

These instructions build AMDAT with the same toolchain on any Linux system using Conda. FFTW and the compiler toolchain come from Conda; the XDR library is vendored in the source tree; Voro++ is built as part of the normal build.

## 1) Prerequisites

- **Conda / Miniconda**
- **Linux** (for other OSes, use a container or adapt the toolchain)

### For CIRCE Users

You can load conda using the following commands.
```bash
module load apps/miniconda/4.7.12
```
 This is a bit of an ancient conda and you should ideally install your own.

## 2) Create the build environment

Navigate to the AMDAT root directory and create a tailored conda environment.
```bash
cd /path/to/AMDAT/
make conda-setup
```
This creates a conda environment called amdat to consistently make AMDAT and all of its libraries. This task has to only be done once.

## 3) Build

Activate the conda environment and build. Use the provided Make targets; they automatically use the Conda compilers and ensure Voro++ is built correctly.
```bash
conda activate amdat
make
```
Or to compile with 10 cores
```bash
make -j10
```

This produces the AMDAT executable in the repo root.

## 4) Run

```bash
./AMDAT -i path/to/input.in
```
Depending on your build, you may need to activate the conda environment before running, or:
```bash
conda run -n amdat ./AMDAT -i path/to/input.in
```
