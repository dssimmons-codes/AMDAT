# AMDAT — Amorphous Molecular Dynamics Analysis Toolkit

AMDAT is a C++ toolkit for analyzing molecular dynamics trajectories, with a focus on amorphous and glassy materials. It loads a trajectory into memory for very fast analysis. It reads LAMMPS formats and (optionally) GROMACS `.xtc`.

## Quick start (Linux, Conda)

### Build with conda environment (Recommended)
```bash
make conda-setup
conda activate amdat
make
```

## Run
```bash
./AMDAT -i path/to/input.in
```

## Documentation
* [Overview](Manual/overview.md)
* [Making AMDAT](Manual/making_amdat.md)
