# AMDAT — Amorphous Molecular Dynamics Analysis Toolkit

AMDAT is a C++ toolkit for post-processing molecular dynamics trajectories, with a focus on amorphous/glassy/polymer materials. It loads a trajectory into memory for very fast analysis. It reads LAMMPS formats and (optionally) GROMACS `.xtc`. It provides high-performance analyses, such as clustering, spatial decomposition, and calculating time-resolved structure factors, mean-square displacements, radial distribution functions, etc.

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

## License
GNU GPLv3.0 (see [LICENSE](LICENSE)
<!-- ## Contributing -->
<!-- We welcome issues, discussions, and PRs. Please read [CONTRIBUTING.md](CONTRIBUTING.md) and the [Code of Conduct](COD first. -->
