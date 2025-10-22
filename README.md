# AMDAT — Amorphous Molecular Dynamics Analysis Toolkit

> High-performance analysis of molecular dynamics (MD) trajectories for amorphous, glass-forming, and polymer systems.

AMDAT is a C++ toolkit that loads trajectories into memory for **rapid** analysis, with robust **atom selection**, **time-resolved** statistics, and modular data objects.
It reads common LAMMPS trajectory formats and GROMACS `.xtc`.
It provides a wide variety of high-performance analyses integral to molecular modeling studies, such as clustering, spatial decomposition, and calculating time-resolved structure factors, mean-square displacements, radial distribution functions, etc.

[![GitHub](https://img.shields.io/badge/GitHub-AMDAT-181717?logo=github&logoColor=white)](https://github.com/dssimmons-codes/AMDAT)
[![Docs](https://img.shields.io/badge/docs-website-blue)](https://dssimmons-codes.github.io/AMDAT/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17417167.svg)](https://doi.org/10.5281/zenodo.17417167)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Group website](https://img.shields.io/badge/Group%20website-Simmons%20Research%20Group-4285F4?logo=google&logoColor=white)](https://sites.google.com/view/simmonsresearchgroup/software)

## Why AMDAT (at a glance)

* Fast, in-memory engine – load once, analyze many time delays without re-reading files.
* Blocked exponential time spacing – efficient long-timescale dynamics across orders of magnitude.
* Modular data abstractions – trajectory / neighbor / multibody / value lists compose into rich workflows.
* Validated analyses –
  * static and time-resolved structure factors,
  * radial distribution functions,
  * mean-square displacements,
  * neighbor correlations,
  * clustering, and more.
* Plain-text outputs – easy post-processing in Python/Matlab/Excel/VMD/OVITO.

*More details in the [Overview](Manual/overview.md).*

## Quick start (Linux, Conda)

### Build with conda environment (Recommended)
```bash
make conda-setup
conda activate amdat
make
```

### Run
```bash
./AMDAT -i path/to/input.in
```

*More details in [Build/Install](Manual/making_amdat.md).*

## Citation

If you use AMDAT, please cite:
> AMDAT — Amorphous Molecular Dynamics Analysis Toolkit
> DOI: 10.5281/zenodo.17417167

*Also see [CITATION.cff](CITATION.cff).*

## Authors & Maintainers
* David S. Simmons (Lead) · Pierre Kawak · William Drayer
* Past contributors include Michael Marvin, Mark Mackura, Daniel Hunsicker.

## License
GNU GPLv3.0 (see [LICENSE](LICENSE)).
It bundles:
* Voro++ (see [Voro++ LICENSE](third_party/voro++-0.4.6/LICENSE))
* xdrfile (LGPL-3.0, see [xdrfile LICENSE](third_party/xdrfile-1.1b/COPYING))

## Contributing
We welcome issues, discussions, and pull requests.
Please skim [CONTRIBUTING.md](CONTRIBUTING.md) first.
