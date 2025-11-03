# AMDAT — Amorphous Molecular Dynamics Analysis Toolkit

> High-performance analysis of molecular dynamics (MD) trajectories for amorphous, glass-forming, and polymer systems.

AMDAT is a C++ toolkit that loads trajectories into memory for **rapid** analysis, with robust **atom selection**, **time-resolved** statistics, and modular data objects.
It reads common LAMMPS trajectory formats and GROMACS `.xtc`.
It provides a wide variety of high-performance analyses integral to molecular modeling studies, such as clustering, spatial decomposition, and calculating time-resolved structure factors, mean-square displacements, radial distribution functions, etc.

[![GitHub](https://img.shields.io/badge/GitHub-AMDAT-181717?logo=github&logoColor=white)](https://github.com/dssimmons-codes/AMDAT)
[![Docs](https://img.shields.io/badge/docs-website-blue)](https://dssimmons-codes.github.io/AMDAT/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17417166.svg)](https://doi.org/10.5281/zenodo.17417166)
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

### Clone the repo from GitHub

Prereqs:
* Conda
* Set up GitHub ssh keys on your machine

```bash
git clone git@github.com:dssimmons-codes/AMDAT.git
cd AMDAT/
```

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

## Highlights

<table>
  <tr>
    <td align="center" width="25%">
      <a href="Manual/displacement_list.md">
        <img src="Manual/assets/highlights/mobile_particles.png" alt="Identification of highly mobile particles" width="100%">
      </a><br/>
      <sub>Identification of highly mobile particles<br/>(image: Pierre Kawak)</sub>
    </td>
    <td align="center" width="25%">
      <a href="Manual/displacement_list.md">
        <img src="Manual/assets/highlights/displacements_2d.png" alt="Particle displacements in binary 2D LJ" width="100%">
      </a><br/>
      <sub>Particle displacements (binary 2D LJ)<br/>(image: Pierre Kawak)</sub>
    </td>
    <td align="center" width="25%">
      <a href="Manual/create_bin_list.md">
        <img src="Manual/assets/highlights/color_by_distance.png" alt="Sorting by distance from nanoparticle/interface" width="100%">
      </a><br/>
      <sub>Sorting by distance from a nanoparticle/interface<br/>(image: Pierre Kawak)</sub>
    </td>
    <td align="center" width="25%">
      <a href="Manual/n_fold.md">
        <img src="Manual/assets/highlights/2d_hexatic.png" alt="6-fold orientational order parameter in binary 2D LJ" width="100%">
      </a><br/>
      <sub>6-fold orientational order parameter (binary 2D LJ)<br/>(image: Daniel Hunsicker)</sub>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td align="center" width="50%">
      <a href="Manual/msd.md">
        <img src="Manual/assets/highlights/msd.jpg" alt="Mean-squre displacments for a bead-spring polymer (image by Sean Hung, adapted from  Hung, Patra, Meenakshisundaram, Mangalara, Simmons, Soft Matter, 15 (2019) 1223-1242. doi: 10.1039/C8SM02051E." width="100%">
      </a><br/>
      <sub>
        Mean-square displacement for bead-spring polymer<br/>
        doi: <a href="https://doi.org/10.1039/C8SM02051E">10.1039/C8SM02051E</a><br/>
        (image: Sean Hung)
      </sub>
    </td>
    <td align="center" width="50%">
      <a href="Manual/isfs.md">
        <img src="Manual/assets/highlights/isfs.gif" alt="Intermediate scattering functions for a bead-spring polymer (image by Sean Hung, adapted from  Hung, Patra, Meenakshisundaram, Mangalara, Simmons, Soft Matter, 15 (2019) 1223-1242. doi: 10.1039/C8SM02051E." width="100%">
      </a><br/>
      <sub>
        Time-resolved structure factor peak for bead-spring polymer<br/>
        doi: <a href="https://doi.org/10.1039/C8SM02051E">10.1039/C8SM02051E</a><br/>
        (image: Sean Hung)
      </sub>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td align="center" width="100%">
      <a href="Manual/overview.md">
        <img src="Manual/assets/highlights/cg_map.png" alt="AMDAT-based post-simulation mapping of atomistic polystyrene repeat units to segmental center of mass calculations (left), and identificiation of string-like cooperative rearrangements (a la doi.org/10.1103/PhysRevLett.80.2338) (middle two), visualized across two timesteps (red and blue in right image).  Images by Sean Hung. Computed for simulations reported in Jui Hsiang Hung, David S Simmons, Do String-like Cooperative Motions Predict Relaxation Times in Glass-Forming Liquids?, Journal of Physical Chemistry B, 124, 1 (2020) 266-276. doi: 10.1021/acs.jpcb.9b09468." width="100%">
      </a><br/>
      <sub>
        Coarse-graining of atomistic polystyrene and identification of string-like cooperative rearrangements<br/>
        doi: <a href="https://doi.org/10.1021/acs.jpcb.9b09468">10.1021/acs.jpcb.9b09468</a>
      </sub>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td align="center" width="100%">
      <a href="Manual/overview.md">
        <img src="Manual/assets/highlights/string_rdf.jpeg" alt="Radial distribution functions for (a) bead−spring polymer; (b) OTP; (c,f,i) binary LJ glass former; (d,g,j) Cu4Ag6; (e,h,k) SiO2. Reproduced from Jui Hsiang Hung, David S Simmons, Do String-like Cooperative Motions Predict Relaxation Times in Glass-Forming Liquids?, Journal of Physical Chemistry B, 124, 1 (2020) 266-276. doi:10.1021/acs.jpcb.9b09468." width="100%">
      </a><br/>
      <sub>
        Radial distribution functions for bead-spring polymer (top left), OTP (top right), binary LJ glass former (left column), Cu4Ag6 (middle column), and SiO2 (right column)<br/>
        doi: <a href="https://doi.org/10.1021/acs.jpcb.9b09468">10.1021/acs.jpcb.9b09468</a>
      </sub>
    </td>
  </tr>
</table>

<table>
  <tr>
    <td align="center" width="100%">
      <a href="Manual/overview.md">
        <img src="Manual/assets/highlights/struct.png" alt="Clockwise from top left, structure factors for Bead-spring polymer, binary Lennard Jones glass-former, binary copper-silver alloy, OTP (atomistic structure factor in red and ring-center-of mass structure factor in blue, see inset), and SiO2. Computed for simulations reported in Jui Hsiang Hung, David S Simmons, Do String-like Cooperative Motions Predict Relaxation Times in Glass-Forming Liquids?, Journal of Physical Chemistry B, 124, 1 (2020) 266-276. doi:10.1021/acs.jpcb.9b09468." width="100%">
      </a><br/>
      <sub>
        Structure factors for (clockwise from top left) bead-spring polymer, binary LJ glass former, binary copper-silver allow, OTP, and SiO2<br/>
        doi: <a href="https://doi.org/10.1021/acs.jpcb.9b09468">10.1021/acs.jpcb.9b09468</a>
      </sub>
    </td>
  </tr>
</table>

## Citation

If you use AMDAT, please cite:
> AMDAT — Amorphous Molecular Dynamics Analysis Toolkit
> DOI: 10.5281/zenodo.17417167

*Also see [CITATION.cff](CITATION.cff).*

## Authors & Maintainers
* Lead: [David S. Simmons](https://github.com/dssimmons)
* Maintainers: [Pierre Kawak](https://github.com/pkawak) · [William Drayer](https://github.com/drayer587)
* Contributors: [Michael Marvin](https://github.com/magico13) · Mark Mackura · Daniel Hunsicker

## License
GNU GPLv3.0 (see [LICENSE](LICENSE)).
It bundles:
* Voro++ (see [Voro++ LICENSE](third_party/voro++-0.4.6/LICENSE))
* xdrfile (LGPL-3.0, see [xdrfile LICENSE](third_party/xdrfile-1.1b/COPYING))

## Contributing
We welcome issues, discussions, and pull requests.
Please skim [CONTRIBUTING.md](CONTRIBUTING.md) first.
