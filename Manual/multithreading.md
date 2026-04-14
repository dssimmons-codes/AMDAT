<h1>AMDAT Multithreading</h1>

## Overview
HPC workflow in general takes a lot of time, and Molecular Dynamic Simulation is no exception. AMDAT facilitates research by supporting multithreading (currently implemented for a subset of functions) using OpenMP. The list of analysis with multithreading enabled to speedup computation is illustrated below.

---
| Analysis | Description | Maximum Speedup Tested (Compared to Serial) |
|----------|-------------|------------------------|
| `md` | Calculate Mean Displacement | 7x |
| `msd` | Calculate Mean Square Displacement | 10x |
| `msd_2d` | Calculate Mean Square Displacement for two dimensional analysis | 9x |
| `isfs` | Calculate Incoherent Scattering Function | 20x |
| `vhs` | Calculate Van Hove Self | 12x |
| `vhd` | Calculate Van Hove Distinct | 20x |
| `baf` | Calculate Bond Autocorrelatoion Function | 15x |
| `ngp` | Calculate Non Gaussian Parameter | 30x |
| `composition` | Calculate Composition | 7x |
| `rdf` | Calculate Radial Distribution Function | 21x |
| `rnf` | Calculate Radial Count | 23x |


## Inconsistency in performance between analysis

AMDAT is currently implemented as highly object-oriented with its underlying structure is AoS (Array of Structure). Hence, data are organized in adjacent memory blocks, which are prone to Cache Coherency due to hardware architecture of fetching blocks of memory into cache line. In parallel programming, this conflict in the memory is called false-sharing. To fix this, AMDAT uses padding to separate each memory location further from each other, reducing Cache Coherency conflicts.
This approach, however, trades spatial locality from serial for alleviated memory conflict, preventing achieving perfect speedup. Therefore, speedups are inconsistent between analysis types due to the nature of how we store their data.

## System Requirements

To use multithreading effectively in AMDAT, make sure the following requirements are met:

- AMDAT is compiled with **OpenMP** support
- Your system has a multi-core CPU
- The OpenMP runtime is available in your environment

## Analysis Requirements

Parallelism is introduced to reduce the runtime of time-based loops. 
To avoid oversubscription and ensure efficient resource utilization, 
the optimal number of threads is allocated can be given by:

$$N_{\text{threads}} = \min \left( N_{\text{iterations}},\; N_{\text{CPU}} \right).$$

where:
- $$\ N_{\text{iterations}} \$$ is the number of independent time-based iterations.
- $$\ N_{\text{CPU}} \$$ is the maximum number of hardware threads (logical CPUs) available on the node.

## Usage

Running AMDAT is very simple with conda build, because OpenMP is already preinstalled. The flag `-n` <Number_of_Threads> will allocate the desired number of parallel regions for your analysis.

### Example with msd

```bash
./AMDAT -n 8 -i ./testfiles/msd.in

