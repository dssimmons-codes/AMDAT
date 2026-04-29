<h1>AMDAT Multithreading</h1>

## Overview
High-performance computing workflows can be computationally expensive, and molecular dynamics analysis is no exception. To accelerate selected analyses, AMDAT supports multithreading with OpenMP. The analyses that currently support multithreading are listed below.

---
| Analysis | Description | Maximum Speedup Tested (Compared to Serial) |
|----------|-------------|------------------------|
| [`md`](mean_displacement.md) | Calculate Mean Displacement | 7x |
| [`msd`](msd.md) | Calculate Mean Square Displacement | 10x |
| [`msd_2d`](msd_2d.md) | Calculate Mean Square Displacement for two-dimensional analysis | 9x |
| [`isfs`](isfs.md) | Calculate Incoherent Scattering Function | 20x |
| [`vhs`](vhs.md) | Calculate self part of Van Hove Function | 12x |
| [`vhd`](vhd.md) | Calculate distinct part of Van Hove Function | 20x |
| `baf` | Calculate Bond Autocorrelation Function | 15x |
| [`ngp`](ngp.md) | Calculate Non Gaussian Parameter | 30x |
| [`composition`](composition.md) | Calculate Composition | 7x |
| [`rdf`](rdf.md) | Calculate Radial Distribution Function | 21x |
| [`rnf`](rnf.md) | Calculate Radial Count | 23x |


## Inconsistency in performance between analysis

AMDAT is currently implemented with a highly object-oriented design and an array-of-structures (AoS) memory layout. Hence, data are organized in adjacent memory blocks, making them more prone to Cache Coherency issue. In parallel system, this conflict results in false-sharing, meaning different threads try to have writing access to one common cache line, despite aiming at independent memory location. To fix performance bottleneck, AMDAT uses padding to separate each memory location further from each other, reducing Cache Coherency conflicts. This approach, however, trades spatial locality from serial implementation for mitigated memory conflict in parallel tasks, making perfect speed-up realistically unachievable. Therefore, speed-ups are inconsistent between analyses due to the differences in how we store the data.

## System Requirements

To use multithreading effectively in AMDAT, make sure the following requirements are met:

- Your system has a multi-core CPU
- The OpenMP runtime is available in your environment (This is already handled by Conda)

## Analysis Requirements

Parallelism is introduced to reduce the runtime of time-based loops. 
However, to avoid oversubscription and ensure efficient resource utilization, 
the optimal number of threads is allocated can be given by:

$$N_{\text{threads}} = \min\left(N_{\text{iterations}},\;N_{\text{CPUs}} \right).$$

where:
- $$\ N_{\text{CPUs}} \$$ is the maximum number of hardware threads (logical CPUs) available on the node.
- $$\ N_{\text{iterations}} \$$ is the number of independent time-based iterations. In AMDAT's input files, it is highlighted as shown below:

> system  
> xyz_log  
> ./testfiles/traj.xyz ./testfiles/log.lammps  
> exponential **`70`** 77 1.05 0 0 .01  
> polymer 50  
>
> ...

## Usage

Running AMDAT is very simple with conda build, because OpenMP is already preinstalled. The flag `-n` <Number_of_Threads> will allocate the desired number of parallel regions for your analysis.

### Example with msd

```bash
./AMDAT -n 8 -i ./testfiles/msd.in

