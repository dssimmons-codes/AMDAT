<h1>Time Schemes</h1>

AMDAT is currently capable of handling three schema for the time-spacing with which snapshots are output from an MD simulation to the system trajectory: _snapshot_, _linear_ and _blocked exponential_. The choice between these two schema has major implications for how dynamical quantities are analyzed. In particular, AMDAT assumes _stationarity_ for most dynamical calculations, such that dynamics are averaged over many start times. The handling of which times are employed as start times is quite different in the cases of the two schema. Formally, if we have some correlation function $C(\Delta t)$, this function is typically averaged over $S$ start times $s_j$, as

$$C(\Delta t)=\sum_{j=1}^{S}(\langle A(s_j)B(s_j+\Delta t) \rangle$$

In addition to inferring different values for the times read in, linear and blocked exponential time schemes lead to different choices of the set of start times s_j to be employed in stationary dynamical calculations.

In AMDAT, the time scheme is specified in the input file header section via the _\<time scheme\>_ line, which is formatted as follows

```
<type> <arguments>
```

Here ``<type>`` is _linear_ or _exponential_. The sections for each of these time schemes below describe the required arguments. 

All trajectory read-in methods begin reading in at the start of the file. There is currently no native way of 'skipping' some number of timesteps at the beginning of the file.

<h2>Snapshot</h2>

_snapshot_ reads in a single initial frame of a trajectory. The syntax is as follows.

```
snapshot
```

No arguments are required or allowed.

<h2>Linear time spacing</h2>

<h3>Syntax</h3>

The time scheme line syntax for linear time schemes is formatted as follows.

```
linear <number of timesteps> <Δ𝜏>
```

Here ``<number of timesteps>`` is the number of frames in the trajectory and ``<Δ𝜏>`` is the time separating sequential frames.

<h3>Frame times</h3>

In linear time spacing, every pair of sequential snapshots is separated by the same $\Delta \tau$. Formally, snapshots are written out at times following the equation
$$t=k\Delta \tau$$
where $k=0,1,2,...T-1$, and were $T$ is the total number of frames.

<h3>Time spacings used in dynamical analysis</h3>

For dynamical quantities, every time is employed as a starting time for dynamical calculations, and every possible pair of snapshots is employed. This means that the number of start times $S$ is a function of the timegap $\Delta t$ considered, as
$$S(\Delta t)=T-\frac{\Delta t}{\Delta \tau}$$
This has two implications. First, in general statistical strength at a minimum drops linearly as per the equation above at larger timegaps simply because of the reduced number of starting times available to average over. Second, the drop is actually somewhat worse than this because time windows for larger time gaps become increasingly overlapping and hence correlated. For most purposes, linear time schemes are therefore not a very efficient way of writing out or analyzing trajectories for the purposes of  dynamical correlation functions in stationary systems.

<h2>Blocked exponential time spacing</h2>

<h3>Syntax</h3>

The time scheme line for blocked exponential time schemes is formatted as follows in terms of the variables in the discussion below.

_exponential_ $< I >$ $< K >$ $< b >$ \<frt\> $< a_0 >$ $< \Delta t >$

In normal language, this can be read as follows.

```
exponential <# of exponential blocks> <size of block> <exp base> <frt> <first exponential index> <time step>
```

The manner in which these parameter choices generate frame times and time gaps is described below.

<h3>Frame times</h3>

In blocked exponential time spacing, frames are written out in blocks, with times in each block spaced exponentially. Generally 

If the underlying timestep of the simulation is $\Delta \tau$, then the times at which frames are written in the first block are broadly speaking governed by an exponential $b^k$, where $b$ is an exponential base and $k$ is an integer iterator. However, this equation does not generally yield intergers, and in some cases leads to spacings between time steps that are less than one timestep apart. This is resolved by, first, taking the floor of $b^k$ and, second, defaulting to linear spacing for k such that $b^k < k$.

Within the first block of a blocked exponential time scheme, the $k^{th}$ frame is written out at a time $t(k)$ given by

$t(k=0)=0$

$t(0 < k \leq K)=\Delta \tau k$ when $\left\lfloor{(b^{k-1+a_0})}\right\rfloor < k$

$t(0 < k \leq K)=\Delta \tau \left\lfloor{(b^{k-1+a_0})}\right\rfloor$ when $\left\lfloor{(b^{k-1+a_0})}\right\rfloor > k$

where $k=0,1,2,...K$ and where $K$ is called the "block size". $a_0$ is the first index of the exponent. It is most commonly set to 0 but may be set to 1.

At the end of this block, the time spacing repeats the exponential algorithm above, such that frames in the next block are written out at times as follows

$t(K < k \leq 2K))=\Delta \tau \left\lfloor{(b^{K-1+a_0})}\right\rfloor + \Delta \tau k$ when $\left\lfloor{(b^{k-K-1+a_0})}\right\rfloor < k - K$

$t(K < k \leq 2K)=\Delta \tau \left\lfloor{(b^{K-1+a_0})}\right\rfloor + \Delta \tau \left\lfloor{(b^{k-K-1+a_0})}\right\rfloor$ when $\left\lfloor{(b^{k-K-1+a_0})}\right\rfloor > k - K$

In a run with $I$ sequential blocks, this is then repeated with the general pattern

$t(k > 0)=\Delta \tau \left\lfloor{(b^{i K-1+a_0})}\right\rfloor + \Delta \tau k$ when $\left\lfloor{(b^{k-i K-1+a_0})}\right\rfloor < k-iK$

$t(k > 0)=\Delta \tau \left\lfloor{(b^{i K-1+a_0})}\right\rfloor + \Delta \tau \left\lfloor{(b^{k-i K-1+a_0})}\right\rfloor$ when $\left\lfloor{(b^{k-i K-1+a_0})}\right\rfloor > k-iK$

where $i = 0, 1, 2, ...I-1$ and where $i$ is the index of the block.



_frt_ should in most cases be set to 0. The value of this quantity does not impact file read-in but does impact which pairs of frames are used for a given time-spacing, as discussed below. When _frt_ = 0, AMDAT essentially treats the last frame of each block as being the same as the first frame of the next block, which is almost always preferrable in maximizing data use. Using _frt_ = 1 instead treats the first frame of a given block as being one frame after the last frame of the prior block, which both reduces the longest time gap accessible by one increment and slightly decreases data efficiency. 

For clarity, below is a segment of pseudocode describing the list of simulation times corresponding to this scheme. Also packaged with AMDAT is an example LAMMPS input file yielding a trajectory file corresponding to this time scheme.

    block_starttime=0;
    for(blockii=0;blockii<#_of_exponentials;blockii++)
    {
      for(expii=1;expii<=timesteps_per_block;expii++)
      {
        timeii++;
        if(pow(exp_base,expii-1+first_exponent) <= expii)
        {
          time[timeii] = block_starttime+expii*time_unit;
        }
        else
        {
          time[timeii] = block_starttime+floor(pow(exp_base,expii-1+first_exponent))*time_unit;
        }
      }
      block_starttime = timelist[timeii];
    }

<h3>Time spacings used in dynamical analysis</h3>

When computing dynamical quantities with blocked time schemes, AMDAT by default employs only the beginning frame of each block as a start time s. This means that the number of start times over which dynamical quantities is averages is held constant (and equal to the number of blocks K) for all time gaps $\Delta t$ that fit within a single block. For these time gaps, AMDAT will compute dynamical quantities at time gaps that can be obtained by the difference between the time of any later frame in the block and the time of the initial frame of the block. This allows for averaging over all $K$ blocks. For time gaps that are outside of a block, AMDAT will compute dynamical quantities only at time gaps equal to an even multiple of the block length, such that the dynamics are computed between beginning (or ending) frames of distinct blocks. This means that, for time gaps outside of a single block, the number of start times s over which results are averaged will linearly drop from K to 1 as the time gap increases to the full duration of the trajectory. Practically, this provides efficient access to reasonable statistical quality over many orders of magnitude in time, while keeping storage and computational requirements relatively low.
