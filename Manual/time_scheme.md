<h1>Time Schemes</h1>

AMDAT is currently capable of handling two schema for the time-spacing with which snapshots are output from an MD simulation to the system trajectory: _linear_ and _blocked exponential_. The choice between these two schema has major implications for how dynamical quantities are analyzed. In particular, AMDAT assumes _stationarity_ for most dynamical calculations, such that dynamics are averaged over many start times. The handling of which times are employed as start times is quite different in the cases of the two schema. Formally, if we have some correlation function $C(\Delta t)$, this function is typically averaged over $S$ start times $s_j$, as

$$C(\Delta t)=\sum_{j=1}^{S}(\langle A(s_j)B(s_j+\Delta t) \rangle$$

In addition to inferring different values for the times read in, linear and blocked exponential time schemes lead to different choices of the set of start times s_j to be employed in stationary dynamical calculations.

In AMDAT, the time scheme is specified in the input file header section via the _\<time scheme\>_ line, which is formatted as follows

_\<type\> \<arguments\>_

Here _\<type\>_ is _linear_ or _exponential_. The sections for each of these time schemes below describe the required arguments. 

<h2>Linear time spacing</h2>

In linear time spacing, every pair of sequential snapshots is separated by the same $\Delta \tau$. Formally, snapshots are written out at times following the equation
$$t=k\Delta \tau$$
where $k=0,1,2,...T-1$, and were $T$ is the total number of frames.

For dynamical quantities, every time is employed as a starting time for dynamical calculations, and every possible pair of snapshots is employed. This means that the number of start times $S$ is a function of the timegap $\Delta t$ considered, as
$$S(\Delta t)=T-\frac{\Delta t}{\Delta \tau}$$
This has two implications. First, in general statistical strength at a minimum drops linearly as per the equation above at larger timegaps simply because of the reduced number of starting times available to average over. Second, the drop is actually somewhat worse than this because time windows for larger time gaps become increasingly overlapping and hence correlated. For most purposes, linear time schemes are therefore not a very efficient way of writing out or analyzing trajectories for the purposes of  dynamical correlation functions in stationary systems.

The time scheme line for linear time schemes is formatted as follows.

_linear \<number of timesteps\>_ $<\Delta t>$

Here _\<number of timesteps\>_ is the number of frames in the trajectory and $<\Delta t>$ is the time separating sequential frames.

<h2>Blocked exponential time spacing</h2>

In blocked exponential time spacing, frames are written out in blocks, with times in each block spaced exponentially. If the underlying timestep of the simulation is $\Delta \tau$, then the times at which frames are written in the first block is given by

frame indices written to the trajectory are broadly speaking governed by an exponential $b^k$, where $b$ is an exponential base and $k$ is an integer iterator. However, this equation does not generally yield intergers, and in some cases leads to spacings between time steps that are less than one timestep apart. This is resolved by, first, taking the floor of $b^k$ and, second, defaulting to linear spacing for k such that $b^k < k$.


The time scheme line for blocked exponential time schemes is formatted as follows.

_Exponential \<# of exponential blocks\> \<timesteps per block\> \<exp base\> \<frt\> \<first exp\> \<time unit\>_

In an exponential timescheme, frame times are given by the following formula:

$$t=

where m is the block number (m = 0,1,2,3,…,M) and k is the frame (k=0,1,2,3… if _frt_ = 0 and k = 1,2,3,… if _frt_ = 1), N is the number of timesteps per block. If _frt_ is zero, then an extra initial time frame at time zero is used at the beginning of the trajectory. If _frt_ is one, then the first time frame is at time (dt)b^_\<first exponent\>_. This has significant implications for how time blocks are handled. First exponent is usually zero or one and is the first value of the exponent for the power law progression in time. Deltat is the time unit. For clarity, below is a segment of pseudocode describing the list of simulation times corresponding to this scheme. Also packaged with AMDAT is an example LAMMPS input file yielding a trajectory file corresponding to this time scheme.

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
  
_Snapshot_ (no args)

| type        | args                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|-------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| linear      | <number of timesteps> <time unit>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Here <number of timesteps> is the number of frames in the trajectory and <time unit> is the time between each frame. ||
| Exponential | <# of exponential blocks> <timesteps per block> <exp base> <frt> <first exp> <time unit>                |
|             | In an exponential timescheme, frame times are given by the following formula:  where m is the block number (m = 0,1,2,3,…,M) and k is the frame (k=0,1,2,3… if frt = 0 and k = 1,2,3,… if frt = 1), N is the number of timesteps per block. If frt is zero, then an extra initial time frame at time zero is used at the beginning of the trajectory. If frt is one, then the first time frame is at time (dt)b^<first exponent>. This has significant implications for how time blocks are handled. First exponent is usually zero or one and is the first value of the exponent for the power law progression in time. Deltat is the time unit. For clarity, below is a segment of pseudocode describing the list of simulation times corresponding to this scheme. Also packaged with AMDAT is an example LAMMPS input file yielding a trajectory file corresponding to this time scheme.  block_starttime=0; for(blockii=0;blockii<#_of_exponentials;blockii++) {   for(expii=1;expii<=timesteps_per_block;expii++)   {     timeii++;     if(pow(exp_base,expii-1+first_exponent) <= expii)     {       time[timeii] = block_starttime+expii*time_unit;     }     else     {       time[timeii] = block_starttime+floor(pow(exp_base,expii-1+first_exponent))*time_unit;     }   }   block_starttime = timelist[timeii]; } |
| Snapshot | no args |
|          | Snapshot is used for a single frame trajectory.|

<h2>Analysis block</h2>
