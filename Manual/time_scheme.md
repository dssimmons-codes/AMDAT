<h1>Time Schemes</h1>

AMDAT is currently capable of handling two schema for the time-spacing with which snapshots are output from an MD simulation to the system trajectory: _linear_ and _blocked exponential_. The choice between these two schema has major implications for how dynamical quantities are analyzed. In particular, AMDAT assumes _stationarity_ for most dynamical calculations, such that dynamics are averaged over many start times. The handling of which times are employed as start times is quite different in the cases of the two schema. Formally, if we have some correlation function $C(\Delta t)$, this function is typically averaged over $S$ start times $s_j$, as

$$C(\Delta t)=\sum_{j=1}^{S}(\langle A(s_j)B(s_j+\Delta t) \rangle$$

In addition to inferring different values for the times read in, linear and blocked exponential time schemes lead to different choices of the set of start times s_j to be employed in stationary dynamical calculations.

<h2>Linear time spacing</h2>

In linear time spacing, every pair of sequential snapshots is separated by the same $\Delta \tau$. Formally, snapshots are written out at times following the equation
$$t=k\Delta \tau$$
where $k=0,1,2,...T-1$, and were $T$ is the total number of frames.

For dynamical quantities, every time is employed as a starting time for dynamical calculations, and every possible pair of snapshots is employed.

<h2>Blocked exponential time spacing</h2>

The _\<time scheme\>_ line is formatted as follows

_\<type\> \<arguments\>_

Formatting for each type is:

typeargs

_linear \<number of timesteps\> \<time unit\>_

Here _\<number of timesteps\>_ is the number of frames in the trajectory and \<time unit\> is the time between each frame.

_Exponential \<# of exponential blocks\> \<timesteps per block\> \<exp base\> \<frt\> \<first exp\> \<time unit\>_

In an exponential timescheme, frame times are given by the following formula:

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
