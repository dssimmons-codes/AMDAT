<h1>strings</h1>

_strings \<output file\> \< timestep index\> \<threshold\> \<sigma matrix file\> \<trajectory\_list\_name\>_

_\<target\>_

Calculates mean string length, mean number of strings, string order parameter, and distribution of string lengths. The target should be a trajectory list of mobile particles. This can be produced via the "find\_fast" command. Threshold sets the threshold (as a multiple of sigmaij for each particle pair) within which one particle is considered to have replaced another. A standard value for this is 0.6. Also creates a trajectory list called \<_trajectory\_list\_name\>_ containing all particles in strings for further analysis.

Sigma matrix file specifies the sigma­­ij for all interaction types as described in the **n\_fold** command.
