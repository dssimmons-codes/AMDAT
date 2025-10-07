<h1>crosscorrelate_value_lists</h1>

<h2>Function</h2>

Computes the normalized time-dependent crosscorrelation function between two value lists. If a quantity $A$ is stored in the first value_list and the a quantity $B$ in the second, then their crosscorrelation function $C_{AB}(\Delta t)$ is given by

$$  C_{AB}(\Delta t) = \frac{1}{\langle A \rangle \langle B \rangle} \sum_{k=1}^{S} \frac{1}{N(s)} \sum_{i=1}^{N(s)} (\{B_i}(s_k+\Delta t)\{A_i}(s_k))   $$

where S is the number of start times employed, N is the number of particles in list A, $\{A}_i(t)$ is the value stored in the list for particle $i$ at time $t$, and $\Delta t$ is a timegap. This tool fundamentally loops over atoms defined in value_list A. Thus, for this tool to work properly, a value in list B must be defined for each particle in list A. 

<h2>Syntax</h2>

Syntax is as follows.  

_crosscorrelate\_value\_lists \<correlation type keyword\> \<output file\> \<name of valuelist A\> \<name of valuelist B>_

_\<correlation type keyword\>_ may be either _dynamic_ or _static_. When dynamic, the full calculation above is performed for all possible $\Delta t$. When static, only $\Delta t = 0$ is computed.
