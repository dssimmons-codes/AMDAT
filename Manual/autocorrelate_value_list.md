<h1>autocorrelate_value_list</h1>

<h2>Function</h2>

Computes the normalized self time autocorrelation function of a value list. If a quantity $A$ is stored in the value_list, then its autocorrelation function $C_A(\Delta t)$ is given by

$$  C_A(\Delta t) = \frac{1}{\langle A \rangle^2} \sum_{k=1}^{S} \frac{1}{N(s)} \sum_{i=1}^{N(s)} (\{A_i}(s_k+\Delta t)\{A_i}(s_k))   $$

where S is the number of start times employed, N is the number of particles, $\{A}_i(t)$ is the value stored in the list for particle $i$ at time $t$, and $\Delta t$ is a timegap.

<h2>Syntax</h2>

Syntax is as follows.  

_autocorrelate\_value\_list \<name of output file\> \<name of value\_list\>_
