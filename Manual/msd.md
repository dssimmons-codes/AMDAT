<h1>msd</h1>
<h2>Function</h2>

Calculates mean square displacement, as follows.

$$  \langle r(\Delta t)^2\rangle =  \sum_{k=1}^{S} \frac{1}{N(s)} \sum_{i=1}^{N(s)} (\vec{r_i}(s_k+\Delta t)-\vec{r_i}(s_k))^2   $$

where S is the number of start times employed, N is the number of particles, $\vec{r}_i(t)$ is the position of particle $i$ at time $t$, and $\Delta t$ is a timegap.

<h2>Syntax</h2>

_msd \<output file\>_

_\<target\>_


