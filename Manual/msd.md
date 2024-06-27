<h1>msd</h1>
<h2>Function</h2>

Calculates mean square displacement, as follows.

$$  \langle r(\Delta t)^2\rangle = \sum_{k=1}^{S} \sum_{i=1}^{N(s)} (\boldsymbol{r_i}(s_k+\Delta t)-\boldsymbol{r_i}(s_k))^2   $$

where S is the number of start times employed, N is the number of particles, $\boldsymbol{r}_i(t)$ is the position of particle $i$ at time $t$, and $\Delta t$ is a timegap.

<h2>Syntax</h2>

_msd \<output file\>_

\<target\>


