<h1>ngp</h1>

<h2>Function</h2>

Calculates non-Gaussian parameter of the mean squared displacement.

$$  ngp(\Delta t) =  \frac{3 \sum_{k=1}^{S} \frac{1}{N(s)} \sum_{i=1}^{N(s)} (\boldsymbol{r_i}(s_k+\Delta t)-\boldsymbol{r_i}(s_k))^4}{5 \langle r(\Delta t)^2\rangle^2} - 1  $$

where S is the number of start times employed, N is the number of particles, $\boldsymbol{r}_i(t)$ is the position of particle $i$ at time $t$, $\Delta t$ is a timegap, and the mean square displacement $\langle r(\Delta t)^2\rangle$ is computed as described in the [msd](msd.md) page.

<h2>Syntax</h2>

_ngp \<output file\>_  
_\<target\>_


