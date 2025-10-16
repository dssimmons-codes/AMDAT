<h1>msd_2d</h1>

<h2>Function</h2>

Calculates mean square displacement in two dimensions

$$  \langle r^{AB}(\Delta t)^2\rangle = \sum_{k=1}^{S} \sum_{i=1}^{N(s)} (\boldsymbol{r}^{AB}_i(s_k+\Delta t)-\boldsymbol{r}^{AB}_i(s_k))^2   $$

where S is the number of start times employed, N is the number of particles, $\boldsymbol{r}^{AB}_i(t)$ is the position of particle $i$ at time $t$ in the $AB$ plan, and $\Delta t$ is a timegap.

<h2>Syntax</h2>

```
msd <output file> <plane:"xy","yz","xz">
<target>
```

Here the choice of plane specifies which plane is employed for $AB$ in the equation above.
