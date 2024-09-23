<h1>raf</h1>

<h2>Function</h2>

Computes the reorientation autocorrelation function for a list of multibodies. All multibodies in the list must consist of exactly 2 bodies so that they define a single vector (for example corresponding to a bond or other some intramolecular vector)

$$  C(\Delta t) =  P_n \left[\sum_{k=1}^{S} \frac{1}{N(s)} \sum_{i=1}^{M(s)} (\hat{r_i}(s_k+\Delta t) \cdot \hat{r_i}(s_k)) \right]   $$

where S is the number of start times employed, M is the number of bonds, $\hat{r}_i(t)$ is the unit vector pointing from the first to second particle within the multibody at time $t$, $\Delta t$ is a timegap, and $P_n$ is the Legendre polynomial of the $n^{th}$ order.

<h2>Syntax</h2>

_raf \<output filename\> \<name of multibody\_list to analyze\> \<Legendre polynomial order to employ: either "1" or "2"\> \<optional: "xyz" or "xy" or "xz" or "yz" or "x" or "y" or "z"\>_

The optional argument specifies which spatial dimensions are including in the calculation (i.e. in the dot product). This defaults to "xyz" if not specified, in which case the full three-dimensional dot product is employed. If, for example, "xy" is specified, only the x and y dimensions of the unit vector are included in the dot product.
