<h1>compare_gaussian</h1>

Find points where the self Van Hove crosses the Gaussian approximation to the self Van Hove.

_compare\_gaussian \<output file\> \< displacement time index (optional)\>_

_\<target\>_

Requires prior calculation of self van hove in the same run. If only two arguments are given, this analysis first calculates the non-Gaussian parameter and used the time spacing corresponding to the peak in the nGP. If a third argument is given, self Van Hove will be compared to the Gaussian approximation at the time index given.
