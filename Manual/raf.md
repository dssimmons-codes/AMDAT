<h1>raf</h1>

Computes the reorientation autocorrelation function for a list of multibodies. All multibodies in the list must consist of exactly 2 bodies so that they define a single vector (for example corresponding to a bond or other some intramolecular vector)

_raf \<output filename\> \<name of multibody\_list to analyze\> \<Legendre polynomial order to employ: either "1" or "2"\> \<optional: "xyz" or "xy" or "xz" or "yz" or "x" or "y" or "z"\>_

The optional argument defaults to "xyz" if not selected. Selects the plane or dimension in which reorientation is calculated. "xyz" gives the usual isotropic reorientation in all dimensions.
