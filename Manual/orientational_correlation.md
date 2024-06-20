<h1>Orientational Correlation</h1>

Calculates the orientational correlation of multibodies with an externally defined vector. All multibodies in the multibody\_list specified must consist of exactly 2 bodies so that they each define a single vector.

_orientational\_correlation_ \<output filnename\> \<name of multibody\_list to analyze\> \<x component of external vector\> \<y component of external vector\> \<z component of external vector\>_

Computes the mean value of the second Legendre Polynomial of the dot product between the multibody vector and the external vector, averaged over all multibodies in the list and over all times.
