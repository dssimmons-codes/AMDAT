<h1>Gyration Radius</h1>

Computes the mean multibody gyration radius of for the multibodies in a specified ``multibody_list``. If a list of multibodies is defined to correspond to the molecules of a given species for example, this will report the mean gyration radius of molecules of that species, averaged over all molecules and all times.

``gyration_radius <output filename> <multibody_list to analyze>``

Note that this analysis method can also be employed to compute the end to end distance of a molecule. To do so, create a multibody consisting of only (two total) end atoms of the molecule, one from each end. For a two particle system of this kind, Rg is simply 1/2 the distance between the two particles.
