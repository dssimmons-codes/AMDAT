<h1>mean_closest_distance</h1>

Requires specification of two target lists. Calculates the average distance from particles in list 1 to the closest member of list 2.

_mean\_closest\_distance \<output file name\> \<timescheme>_

\<timescheme\> selects which times are averaged over. If timescheme = -1, loop over all times. If timescheme is zero or positive, only use one time per block, with the value setting the time index offset from the beginning of the block.
