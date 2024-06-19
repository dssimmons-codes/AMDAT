<h1>threshold_multibody_list</h1>

Creates a new multibody list by downselecting a multibody list based on the number of bodies in each multibody.

_threshold\_multibody\_list \<name of new multibody\_list to create\> \<name of existing multibody\_list from which to down-select\> \<comparison keyword\> \<comparison value\>_

Comparison keyword is either "greater" or "less", and "value" sets what number of bodies is used as a comparator for this keyword. For example, "greater 6" would select only those multibodies comprised of 7 or more particles.
