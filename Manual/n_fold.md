<h1>n_fold</h1>

Calculates mean and time dependence of n-fold order parameter and writes maps to files.

_nfold \<output file name\> \<n-fold order\> \<plane:"xy","yz","xz"\> \<neighbor cutoff dist.\> \<sigma matrix file\> \<analysis\_value\_list name\> \<optional 2: map file name stem\> \<optional 2: first map time\> \<optional 2: last map time\>_

_\<target\>_

_\<sigma matrix file\>_ refers to a file that describes the cross interaction σ values. It is formatted as described below:

_\<arbitraty name for species 1\> \<σ <sub>11</sub> \> \<σ <sub>12</sub> \> \<σ <sub>13</sub> \> … \<σ <sub>1n</sub> \>_

_\<arbitraty name for species 2\> \<σ <sub>21</sub> \> \<σ <sub>22</sub> \> \<σ <sub>23</sub> \> … \<σ <sub>2n</sub> \>_

_._

_._

_._

_\<arbitrary name for species n\> \<σ <sub>n1</sub> \> \<σ <sub>n2</sub> \> \<σ <sub>n3</sub> \> … \<σ <sub>nn</sub> \>_
