<h1>n_fold</h1>

Calculates mean and time dependence of n-fold order parameter and writes maps to files.

```
nfold <output file name> <n-fold order> <sigma matrix file> <plane:"xy","yz","xz"> <neighbor cutoff distance> <analysis_value_list name> <optional 2: map file name stem> <optional 2: first map time> <optional 2: last map time>
<target>
```

`<sigma matrix file>` refers to a file that describes the cross interaction $\sigma$ values. It is formatted as described below:

```
<arbitraty name for species 1> <sigma_11> <sigma_12> <sigma_13> … <sigma_1n>
<arbitraty name for species 2> <sigma_21> <sigma_22> <sigma_23> … <sigma_2n>

...

<arbitrary name for species n> <sigma_n1> <sigma_n2> <sigma_n3> … <sigma_nn>
```
