<h1>isf_list</h1>

Calculate intermediate scattering function using a trajectory list.

```
isf_list <output file> <list name> <plane> <first wavenumber index> <last wavenumber index>
```

Options for `<plane>` are `xyz`, `xy`, `xz`, `yz`, `x`, `y`, and `z`. This chooses which dimensions in $k$-space to include in the calculation of the intermediate scattering function. `xyz` computes the full radial three dimensional isf, `xy`, `yz`, and `xz` calculate two-dimensional in-plane radial isf's, and `x`, `y`, and `z` compute one-dimensional isf's.
