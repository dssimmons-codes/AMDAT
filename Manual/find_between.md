<h1>find_between</h1>

`find_between` creates a trajectory of particles (from the first target) that reside between two other particles (from the second target) subject to constraints on the total distance and angle formed between the two vectors.

```
find_between <name of trajectory list to create> <dist_cutoff> <costheta_cutoff> <is_only_diff>
<target1>
<target2>
```

<h3 align="center">2D illustration of find_between command with two particles separated by a distance of 3</a></h3>

<table>
  <tr>
    <td align="center" width="25%">
      <img src="assets/commands/find_between/fig-find_between_schematic_dc4.0_tc0.0_eg1.png" alt="Two-dimensional illustration of how AMDAT's find_between command works." width="100%">
      <br/>
      <sub><code>dist_cutoff</code> = 4.0;<br/><code>costheta_cutoff</code> = 0.0;<br/> Vectors and angle for a sample point</sub>
    </td>
    <td align="center" width="25%">
      <img src="assets/commands/find_between/fig-find_between_schematic_dc4.0_tc-0.5_eg0.png" alt="Two-dimensional illustration of how AMDAT's find_between command works." width="100%">
      <br/>
      <sub><code>dist_cutoff</code> = 4.0;<br/><code>costheta_cutoff</code> = -0.5</sub>
    </td>
    <td align="center" width="25%">
      <img src="assets/commands/find_between/fig-find_between_schematic_dc5.0_tc0.0_eg0.png" alt="Two-dimensional illustration of how AMDAT's find_between command works." width="100%">
      <br/>
      <sub><code>dist_cutoff</code> = 5.0;<br/><code>costheta_cutoff</code> = 0.0</sub>
    </td>
    <td align="center" width="25%">
      <img src="assets/commands/find_between/fig-find_between_schematic_dc5.0_tc0.5_eg0.png" alt="Two-dimensional illustration of how AMDAT's find_between command works." width="100%">
      <br/>
      <sub><code>dist_cutoff</code> = 5.0;<br/><code>costheta_cutoff</code> = 0.5</sub>
    </td>
  </tr>
</table>

Developed for Kawak, Bhapkar, Simmons. _Macromolecules_ 57 (2024). doi: <a href="https://pubs.acs.org/doi/10.1021/acs.macromol.4c00489#fig1">10.1021/acs.macromol.4c00489

As can be seen by the illustration above, AMDAT will compute two vectors with an atom from the first target (`<target1>`) as the center with two atoms from the second target (`<target2>`).
An atom from the first target will be added to the resulting trajectory list if the two vector magnitudes ($d_1$ and $d_2$) and the angle formed between them obey ($\theta$) the following criteria:

$$
\frac{d_1+d_2}{2} < d_c
$$
$$
\cos \theta < \cos \theta_c
$$

where $d_c$ is `<dist_cutoff>` and $\cos \theta_c$ is `<costheta_cutoff>`. The final option in this command is `<is_only_diff>` and it ensures that the two atoms from the second target do not belong to the same molecule. This last option helps distinguish inter-molecular phenomena from intra-molecular phenomena.
