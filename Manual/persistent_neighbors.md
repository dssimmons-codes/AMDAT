<h1>Persistent Neighbors</h1>

Create a new multibody list (and new multibodies therein) comprised of clusters of particles that remain neighbors after a specified time spacing. Generates a distinct ``multibody_list`` for each time block. Syntax is as follows.

```
persistent_neighbors <time spacing at which to perform analysis> <name of neighbor list to use> <name of new multibody list and new trajectory list to create> <"centroid" or "COM">
```

This should be followed by a ``trajectory_list`` containing the particles to be analyzed.
