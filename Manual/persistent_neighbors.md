<h1>Persistent Neighbors</h1>

Create a new multibody list (and new multibodies thereing) comprised of clusters of particles that remain neighbors after a specified time spacing. Generates a distinct multibody_list for each time block. Syntax is as follows.

_persistent\_neighbors \<time spacing at which to performan analysis\> \<name neighbor list to use\> \<name of new multibody list and new trajectory list to create\> <"centroid" or "COM">_

This should be followed by a trajectory_list containing the particles to be analyzed.
