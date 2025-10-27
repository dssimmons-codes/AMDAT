<h1>value list command</h1>

Versatile command that can do several things based on the keyword that follows. The general structure is as follows.

``value_list <keyword> <arguments>``

A list of keywords and associated arguments is as follows.

``value_list threshold_value <value list name> <new trajectory list to create> <threshold keyword> <threshold1> <optional: threshold 2>``

Constructs a new trajectory list, by selecting only those particles with values in the selected value list that are above, below, or between selected threshold values.
``<threshold keyword>`` may be ``greater``, ``less``, or ``between``. When ``greater``, ``<threshold1>`` sets the value above which particles are selected and ``<threshold2>`` is not used. When ``less``, ``<threshold1>`` sets the value below which particles are selected and ``<threshold2>`` is not used. When ``between``, ``<threshold1>`` sets the value above which particles are selected and ``<threshold2>`` sets the value below which particles are selected.

``value_list threshold_percentile <value list name> <new trajectory list to create> <threshold keyword> <percentile threshold1> <optional: percentile threshold 2>``

Constructs a new trajectory list, by selecting only those particles whose values in the selected value list are in a specified percentile range (among all values in the list), either above, below, or between selected threshold percentile values.
``<threshold keyword>`` may be ``greater``, ``less``, or ``between``. When ``greater``, ``<threshold1>`` sets the percentile value above which particles are selected and ``<threshold2>`` is not used. When ``less``, ``<threshold1>`` sets the percentile value below which particles are selected and ``<threshold2>`` is not used. When ``between``, ``<threshold1>`` sets the value above which particles are selected and ``<threshold2>`` sets the percentile value below which particles are selected.

``value_list write_pdb <value list name> <output file name stem> <value time index> <position time index> <optional: maximum value to write out>``

Writes out a .pdb file in which the value in the beta column is drawn from a value list. This can readily allow visualization of the system with particles colored by a computed per-particle value in ``VMD`` or ``OVITO``. ``<output file name stem>`` is the filename prior to the .pdb suffix, which will be added by AMDAT before writing. ``<value time index>`` specifies the time at which the value will be drawn from the value list. ``<position time index>`` specifies the time at which the configuration will be written. The optional maximum value argument will set a largest value to be written to the beta value in the pdb file; any value higher than this will be reduced to this value. This can help in visualization work where a small number of high values can 'blow out' the color scale.
