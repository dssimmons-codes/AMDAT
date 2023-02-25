CXX=g++ -std=gnu++0x -std=c++0x 
##Normal flags
#CFLAGS=-c -fopenmp #-O2
#CFLAGSmain= -fopenmp# -O2
##Debugging flags
CFLAGS=-Wall -Wextra -Wcast-qual -Wcast-align -O0 -ggdb -g3 -fstack-protector-all -fno-inline -c -fopenmp #-O2
CFLAGSmain=-Wall -Wextra -Wcast-qual -Wcast-align -O0 -ggdb -g3 -fstack-protector-all -fno-inline -fopenmp #-O2
SERVER=TACC
WAVEVECTORS3D="\"./src/qvectors/qvectors3d/qvector\""
WAVEVECTORS2D="\"./src/qvectors/qvectors2d/qvector\""
WAVEVECTORS1D="\"./src/qvectors/qvectors1d/qvector\""

OBJDIR = ./build

OBJECTS= $(OBJDIR)/amdat.o $(OBJDIR)/tokenize.o $(OBJDIR)/coordinate.o $(OBJDIR)/trajectory.o $(OBJDIR)/atom_trajectory.o $(OBJDIR)/molecule.o $(OBJDIR)/system.o $(OBJDIR)/analysis.o $(OBJDIR)/space-time_correlation_function.o $(OBJDIR)/van_hove_self.o $(OBJDIR)/progress.o $(OBJDIR)/mean_square_displacement.o $(OBJDIR)/van_hove_distinct.o $(OBJDIR)/control.o \
$(OBJDIR)/wave_vectors.o $(OBJDIR)/wave_density.o $(OBJDIR)/intermediate_scattering_function.o $(OBJDIR)/correlation_2d.o $(OBJDIR)/incoherent_scattering_function.o $(OBJDIR)/debyewaller_dist.o $(OBJDIR)/stiffness_dist.o $(OBJDIR)/non_gaussian_parameter.o \
$(OBJDIR)/gaussian_comparison.o $(OBJDIR)/radial_debye_waller.o $(OBJDIR)/mean_square_displacement_2d.o $(OBJDIR)/velocity_autocorrelation.o $(OBJDIR)/strings.o $(OBJDIR)/trajectory_list.o $(OBJDIR)/static_trajectory_list.o $(OBJDIR)/exptime_trajectory_list.o $(OBJDIR)/rgtensor.o $(OBJDIR)/trajmath.o $(OBJDIR)/rgtensor_stats.o \
$(OBJDIR)/displacement_distribution.o $(OBJDIR)/boolean_list.o $(OBJDIR)/fast_particles.o $(OBJDIR)/displacement_map.o $(OBJDIR)/composition.o $(OBJDIR)/n_fold_order_parameter.o $(OBJDIR)/trajectory_list_bins.o $(OBJDIR)/structure_factor.o $(OBJDIR)/clustered_list.o $(OBJDIR)/trajectory_list_decay.o \
$(OBJDIR)/vector_autocorrelation.o $(OBJDIR)/error.o $(OBJDIR)/mean_displacement.o $(OBJDIR)/multibody.o $(OBJDIR)/multibody_set.o $(OBJDIR)/multibody_list.o $(OBJDIR)/multibody_analysis.o $(OBJDIR)/gyration_radius.o $(OBJDIR)/trajectory_set.o $(OBJDIR)/edge_detector_timedependent.o $(OBJDIR)/mean_velocity_unsteady.o \
$(OBJDIR)/mean_unsteady_displacement.o $(OBJDIR)/analysis_onetime.o $(OBJDIR)/radial_distribution_function.o $(OBJDIR)/bond_autocorrelation_function.o $(OBJDIR)/displacement_list.o $(OBJDIR)/orientational_correlation.o $(OBJDIR)/size_statistics.o $(OBJDIR)/multibody_region.o $(OBJDIR)/provisional_multibodies.o \
$(OBJDIR)/dynamic_cluster_multibodies.o $(OBJDIR)/string_multibodies.o $(OBJDIR)/comover_multibodies.o $(OBJDIR)/relative_displacement_strings.o $(OBJDIR)/neighbor_list.o $(OBJDIR)/distance_neighbor_list.o $(OBJDIR)/persistent_neighbors.o $(OBJDIR)/voronoi_neighbor_list.o $(OBJDIR)/neighbor_decorrelation_function.o $(OBJDIR)/radial_count.o $(OBJDIR)/mean_closest_distance.o $(OBJDIR)/particles_between.o

CONTROLHEADERS=src/control.h src/system.h src/van_hove_self.h src/mean_square_displacement.h src/van_hove_distinct.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/debyewaller_dist.h src/stiffness_dist.h src/non_gaussian_parameter.h \
src/gaussian_comparison.h src/fast_particles.h src/tokenize.h src/radial_debye_waller.h src/mean_square_displacement_2d.h src/velocity_autocorrelation.h src/strings.h src/rgtensor_stats.h src/displacement_map.h src/trajectory_list_bins.h src/bin_dynamics_analysis.h \
src/bin_static_analysis.h src/composition.h src/n_fold_order_parameter.h src/trajectory_list_decay.h src/multibody_set.h src/multibody.h src/multibody_list.h src/multibody_analysis.h src/gyration_radius.h src/trajectory_set.h src/edge_detector_timedependent.h \
src/mean_velocity_unsteady.h src/mean_unsteady_displacement.h src/radial_distribution_function.h src/bond_autocorrelation_function.h src/orientational_correlation.h src/size_statistics.h src/provisional_multibodies.h src/dynamic_cluster_multibodies.h \
src/string_multibodies.h src/comover_multibodies.h src/relative_displacement_strings.h src/neighbor_list.h src/distance_neighbor_list.h src/persistent_neighbors.h src/voronoi_neighbor_list.h src/neighbor_decorrelation_function.h src/radial_count.h src/mean_closest_distance.h src/particles_between.h

## Voronoi needs to be compiled at ./src/voro++-0.4.6/ level before compiling AMDAT
VPATH = ./src/voro++-0.4.6/src

ifeq ($(SERVER),TACC)
AMDAT: $(OBJECTS) src/version.h
	$(CXX) $(CFLAGSmain) $(OBJECTS) -o AMDAT -I./src/voro++-0.4.6/src -L./src/voro++-0.4.6/src -lvoro++

$(OBJDIR)/system.o: src/system.cpp src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/trajectory_set.h
	$(CXX) $(CFLAGS) src/system.cpp -DTACC -o $@

$(OBJDIR)/velocity_autocorrelation.o: src/velocity_autocorrelation.cpp src/velocity_autocorrelation.h src/mean_square_displacement.h src/system.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/trajectory.h
	$(CXX) $(CFLAGS) src/velocity_autocorrelation.cpp -DTACC -o $@

$(OBJDIR)/control.o: src/control.cpp $(CONTROLHEADERS)
	$(CXX) $(CFLAGS) src/control.cpp -DTACC -I./src/voro++-0.4.6/src -o $@
else
AMDAT: $(OBJECTS) src/version.h
	$(CXX) $(CFLAGSmain) $(OBJECTS) -o AMDAT -lfftw3 -lm -I./src/voro++-0.4.6/src -L./src/voro++-0.4.6/src -lvoro++ -o $@

$(OBJDIR)/system.o: src/system.cpp src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/trajectory_set.h
	$(CXX) $(CFLAGS) src/system.cpp -o $@

$(OBJDIR)/velocity_autocorrelation.o: src/velocity_autocorrelation.cpp src/velocity_autocorrelation.h src/mean_square_displacement.h src/system.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/trajectory.h
	$(CXX) $(CFLAGS) src/velocity_autocorrelation.cpp -lfftw3 -lm -o $@

$(OBJDIR)/control.o: src/control.cpp $(CONTROLHEADERS)
	$(CXX) $(CFLAGS) src/control.cpp -I./src/voro++-0.4.6/src -o $@
endif

$(OBJDIR)/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/amdat.o: src/amdat.cpp src/system.h src/van_hove_self.h src/mean_square_displacement.h src/van_hove_distinct.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/control.h

$(OBJDIR)/coordinate.o: src/coordinate.h src/coordinate.cpp

$(OBJDIR)/trajectory.o: src/trajectory.cpp src/trajectory.h src/coordinate.h

$(OBJDIR)/atom_trajectory.o: src/atom_trajectory.cpp src/atom_trajectory.h src/trajectory.h src/coordinate.h

$(OBJDIR)/molecule.o: src/molecule.cpp src/molecule.h src/atom_trajectory.h src/coordinate.h src/trajectory.h src/multibody.h

$(OBJDIR)/analysis.o: src/analysis.cpp src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h

$(OBJDIR)/analysis_onetime.o: src/analysis_onetime.cpp src/analysis_onetime.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h

$(OBJDIR)/van_hove_self.o: src/van_hove_self.cpp src/van_hove_self.h src/space-time_correlation_function.h src/system.h src/coordinate.h src/molecule.h src/atom_trajectory.h src/analysis.h src/trajectory.h src/trajectory_list.h

$(OBJDIR)/progress.o: src/progress.cpp src/progress.h

$(OBJDIR)/mean_square_displacement.o: src/mean_square_displacement.cpp src/mean_square_displacement.h src/system.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/trajectory.h src/trajectory_list.h

$(OBJDIR)/mean_displacement.o: src/mean_displacement.cpp src/mean_displacement.h src/system.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/trajectory.h src/trajectory_list.h

#msd_listprint.o: src/msd_listprint.cpp src/msd_listprint.h src/system.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/trajectory.h src/trajectory_list.h

$(OBJDIR)/van_hove_distinct.o: src/van_hove_distinct.cpp src/van_hove_distinct.h src/space-time_correlation_function.h src/system.h src/coordinate.h src/molecule.h src/analysis.h src/atom_trajectory.h src/trajectory.h

$(OBJDIR)/wave_vectors.o: src/wave_vectors.h src/wave_vectors.cpp src/coordinate.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/trajectory.h
	$(CXX) $(CFLAGS) src/wave_vectors.cpp -DWV3D=${WAVEVECTORS3D} -DWV2D=${WAVEVECTORS2D} -DWV1D=${WAVEVECTORS1D} -o $@

$(OBJDIR)/wave_density.o: src/wave_density.h src/wave_density.cpp src/system.h src/wave_vectors.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/trajectory.h

$(OBJDIR)/intermediate_scattering_function.o: src/intermediate_scattering_function.h src/intermediate_scattering_function.cpp src/wave_density.h src/wave_vectors.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/correlation_2d.h src/trajectory.h

$(OBJDIR)/correlation_2d.o: src/correlation_2d.h src/correlation_2d.cpp src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h

$(OBJDIR)/incoherent_scattering_function.o: src/incoherent_scattering_function.h src/incoherent_scattering_function.cpp src/correlation_2d.h src/wave_vectors.h src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h

$(OBJDIR)/debyewaller_dist.o: src/debyewaller_dist.h src/debyewaller_dist.cpp src/analysis.h src/coordinate.h src/molecule.h src/system.h src/atom_trajectory.h src/trajectory.h

$(OBJDIR)/stiffness_dist.o: src/stiffness_dist.h src/stiffness_dist.cpp src/analysis.h src/coordinate.h src/molecule.h src/system.h src/atom_trajectory.h src/trajectory.h

$(OBJDIR)/displacement_distribution.o: src/displacement_distribution.h src/displacement_distribution.cpp src/analysis.h src/coordinate.h src/molecule.h src/system.h src/atom_trajectory.h src/trajectory.h

$(OBJDIR)/non_gaussian_parameter.o: src/non_gaussian_parameter.h src/non_gaussian_parameter.cpp src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/mean_square_displacement.h src/trajectory.h

$(OBJDIR)/gaussian_comparison.o: src/gaussian_comparison.h src/gaussian_comparison.cpp src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/non_gaussian_parameter.h src/mean_square_displacement.h src/van_hove_self.h src/space-time_correlation_function.h src/trajectory.h

$(OBJDIR)/radial_debye_waller.o: src/radial_debye_waller.h src/radial_debye_waller.cpp src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h

$(OBJDIR)/tokenize.o: src/tokenize.h src/tokenize.cpp 

$(OBJDIR)/mean_square_displacement_2d.o: src/mean_square_displacement_2d.cpp src/mean_square_displacement_2d.h src/system.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/trajectory.h

$(OBJDIR)/space-time_correlation_function.o: src/space-time_correlation_function.cpp src/space-time_correlation_function.h src/system.h src/coordinate.h src/molecule.h src/analysis.h src/atom_trajectory.h src/trajectory.h

$(OBJDIR)/strings.o: src/strings.h src/strings.cpp src/analysis.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h src/progress.h src/system.h

$(OBJDIR)/trajectory_list.o: src/trajectory_list.h src/trajectory_list.cpp src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h

$(OBJDIR)/static_trajectory_list.o: src/static_trajectory_list.cpp src/static_trajectory_list.h src/trajectory_list.h src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h src/trajectory_set.h src/multibody_set.h

$(OBJDIR)/exptime_trajectory_list.o: src/exptime_trajectory_list.cpp src/exptime_trajectory_list.h src/trajectory_list.h src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h 

$(OBJDIR)/rgtensor.o: src/rgtensor.cpp src/rgtensor.h src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h

$(OBJDIR)/trajmath.o: src/trajmath.h src/trajmath.cpp src/coordinate.h 

$(OBJDIR)/rgtensor_stats.o: src/rgtensor_stats.cpp src/rgtensor_stats.h src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h

$(OBJDIR)/boolean_list.o: src/boolean_list.cpp src/boolean_list.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h

$(OBJDIR)/fast_particles.o: src/fast_particles.h src/fast_particles.cpp src/analysis.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/system.h src/gaussian_comparison.h src/trajectory.h src/exptime_trajectory_list.h src/trajectory_list.h 

$(OBJDIR)/n_fold_order_parameter.o: src/n_fold_order_parameter.h src/n_fold_order_parameter.cpp src/version.h src/value_list.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h src/boolean_list.h

$(OBJDIR)/composition.o: src/composition.h src/composition.cpp src/analysis.h src/version.h src/system.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/trajectory_list.h

$(OBJDIR)/displacement_map.o: src/displacement_map.h src/displacement_map.cpp src/version.h src/value_list.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h src/boolean_list.h

$(OBJDIR)/trajectory_list_bins.o: src/trajectory_list_bins.cpp src/trajectory_list_bins.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/trajectory_list.h src/boolean_list.h src/progress.h

$(OBJDIR)/structure_factor.o: src/structure_factor.cpp src/structure_factor.h src/system.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h src/trajectory_list.h src/wave_vectors.h src/version.h src/tokenize.h

$(OBJDIR)/clustered_list.o: src/clustered_list.cpp src/clustered_list.h src/trajectory_list.h src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h

$(OBJDIR)/vector_autocorrelation.o: src/vector_autocorrelation.cpp src/vector_autocorrelation.h src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h src/version.h src/tokenize.h

$(OBJDIR)/trajectory_list_decay.o: src/trajectory_list_decay.cpp src/trajectory_list_decay.h src/analysis.h src/system.h src/trajectory_list.h src/version.h src/boolean_list.h src/trajectory.h src/molecule.h src/coordinate.h src/atom_trajectory.h

$(OBJDIR)/error.o: src/error.cpp src/error.h src/control.cpp src/control.h

$(OBJDIR)/multibody.o: src/multibody.cpp src/multibody.h src/trajectory.h src/coordinate.h src/system.h src/analysis.h src/atom_trajectory.h src/tokenize.h 

$(OBJDIR)/multibody_set.o: src/multibody_set.cpp src/multibody_set.h src/trajectory.h src/coordinate.h src/system.h src/analysis.h src/atom_trajectory.h src/tokenize.h src/trajectory_set.h

$(OBJDIR)/multibody_list.o: src/multibody_list.cpp src/multibody_list.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/multibody.h src/multibody_set.h src/multibody_analysis.h src/multibody.h

$(OBJDIR)/multibody_analysis.o: src/multibody_analysis.cpp src/multibody_analysis.h src/multibody_list.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/multibody.h src/multibody_set.h src/version.h

$(OBJDIR)/trajectory_set.o: src/trajectory_set.cpp src/trajectory_set.h src/trajectory.h src/coordinate.h src/system.h src/analysis.h src/atom_trajectory.h src/tokenize.h 

$(OBJDIR)/gyration_radius.o: src/gyration_radius.cpp src/gyration_radius.h src/multibody_analysis.h src/multibody_list.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/multibody.h src/multibody_set.h src/version.h

$(OBJDIR)/bond_autocorrelation_function.o: src/bond_autocorrelation_function.cpp src/bond_autocorrelation_function.h src/multibody_analysis.h src/multibody_list.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/multibody.h src/multibody_set.h src/version.h

$(OBJDIR)/edge_detector_timedependent.o: src/edge_detector_timedependent.cpp src/edge_detector_timedependent.h src/system.h src/coordinate.h src/analysis.h src/version.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/trajectory_set.h

$(OBJDIR)/mean_velocity_unsteady.o: src/mean_velocity_unsteady.cpp src/mean_velocity_unsteady.h src/system.h src/coordinate.h src/version.h src/static_trajectory_list.h src/analysis.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/trajectory_set.h

$(OBJDIR)/mean_unsteady_displacement.o: src/mean_unsteady_displacement.cpp src/mean_unsteady_displacement.h src/system.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/trajectory.h src/trajectory_list.h

$(OBJDIR)/radial_distribution_function.o: src/radial_distribution_function.cpp src/radial_distribution_function.h src/analysis_onetime.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h src/version.h

$(OBJDIR)/particles_between.o: src/particles_between.cpp src/particles_between.h src/analysis_onetime.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h src/version.h

$(OBJDIR)/mean_closest_distance.o: src/mean_closest_distance.cpp src/mean_closest_distance.h src/analysis_onetime.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h src/version.h

$(OBJDIR)/radial_count.o: src/radial_count.cpp src/radial_count.h src/analysis_onetime.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/analysis.h src/trajectory.h src/trajectory_list.h src/version.h

$(OBJDIR)/displacement_list.o: src/displacement_list.cpp src/displacement_list.h src/static_trajectory_list.h src/value_list.h src/analysis.h src/boolean_list.h src/trajectory_list.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/trajectory_set.h

$(OBJDIR)/orientational_correlation.o: src/orientational_correlation.cpp src/orientational_correlation.h src/multibody_analysis.h src/multibody_list.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/multibody.h src/multibody_set.h src/version.h

$(OBJDIR)/size_statistics.o: src/size_statistics.cpp src/size_statistics.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/multibody.h src/multibody_set.h src/multibody_analysis.h src/version.h

$(OBJDIR)/multibody_region.o: src/multibody_region.cpp src/multibody_region.h src/system.h src/molecule.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/tokenize.h src/trajectory.h src/multibody.h src/multibody_set.h src/multibody_analysis.h src/version.h

$(OBJDIR)/provisional_multibodies.o: src/provisional_multibodies.cpp src/provisional_multibodies.h src/multibody_set.h src/trajectory.h src/coordinate.h src/system.h src/analysis.h src/atom_trajectory.h src/tokenize.h src/trajectory_set.h src/multibody.h src/control.h src/system.h src/van_hove_self.h src/mean_square_displacement.h src/van_hove_distinct.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/debyewaller_dist.h src/stiffness_dist.h src/non_gaussian_parameter.h src/gaussian_comparison.h src/fast_particles.h src/tokenize.h src/radial_debye_waller.h src/mean_square_displacement_2d.h src/velocity_autocorrelation.h src/strings.h src/rgtensor_stats.h src/displacement_map.h src/trajectory_list_bins.h src/bin_dynamics_analysis.h src/bin_static_analysis.h src/composition.h src/n_fold_order_parameter.h src/trajectory_list_decay.h src/multibody_set.h src/multibody.h src/multibody_list.h src/multibody_analysis.h src/gyration_radius.h src/trajectory_set.h src/edge_detector_timedependent.h src/mean_velocity_unsteady.h src/mean_unsteady_displacement.h\
 src/radial_distribution_function.h src/bond_autocorrelation_function.h src/orientational_correlation.h src/multibody_region.h src/coordinate.h src/size_statistics.h

$(OBJDIR)/dynamic_cluster_multibodies.o: src/dynamic_cluster_multibodies.h src/dynamic_cluster_multibodies.cpp src/provisional_multibodies.h src/provisional_multibodies.h src/multibody_set.h src/trajectory.h src/coordinate.h src/system.h src/analysis.h src/atom_trajectory.h src/tokenize.h src/trajectory_set.h src/multibody.h src/control.h src/system.h src/van_hove_self.h src/mean_square_displacement.h src/van_hove_distinct.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/debyewaller_dist.h src/stiffness_dist.h src/non_gaussian_parameter.h src/gaussian_comparison.h src/fast_particles.h src/tokenize.h src/radial_debye_waller.h src/mean_square_displacement_2d.h src/velocity_autocorrelation.h src/strings.h src/rgtensor_stats.h src/displacement_map.h src/trajectory_list_bins.h src/bin_dynamics_analysis.h src/bin_static_analysis.h src/composition.h src/n_fold_order_parameter.h src/trajectory_list_decay.h src/multibody_set.h src/multibody.h src/multibody_list.h src/multibody_analysis.h src/gyration_radius.h src/trajectory_set.h src/edge_detector_timedependent.h src/mean_velocity_unsteady.h src/mean_unsteady_displacement.h\
 src/radial_distribution_function.h src/radial_count.h src/bond_autocorrelation_function.h src/orientational_correlation.h src/multibody_region.h src/coordinate.h src/size_statistics.h

$(OBJDIR)/string_multibodies.o: src/string_multibodies.h src/string_multibodies.cpp src/dynamic_cluster_multibodies.h src/provisional_multibodies.h src/provisional_multibodies.h src/multibody_set.h src/trajectory.h src/coordinate.h src/system.h src/analysis.h src/atom_trajectory.h src/tokenize.h src/trajectory_set.h src/multibody.h src/control.h src/system.h src/van_hove_self.h src/mean_square_displacement.h src/van_hove_distinct.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/debyewaller_dist.h src/stiffness_dist.h src/non_gaussian_parameter.h src/gaussian_comparison.h src/fast_particles.h src/tokenize.h src/radial_debye_waller.h src/mean_square_displacement_2d.h src/velocity_autocorrelation.h src/strings.h src/rgtensor_stats.h src/displacement_map.h src/trajectory_list_bins.h src/bin_dynamics_analysis.h src/bin_static_analysis.h src/composition.h src/n_fold_order_parameter.h src/trajectory_list_decay.h src/multibody_set.h src/multibody.h src/multibody_list.h src/multibody_analysis.h src/gyration_radius.h src/trajectory_set.h src/edge_detector_timedependent.h src/mean_velocity_unsteady.h src/mean_unsteady_displacement.h\
 src/radial_distribution_function.h src/bond_autocorrelation_function.h src/orientational_correlation.h src/multibody_region.h src/coordinate.h src/size_statistics.h src/radial_count.h

$(OBJDIR)/comover_multibodies.o: src/comover_multibodies.h src/comover_multibodies.cpp src/dynamic_cluster_multibodies.h src/provisional_multibodies.h src/provisional_multibodies.h src/multibody_set.h src/trajectory.h src/coordinate.h src/system.h src/analysis.h src/atom_trajectory.h src/tokenize.h src/trajectory_set.h src/multibody.h src/control.h src/system.h src/van_hove_self.h src/mean_square_displacement.h src/van_hove_distinct.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/debyewaller_dist.h src/stiffness_dist.h src/non_gaussian_parameter.h src/gaussian_comparison.h src/fast_particles.h src/tokenize.h src/radial_debye_waller.h src/mean_square_displacement_2d.h src/velocity_autocorrelation.h src/strings.h src/rgtensor_stats.h src/displacement_map.h src/trajectory_list_bins.h src/bin_dynamics_analysis.h src/bin_static_analysis.h src/composition.h src/n_fold_order_parameter.h src/trajectory_list_decay.h src/multibody_set.h src/multibody.h src/multibody_list.h src/multibody_analysis.h src/gyration_radius.h src/trajectory_set.h src/edge_detector_timedependent.h src/mean_velocity_unsteady.h src/mean_unsteady_displacement.h\
 src/radial_distribution_function.h src/bond_autocorrelation_function.h src/orientational_correlation.h src/multibody_region.h src/coordinate.h src/size_statistics.h src/radial_count.h

$(OBJDIR)/relative_displacement_strings.o: src/relative_displacement_strings.h src/relative_displacement_strings.cpp src/dynamic_cluster_multibodies.h src/provisional_multibodies.h src/provisional_multibodies.h src/multibody_set.h src/trajectory.h src/coordinate.h src/system.h src/analysis.h src/atom_trajectory.h src/tokenize.h src/trajectory_set.h src/multibody.h src/control.h src/system.h src/van_hove_self.h src/mean_square_displacement.h src/van_hove_distinct.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/debyewaller_dist.h src/stiffness_dist.h src/non_gaussian_parameter.h src/gaussian_comparison.h src/fast_particles.h src/tokenize.h src/radial_debye_waller.h src/mean_square_displacement_2d.h src/velocity_autocorrelation.h src/strings.h src/rgtensor_stats.h src/displacement_map.h src/trajectory_list_bins.h src/bin_dynamics_analysis.h src/bin_static_analysis.h src/composition.h src/n_fold_order_parameter.h src/trajectory_list_decay.h src/multibody_set.h src/multibody.h src/multibody_list.h src/multibody_analysis.h src/gyration_radius.h src/trajectory_set.h src/edge_detector_timedependent.h src/mean_velocity_unsteady.h src/mean_unsteady_displacement.h\
 src/radial_distribution_function.h src/bond_autocorrelation_function.h src/orientational_correlation.h src/multibody_region.h src/coordinate.h src/size_statistics.h src/radial_count.h

$(OBJDIR)/neighbor_list.o: src/neighbor_list.cpp src/neighbor_list.h src/trajectory_list.h src/system.h src/analysis.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h src/value_list.h

$(OBJDIR)/distance_neighbor_list.o: src/distance_neighbor_list.cpp src/distance_neighbor_list.h src/neighbor_list.h src/trajectory_list.h src/system.h src/trajectory_list.h src/trajectory_list.cpp src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h src/analysis_onetime.h src/coordinate.h src/molecule.h src/value_list.h

$(OBJDIR)/voronoi_neighbor_list.o: src/voronoi_neighbor_list.cpp src/voronoi_neighbor_list.h src/neighbor_list.h src/trajectory_list.h src/system.h src/trajectory_list.h src/trajectory_list.cpp src/analysis.h src/system.h src/atom_trajectory.h src/coordinate.h src/molecule.h src/trajectory.h src/analysis_onetime.h src/coordinate.h src/molecule.h src/value_list.h
	$(CXX) $(CFLAGS) src/voronoi_neighbor_list.cpp -I./src/voro++-0.4.6/src -o $@

$(OBJDIR)/persistent_neighbors.o: src/persistent_neighbors.cpp src/persistent_neighbors.h src/dynamic_cluster_multibodies.h src/provisional_multibodies.h src/provisional_multibodies.h src/multibody_set.h src/trajectory.h src/coordinate.h src/system.h src/analysis.h src/atom_trajectory.h src/tokenize.h src/trajectory_set.h src/multibody.h src/control.h src/system.h src/van_hove_self.h src/mean_square_displacement.h src/van_hove_distinct.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/debyewaller_dist.h src/stiffness_dist.h src/non_gaussian_parameter.h src/gaussian_comparison.h src/fast_particles.h src/tokenize.h src/radial_debye_waller.h src/mean_square_displacement_2d.h src/velocity_autocorrelation.h src/strings.h src/rgtensor_stats.h src/displacement_map.h src/trajectory_list_bins.h src/bin_dynamics_analysis.h src/bin_static_analysis.h src/composition.h src/n_fold_order_parameter.h src/trajectory_list_decay.h src/multibody_set.h src/multibody.h src/multibody_list.h src/multibody_analysis.h src/gyration_radius.h src/trajectory_set.h src/edge_detector_timedependent.h src/mean_velocity_unsteady.h src/mean_unsteady_displacement.h\
 src/radial_distribution_function.h src/bond_autocorrelation_function.h src/orientational_correlation.h src/multibody_region.h src/coordinate.h src/size_statistics.h src/neighbor_list.h src/radial_count.h

$(OBJDIR)/neighbor_decorrelation_function.o: src/neighbor_decorrelation_function.cpp src/neighbor_decorrelation_function.h src/system.h src/molecule.h src/atom_trajectory.h src/coordinate.h src/analysis.h src/trajectory.h src/trajectory_list.h src/neighbor_list.h src/value_list.h

clean: 
	rm -f $(OBJECTS)
