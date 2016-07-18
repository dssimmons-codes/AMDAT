CXX=g++ -std=gnu++0x -std=c++0x 
##Normal flags
CFLAGS=-c -fopenmp #-O2
CFLAGSmain= -fopenmp# -O2
##Debugging flags
#CFLAGS=-Wall -Wextra -Wcast-qual -Wcast-align -O0 -ggdb -g3 -fstack-protector-all -fno-inline -c -fopenmp #-O2
#CFLAGSmain=-Wall -Wextra -Wcast-qual -Wcast-align -O0 -ggdb -g3 -fstack-protector-all -fno-inline -fopenmp #-O2
SERVER=TACC
WAVEVECTORS3D="\"./qvectors/qvectors3d/qvector\""
WAVEVECTORS2D="\"./qvectors/qvectors2d/qvector\""
WAVEVECTORS1D="\"./qvectors/qvectors1d/qvector\""

OBJECTS= amdat.o tokenize.o coordinate.o trajectory.o atom_trajectory.o molecule.o system.o analysis.o space-time_correlation_function.o van_hove_self.o progress.o mean_square_displacement.o van_hove_distinct.o control.o \
wave_vectors.o wave_density.o intermediate_scattering_function.o correlation_2d.o incoherent_scattering_function.o debyewaller_dist.o stiffness_dist.o non_gaussian_parameter.o \
gaussian_comparison.o radial_debye_waller.o mean_square_displacement_2d.o velocity_autocorrelation.o strings.o trajectory_list.o static_trajectory_list.o exptime_trajectory_list.o rgtensor.o trajmath.o rgtensor_stats.o \
displacement_distribution.o boolean_list.o fast_particles.o   displacement_map.o composition.o n_fold_order_parameter.o trajectory_list_bins.o structure_factor.o clustered_list.o trajectory_list_decay.o \
vector_autocorrelation.o error.o mean_displacement.o multibody.o multibody_set.o multibody_list.o multibody_analysis.o gyration_radius.o trajectory_set.o edge_detector_timedependent.o mean_velocity_unsteady.o \
mean_unsteady_displacement.o analysis_onetime.o radial_distribution_function.o bond_autocorrelation_function.o displacement_list.o orientational_correlation.o size_statistics.o multibody_region.o provisional_multibodies.o \
dynamic_cluster_multibodies.o string_multibodies.o comover_multibodies.o relative_displacement_strings.o neighbor_list.o distance_neighbor_list.o persistent_neighbors.o voronoi_neighbor_list.o

CONTROLHEADERS=control.h system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h debyewaller_dist.h stiffness_dist.h non_gaussian_parameter.h \
gaussian_comparison.h fast_particles.h tokenize.h radial_debye_waller.h mean_square_displacement_2d.h velocity_autocorrelation.h strings.h rgtensor_stats.h displacement_map.h trajectory_list_bins.h bin_dynamics_analysis.h \
bin_static_analysis.h composition.h n_fold_order_parameter.h trajectory_list_decay.h multibody_set.h multibody.h multibody_list.h multibody_analysis.h gyration_radius.h trajectory_set.h edge_detector_timedependent.h \
mean_velocity_unsteady.h mean_unsteady_displacement.h radial_distribution_function.h bond_autocorrelation_function.h orientational_correlation.h size_statistics.h provisional_multibodies.h dynamic_cluster_multibodies.h \
string_multibodies.h comover_multibodies.h relative_displacement_strings.h neighbor_list.h distance_neighbor_list.h persistent_neighbors.h voronoi_neighbor_list.h

## Voronoi needs to be compiled at ./voro++-0.4.6/ level before compiling AMDAT
VPATH = ./voro++-0.4.6/src

ifeq ($(SERVER),TACC)
AMDAT: $(OBJECTS) version.h
	$(CXX) $(CFLAGSmain) $(OBJECTS) -o AMDAT -I./voro++-0.4.6/src -L./voro++-0.4.6/src -lvoro++

system.o: system.cpp system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h trajectory_set.h
	$(CXX) $(CFLAGS) system.cpp -DTACC

velocity_autocorrelation.o: velocity_autocorrelation.cpp velocity_autocorrelation.h mean_square_displacement.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h trajectory.h
	$(CXX) $(CFLAGS) velocity_autocorrelation.cpp -DTACC

control.o: control.cpp  $(CONTROLHEADERS)
	$(CXX) $(CFLAGS) control.cpp -DTACC -I./voro++-0.4.6/src
else
AMDAT: $(OBJECTS) version.h
	$(CXX) $(CFLAGSmain) $(OBJECTS) -o AMDAT -lfftw3 -lm -I./voro++-0.4.6/src -L./voro++-0.4.6/src -lvoro++

system.o: system.cpp system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h trajectory_set.h
	$(CXX) $(CFLAGS) system.cpp

velocity_autocorrelation.o: velocity_autocorrelation.cpp velocity_autocorrelation.h mean_square_displacement.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h trajectory.h
	$(CXX) $(CFLAGS) velocity_autocorrelation.cpp -lfftw3 -lm

control.o: control.cpp $(CONTROLHEADERS)
	$(CXX) $(CFLAGS) control.cpp -I./voro++-0.4.6/src
endif

%.o: %.cpp
	$(CXX) $(CFLAGS) $< -o $@

amdat.o: amdat.cpp system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h control.h

coordinate.o: coordinate.h coordinate.cpp

trajectory.o: trajectory.cpp trajectory.h coordinate.h

atom_trajectory.o: atom_trajectory.cpp atom_trajectory.h trajectory.h coordinate.h

molecule.o: molecule.cpp molecule.h atom_trajectory.h coordinate.h trajectory.h multibody.h

analysis.o: analysis.cpp system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h

analysis_onetime.o: analysis_onetime.cpp analysis_onetime.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h

van_hove_self.o: van_hove_self.cpp van_hove_self.h space-time_correlation_function.h system.h coordinate.h molecule.h atom_trajectory.h analysis.h trajectory.h trajectory_list.h

progress.o: progress.cpp progress.h

mean_square_displacement.o: mean_square_displacement.cpp mean_square_displacement.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h  trajectory.h trajectory_list.h

mean_displacement.o: mean_displacement.cpp mean_displacement.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h  trajectory.h trajectory_list.h

van_hove_distinct.o: van_hove_distinct.cpp van_hove_distinct.h space-time_correlation_function.h system.h coordinate.h molecule.h analysis.h atom_trajectory.h trajectory.h

wave_vectors.o: wave_vectors.h wave_vectors.cpp coordinate.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h trajectory.h
	$(CXX) $(CFLAGS) wave_vectors.cpp -DWV3D=${WAVEVECTORS3D} -DWV2D=${WAVEVECTORS2D} -DWV1D=${WAVEVECTORS1D}

wave_density.o: wave_density.h wave_density.cpp system.h wave_vectors.h molecule.h analysis.h atom_trajectory.h coordinate.h trajectory.h

intermediate_scattering_function.o: intermediate_scattering_function.h intermediate_scattering_function.cpp wave_density.h wave_vectors.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h correlation_2d.h trajectory.h

correlation_2d.o: correlation_2d.h correlation_2d.cpp analysis.h system.h atom_trajectory.h coordinate.h molecule.h  trajectory.h

incoherent_scattering_function.o: incoherent_scattering_function.h incoherent_scattering_function.cpp correlation_2d.h wave_vectors.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h

debyewaller_dist.o: debyewaller_dist.h debyewaller_dist.cpp analysis.h coordinate.h molecule.h system.h atom_trajectory.h trajectory.h

stiffness_dist.o: stiffness_dist.h stiffness_dist.cpp analysis.h coordinate.h molecule.h system.h atom_trajectory.h trajectory.h

displacement_distribution.o: displacement_distribution.h displacement_distribution.cpp analysis.h coordinate.h molecule.h system.h atom_trajectory.h trajectory.h

non_gaussian_parameter.o: non_gaussian_parameter.h non_gaussian_parameter.cpp system.h atom_trajectory.h coordinate.h molecule.h analysis.h mean_square_displacement.h trajectory.h

gaussian_comparison.o: gaussian_comparison.h gaussian_comparison.cpp system.h molecule.h analysis.h atom_trajectory.h coordinate.h non_gaussian_parameter.h mean_square_displacement.h van_hove_self.h space-time_correlation_function.h trajectory.h

radial_debye_waller.o: radial_debye_waller.h radial_debye_waller.cpp analysis.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h

tokenize.o: tokenize.h tokenize.cpp 

mean_square_displacement_2d.o: mean_square_displacement_2d.cpp mean_square_displacement_2d.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h  trajectory.h

space-time_correlation_function.o: space-time_correlation_function.cpp space-time_correlation_function.h system.h coordinate.h molecule.h analysis.h atom_trajectory.h trajectory.h

strings.o: strings.h strings.cpp analysis.h atom_trajectory.h coordinate.h molecule.h trajectory.h progress.h system.h

trajectory_list.o: trajectory_list.h trajectory_list.cpp analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h

static_trajectory_list.o: static_trajectory_list.cpp static_trajectory_list.h trajectory_list.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h trajectory_set.h multibody_set.h

exptime_trajectory_list.o: exptime_trajectory_list.cpp exptime_trajectory_list.h trajectory_list.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h 

rgtensor.o: rgtensor.cpp rgtensor.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h

trajmath.o: trajmath.h trajmath.cpp coordinate.h 

rgtensor_stats.o: rgtensor_stats.cpp rgtensor_stats.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h

boolean_list.o: boolean_list.cpp boolean_list.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h

fast_particles.o: fast_particles.h fast_particles.cpp analysis.h atom_trajectory.h coordinate.h molecule.h system.h gaussian_comparison.h  trajectory.h exptime_trajectory_list.h trajectory_list.h 

n_fold_order_parameter.o: n_fold_order_parameter.h n_fold_order_parameter.cpp version.h value_list.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h boolean_list.h

composition.o: composition.h composition.cpp analysis.h version.h system.h molecule.h atom_trajectory.h coordinate.h tokenize.h trajectory.h trajectory_list.h

displacement_map.o: displacement_map.h displacement_map.cpp version.h value_list.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h boolean_list.h

trajectory_list_bins.o: trajectory_list_bins.cpp trajectory_list_bins.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h trajectory_list.h boolean_list.h progress.h

structure_factor.o: structure_factor.cpp structure_factor.h system.h analysis.h atom_trajectory.h coordinate.h molecule.h  trajectory.h trajectory_list.h wave_vectors.h version.h tokenize.h

clustered_list.o: clustered_list.cpp clustered_list.h trajectory_list.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h

vector_autocorrelation.o: vector_autocorrelation.cpp vector_autocorrelation.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h version.h tokenize.h

trajectory_list_decay.o: trajectory_list_decay.cpp trajectory_list_decay.h analysis.h system.h trajectory_list.h version.h boolean_list.h trajectory.h molecule.h coordinate.h atom_trajectory.h

error.o: error.cpp error.h control.cpp control.h

multibody.o: multibody.cpp multibody.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h 

multibody_set.o: multibody_set.cpp multibody_set.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h trajectory_set.h

multibody_list.o: multibody_list.cpp multibody_list.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h multibody.h multibody_set.h multibody_analysis.h multibody.h

multibody_analysis.o: multibody_analysis.cpp multibody_analysis.h multibody_list.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h multibody.h multibody_set.h version.h

trajectory_set.o: trajectory_set.cpp trajectory_set.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h 

gyration_radius.o: gyration_radius.cpp gyration_radius.h multibody_analysis.h multibody_list.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h multibody.h multibody_set.h version.h

bond_autocorrelation_function.o: bond_autocorrelation_function.cpp bond_autocorrelation_function.h multibody_analysis.h multibody_list.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h multibody.h multibody_set.h version.h

edge_detector_timedependent.o: edge_detector_timedependent.cpp edge_detector_timedependent.h system.h coordinate.h analysis.h version.h molecule.h  atom_trajectory.h coordinate.h tokenize.h trajectory.h trajectory_set.h

mean_velocity_unsteady.o: mean_velocity_unsteady.cpp mean_velocity_unsteady.h system.h coordinate.h version.h static_trajectory_list.h analysis.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h trajectory_set.h

mean_unsteady_displacement.o: mean_unsteady_displacement.cpp mean_unsteady_displacement.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h  trajectory.h trajectory_list.h

radial_distribution_function.o: radial_distribution_function.cpp radial_distribution_function.h analysis_onetime.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h version.h

displacement_list.o: displacement_list.cpp displacement_list.h static_trajectory_list.h value_list.h analysis.h boolean_list.h trajectory_list.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h trajectory_set.h

orientational_correlation.o: orientational_correlation.cpp orientational_correlation.h multibody_analysis.h multibody_list.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h multibody.h multibody_set.h version.h

size_statistics.o:size_statistics.cpp size_statistics.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h multibody.h multibody_set.h multibody_analysis.h version.h

multibody_region.o: multibody_region.cpp multibody_region.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h multibody.h multibody_set.h multibody_analysis.h version.h

provisional_multibodies.o: provisional_multibodies.cpp provisional_multibodies.h multibody_set.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h trajectory_set.h multibody.h control.h system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h debyewaller_dist.h stiffness_dist.h non_gaussian_parameter.h gaussian_comparison.h fast_particles.h   tokenize.h radial_debye_waller.h mean_square_displacement_2d.h velocity_autocorrelation.h strings.h rgtensor_stats.h displacement_map.h trajectory_list_bins.h bin_dynamics_analysis.h bin_static_analysis.h composition.h n_fold_order_parameter.h trajectory_list_decay.h multibody_set.h multibody.h multibody_list.h multibody_analysis.h gyration_radius.h trajectory_set.h edge_detector_timedependent.h mean_velocity_unsteady.h mean_unsteady_displacement.h\
 radial_distribution_function.h bond_autocorrelation_function.h orientational_correlation.h multibody_region.h coordinate.h size_statistics.h

dynamic_cluster_multibodies.o: dynamic_cluster_multibodies.h dynamic_cluster_multibodies.cpp provisional_multibodies.h provisional_multibodies.h multibody_set.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h trajectory_set.h multibody.h control.h system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h debyewaller_dist.h stiffness_dist.h non_gaussian_parameter.h gaussian_comparison.h fast_particles.h   tokenize.h radial_debye_waller.h mean_square_displacement_2d.h velocity_autocorrelation.h strings.h rgtensor_stats.h displacement_map.h trajectory_list_bins.h bin_dynamics_analysis.h bin_static_analysis.h composition.h n_fold_order_parameter.h trajectory_list_decay.h multibody_set.h multibody.h multibody_list.h multibody_analysis.h gyration_radius.h trajectory_set.h edge_detector_timedependent.h mean_velocity_unsteady.h mean_unsteady_displacement.h\
 radial_distribution_function.h bond_autocorrelation_function.h orientational_correlation.h multibody_region.h coordinate.h size_statistics.h

string_multibodies.o: string_multibodies.h string_multibodies.cpp dynamic_cluster_multibodies.h provisional_multibodies.h provisional_multibodies.h multibody_set.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h trajectory_set.h multibody.h control.h system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h debyewaller_dist.h stiffness_dist.h non_gaussian_parameter.h gaussian_comparison.h fast_particles.h   tokenize.h radial_debye_waller.h mean_square_displacement_2d.h velocity_autocorrelation.h strings.h rgtensor_stats.h displacement_map.h trajectory_list_bins.h bin_dynamics_analysis.h bin_static_analysis.h composition.h n_fold_order_parameter.h trajectory_list_decay.h multibody_set.h multibody.h multibody_list.h multibody_analysis.h gyration_radius.h trajectory_set.h edge_detector_timedependent.h mean_velocity_unsteady.h mean_unsteady_displacement.h\
 radial_distribution_function.h bond_autocorrelation_function.h orientational_correlation.h multibody_region.h coordinate.h size_statistics.h

comover_multibodies.o: comover_multibodies.h comover_multibodies.cpp dynamic_cluster_multibodies.h provisional_multibodies.h provisional_multibodies.h multibody_set.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h trajectory_set.h multibody.h control.h system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h debyewaller_dist.h stiffness_dist.h non_gaussian_parameter.h gaussian_comparison.h fast_particles.h   tokenize.h radial_debye_waller.h mean_square_displacement_2d.h velocity_autocorrelation.h strings.h rgtensor_stats.h displacement_map.h trajectory_list_bins.h bin_dynamics_analysis.h bin_static_analysis.h composition.h n_fold_order_parameter.h trajectory_list_decay.h multibody_set.h multibody.h multibody_list.h multibody_analysis.h gyration_radius.h trajectory_set.h edge_detector_timedependent.h mean_velocity_unsteady.h mean_unsteady_displacement.h\
 radial_distribution_function.h bond_autocorrelation_function.h orientational_correlation.h multibody_region.h coordinate.h size_statistics.h

relative_displacement_strings.o: relative_displacement_strings.h relative_displacement_strings.cpp dynamic_cluster_multibodies.h provisional_multibodies.h provisional_multibodies.h multibody_set.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h trajectory_set.h multibody.h control.h system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h debyewaller_dist.h stiffness_dist.h non_gaussian_parameter.h gaussian_comparison.h fast_particles.h   tokenize.h radial_debye_waller.h mean_square_displacement_2d.h velocity_autocorrelation.h strings.h rgtensor_stats.h displacement_map.h trajectory_list_bins.h bin_dynamics_analysis.h bin_static_analysis.h composition.h n_fold_order_parameter.h trajectory_list_decay.h multibody_set.h multibody.h multibody_list.h multibody_analysis.h gyration_radius.h trajectory_set.h edge_detector_timedependent.h mean_velocity_unsteady.h mean_unsteady_displacement.h\
 radial_distribution_function.h bond_autocorrelation_function.h orientational_correlation.h multibody_region.h coordinate.h size_statistics.h

neighbor_list.o: neighbor_list.cpp neighbor_list.h trajectory_list.h system.h trajectory_list.h trajectory_list.cpp analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h value_list.h

distance_neighbor_list.o: distance_neighbor_list.cpp distance_neighbor_list.h neighbor_list.h trajectory_list.h system.h trajectory_list.h trajectory_list.cpp analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h analysis_onetime.h coordinate.h molecule.h value_list.h

voronoi_neighbor_list.o: voronoi_neighbor_list.cpp voronoi_neighbor_list.h neighbor_list.h trajectory_list.h system.h trajectory_list.h trajectory_list.cpp analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h analysis_onetime.h coordinate.h molecule.h value_list.h
	$(CXX) $(CFLAGS) voronoi_neighbor_list.cpp -I./voro++-0.4.6/src

persistent_neighbors.o:	persistent_neighbors.cpp persistent_neighbors.h dynamic_cluster_multibodies.h provisional_multibodies.h provisional_multibodies.h multibody_set.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h trajectory_set.h multibody.h control.h system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h debyewaller_dist.h stiffness_dist.h non_gaussian_parameter.h gaussian_comparison.h fast_particles.h   tokenize.h radial_debye_waller.h mean_square_displacement_2d.h velocity_autocorrelation.h strings.h rgtensor_stats.h displacement_map.h trajectory_list_bins.h bin_dynamics_analysis.h bin_static_analysis.h composition.h n_fold_order_parameter.h trajectory_list_decay.h multibody_set.h multibody.h multibody_list.h multibody_analysis.h gyration_radius.h trajectory_set.h edge_detector_timedependent.h mean_velocity_unsteady.h mean_unsteady_displacement.h\
 radial_distribution_function.h bond_autocorrelation_function.h orientational_correlation.h multibody_region.h coordinate.h size_statistics.h neighbor_list.h


clean: 
	rm -f $(OBJECTS)