CXX=g++ -std=c++0x -std=gnu++0x
##Normal flags
CFLAGS=-c -fopenmp #-O2
CFLAGSmain= -fopenmp# -O2
##Debugging flags
#CFLAGS=-Wall -Wextra -Wcast-qual -Wcast-align -O0 -ggdb -g3 -fstack-protector-all -fno-inline -c -fopenmp #-O2
#CFLAGSmain=-Wall -Wextra -Wcast-qual -Wcast-align -O0 -ggdb -g3 -fstack-protector -fno-inline -all -fopenmp #-O2
SERVER=
WAVEVECTORS3D="\"/home/david/qvectors/qvectors3d/qvector\""
WAVEVECTORS2D="\"/home/david/qvectors/qvectors2d/qvector\""
WAVEVECTORS1D="\"/home/david/qvectors/qvectors1d/qvector\""

ifeq ($(SERVER),TACC)
AMDAT: amdat.o coordinate.o trajectory.o atom_trajectory.o molecule.o system.o analysis.o space-time_correlation_function.o van_hove_self.o progress.o mean_square_displacement.o van_hove_distinct.o control.o wave_vectors.o wave_density.o intermediate_scattering_function.o correlation_2d.o incoherent_scattering_function.o debyewaller_dist.o stiffness_dist.o gaussian_comparison.o fast_particles.o non_gaussian_parameter.o  gaussian_comparison.o radial_debye_waller.o tokenize.o mean_square_displacement_2d.o velocity_autocorrelation.o strings.o trajectory_list.o static_trajectory_list.o exptime_trajectory_list.o rgtensor.o trajmath.o rgtensor_stats.o displacement_distribution.o version.h boolean_list.o fast_particles.o   displacement_map.o composition.o n_fold_order_parameter.o trajectory_list_bins.o structure_factor.o clustered_list.o trajectory_list_decay.o vector_autocorrelation.o error.o mean_displacement.o multibody.o multibody_set.o
	${CXX} $(CFLAGSmain) amdat.o tokenize.o coordinate.o trajectory.o atom_trajectory.o molecule.o system.o analysis.o space-time_correlation_function.o van_hove_self.o progress.o mean_square_displacement.o van_hove_distinct.o control.o wave_vectors.o  wave_density.o intermediate_scattering_function.o correlation_2d.o incoherent_scattering_function.o debyewaller_dist.o stiffness_dist.o non_gaussian_parameter.o gaussian_comparison.o radial_debye_waller.o mean_square_displacement_2d.o velocity_autocorrelation.o strings.o trajectory_list.o  static_trajectory_list.o clustered_list.o exptime_trajectory_list.o rgtensor.o trajmath.o rgtensor_stats.o displacement_distribution.o boolean_list.o fast_particles.o   composition.o n_fold_order_parameter.o displacement_map.o trajectory_list_bins.o structure_factor.o vector_autocorrelation.o trajectory_list_decay.o error.o mean_displacement.o multibody.o multibody_set.o -o AMDAT 

system.o: system.cpp system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h
	${CXX} $(CFLAGS) system.cpp -DTACC

space-time_correlation_function.o: space-time_correlation_function.cpp space-time_correlation_function.h system.h coordinate.h molecule.h analysis.h atom_trajectory.h trajectory.h
	${CXX} $(CFLAGS) space-time_correlation_function.cpp -DTACC

velocity_autocorrelation.o: velocity_autocorrelation.cpp velocity_autocorrelation.h mean_square_displacement.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h trajectory.h
	${CXX} $(CFLAGS) velocity_autocorrelation.cpp -DTACC

control.o: control.cpp control.h system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h debyewaller_dist.h stiffness_dist.h non_gaussian_parameter.h gaussian_comparison.h fast_particles.h   tokenize.h radial_debye_waller.h mean_square_displacement_2d.h velocity_autocorrelation.h strings.h rgtensor_stats.h displacement_map.h trajectory_list_bins.h bin_dynamics_analysis.h bin_static_analysis.h composition.h n_fold_order_parameter.h trajectory_list_decay.h multibody_set.h multibody.h
	${CXX} $(CFLAGS) control.cpp -DTACC
else
AMDAT: amdat.o coordinate.o trajectory.o atom_trajectory.o molecule.o system.o analysis.o space-time_correlation_function.o van_hove_self.o progress.o mean_square_displacement.o van_hove_distinct.o control.o wave_vectors.o wave_density.o intermediate_scattering_function.o correlation_2d.o incoherent_scattering_function.o debyewaller_dist.o stiffness_dist.o gaussian_comparison.o fast_particles.o non_gaussian_parameter.o  gaussian_comparison.o radial_debye_waller.o tokenize.o mean_square_displacement_2d.o velocity_autocorrelation.o strings.o trajectory_list.o static_trajectory_list.o clustered_list.o exptime_trajectory_list.o rgtensor.o trajmath.o rgtensor_stats.o displacement_distribution.o version.h boolean_list.o fast_particles.o   displacement_map.o composition.o n_fold_order_parameter.o trajectory_list_bins.o structure_factor.o trajectory_list_decay.o vector_autocorrelation.o error.o mean_displacement.o multibody.o multibody_set.o
	${CXX} $(CFLAGSmain) amdat.o tokenize.o coordinate.o trajectory.o atom_trajectory.o molecule.o system.o analysis.o space-time_correlation_function.o van_hove_self.o progress.o mean_square_displacement.o van_hove_distinct.o control.o wave_vectors.o wave_density.o intermediate_scattering_function.o correlation_2d.o incoherent_scattering_function.o debyewaller_dist.o stiffness_dist.o non_gaussian_parameter.o gaussian_comparison.o radial_debye_waller.o mean_square_displacement_2d.o velocity_autocorrelation.o strings.o trajectory_list.o  static_trajectory_list.o clustered_list.o exptime_trajectory_list.o rgtensor.o trajmath.o rgtensor_stats.o displacement_distribution.o boolean_list.o fast_particles.o   composition.o n_fold_order_parameter.o displacement_map.o trajectory_list_bins.o structure_factor.o vector_autocorrelation.o trajectory_list_decay.o error.o mean_displacement.o multibody.o multibody_set.o -o AMDAT -lfftw3 -lm -lxdrfile 

system.o: system.cpp system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h
	${CXX} $(CFLAGS) system.cpp -lxdrfile

space-time_correlation_function.o: space-time_correlation_function.cpp space-time_correlation_function.h system.h coordinate.h molecule.h analysis.h atom_trajectory.h trajectory.h
	${CXX} $(CFLAGS) space-time_correlation_function.cpp -lfftw3 -lm 

velocity_autocorrelation.o: velocity_autocorrelation.cpp velocity_autocorrelation.h mean_square_displacement.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h trajectory.h
	${CXX} $(CFLAGS) velocity_autocorrelation.cpp -lfftw3 -lm

control.o: control.cpp control.h system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h debyewaller_dist.h stiffness_dist.h non_gaussian_parameter.h gaussian_comparison.h fast_particles.h   tokenize.h radial_debye_waller.h mean_square_displacement_2d.h velocity_autocorrelation.h strings.h rgtensor_stats.h displacement_map.h trajectory_list_bins.h bin_dynamics_analysis.h bin_static_analysis.h composition.h n_fold_order_parameter.h trajectory_list_decay.h multibody_set.h multibody.h
	${CXX} $(CFLAGS) control.cpp
endif

amdat.o: amdat.cpp system.h van_hove_self.h mean_square_displacement.h van_hove_distinct.h molecule.h atom_trajectory.h coordinate.h analysis.h control.h
	${CXX} $(CFLAGS) amdat.cpp

coordinate.o: coordinate.h coordinate.cpp
	${CXX} $(CFLAGS) coordinate.cpp
	
trajectory.o: trajectory.cpp trajectory.h coordinate.h
	${CXX} $(CFLAGS) trajectory.cpp

atom_trajectory.o: atom_trajectory.cpp atom_trajectory.h trajectory.h coordinate.h
	${CXX} $(CFLAGS) atom_trajectory.cpp
	
molecule.o: molecule.cpp molecule.h atom_trajectory.h coordinate.h trajectory.h multibody.h
	${CXX} $(CFLAGS) molecule.cpp
	
analysis.o: analysis.cpp system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h
	${CXX} $(CFLAGS) analysis.cpp
	
van_hove_self.o: van_hove_self.cpp van_hove_self.h space-time_correlation_function.h system.h coordinate.h molecule.h atom_trajectory.h analysis.h trajectory.h trajectory_list.h
	${CXX} $(CFLAGS) van_hove_self.cpp
	
progress.o: progress.cpp progress.h
	${CXX} $(CFLAGS) progress.cpp
	
mean_square_displacement.o: mean_square_displacement.cpp mean_square_displacement.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h  trajectory.h trajectory_list.h
	${CXX} $(CFLAGS) mean_square_displacement.cpp

mean_displacement.o: mean_displacement.cpp mean_displacement.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h  trajectory.h trajectory_list.h
	${CXX} $(CFLAGS) mean_displacement.cpp
	
van_hove_distinct.o: van_hove_distinct.cpp van_hove_distinct.h space-time_correlation_function.h system.h coordinate.h molecule.h analysis.h atom_trajectory.h trajectory.h
	${CXX} $(CFLAGS) van_hove_distinct.cpp

wave_vectors.o: wave_vectors.h wave_vectors.cpp coordinate.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h trajectory.h
	${CXX} $(CFLAGS) wave_vectors.cpp -DWV3D=${WAVEVECTORS3D} -DWV2D=${WAVEVECTORS2D} -DWV1D=${WAVEVECTORS1D}
	
wave_density.o: wave_density.h wave_density.cpp system.h wave_vectors.h molecule.h analysis.h atom_trajectory.h coordinate.h trajectory.h
	${CXX} $(CFLAGS) wave_density.cpp
	
intermediate_scattering_function.o: intermediate_scattering_function.h intermediate_scattering_function.cpp wave_density.h wave_vectors.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h correlation_2d.h trajectory.h
	${CXX} $(CFLAGS) intermediate_scattering_function.cpp

correlation_2d.o: correlation_2d.h correlation_2d.cpp analysis.h system.h atom_trajectory.h coordinate.h molecule.h  trajectory.h
	${CXX} $(CFLAGS) correlation_2d.cpp

incoherent_scattering_function.o: incoherent_scattering_function.h incoherent_scattering_function.cpp correlation_2d.h wave_vectors.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h
	${CXX} $(CFLAGS) incoherent_scattering_function.cpp

debyewaller_dist.o: debyewaller_dist.h debyewaller_dist.cpp analysis.h coordinate.h molecule.h system.h atom_trajectory.h trajectory.h
	${CXX} $(CFLAGS) debyewaller_dist.cpp

stiffness_dist.o: stiffness_dist.h stiffness_dist.cpp analysis.h coordinate.h molecule.h system.h atom_trajectory.h trajectory.h
	${CXX} $(CFLAGS) stiffness_dist.cpp

displacement_distribution.o: displacement_distribution.h displacement_distribution.cpp analysis.h coordinate.h molecule.h system.h atom_trajectory.h trajectory.h
	${CXX} $(CFLAGS) displacement_distribution.cpp

non_gaussian_parameter.o: non_gaussian_parameter.h non_gaussian_parameter.cpp system.h atom_trajectory.h coordinate.h molecule.h analysis.h mean_square_displacement.h trajectory.h
	${CXX} $(CFLAGS) non_gaussian_parameter.cpp

gaussian_comparison.o: gaussian_comparison.h gaussian_comparison.cpp system.h molecule.h analysis.h atom_trajectory.h coordinate.h non_gaussian_parameter.h mean_square_displacement.h van_hove_self.h space-time_correlation_function.h trajectory.h
	${CXX} $(CFLAGS) gaussian_comparison.cpp

radial_debye_waller.o: radial_debye_waller.h radial_debye_waller.cpp analysis.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h
	${CXX} $(CFLAGS) radial_debye_waller.cpp
	
tokenize.o: tokenize.h tokenize.cpp 
	${CXX} $(CFLAGS) tokenize.cpp
	
mean_square_displacement_2d.o: mean_square_displacement_2d.cpp mean_square_displacement_2d.h system.h molecule.h atom_trajectory.h coordinate.h analysis.h  trajectory.h
	${CXX} $(CFLAGS) mean_square_displacement_2d.cpp

strings.o: strings.h strings.cpp analysis.h atom_trajectory.h coordinate.h molecule.h trajectory.h progress.h system.h
	${CXX} $(CFLAGS) strings.cpp
	
trajectory_list.o: trajectory_list.h trajectory_list.cpp analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h
	${CXX} $(CFLAGS) trajectory_list.cpp
	
static_trajectory_list.o: static_trajectory_list.cpp static_trajectory_list.h trajectory_list.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h
	${CXX} $(CFLAGS) static_trajectory_list.cpp
	
exptime_trajectory_list.o: exptime_trajectory_list.cpp exptime_trajectory_list.h trajectory_list.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h 
	${CXX} $(CFLAGS) exptime_trajectory_list.cpp

rgtensor.o: rgtensor.cpp rgtensor.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h
	${CXX} $(CFLAGS) rgtensor.cpp

trajmath.o: trajmath.h trajmath.cpp coordinate.h 
	${CXX} $(CFLAGS) trajmath.cpp 

rgtensor_stats.o: rgtensor_stats.cpp rgtensor_stats.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h
	${CXX} $(CFLAGS) rgtensor_stats.cpp 

boolean_list.o: boolean_list.cpp boolean_list.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h
	${CXX} $(CFLAGS) boolean_list.cpp 

fast_particles.o: fast_particles.h fast_particles.cpp analysis.h atom_trajectory.h coordinate.h molecule.h system.h gaussian_comparison.h  trajectory.h exptime_trajectory_list.h trajectory_list.h 
	${CXX} $(CFLAGS) fast_particles.cpp

n_fold_order_parameter.o: n_fold_order_parameter.h n_fold_order_parameter.cpp version.h value_list.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h boolean_list.h
	${CXX} $(CFLAGS) n_fold_order_parameter.cpp

composition.o: composition.h composition.cpp analysis.h version.h system.h molecule.h atom_trajectory.h coordinate.h tokenize.h trajectory.h trajectory_list.h
	${CXX} $(CFLAGS) composition.cpp

displacement_map.o: displacement_map.h displacement_map.cpp version.h value_list.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h boolean_list.h
	${CXX} $(CFLAGS) displacement_map.cpp

trajectory_list_bins.o: trajectory_list_bins.cpp trajectory_list_bins.h system.h molecule.h analysis.h atom_trajectory.h coordinate.h tokenize.h trajectory.h trajectory_list.h boolean_list.h progress.h
	${CXX} $(CFLAGS) trajectory_list_bins.cpp

structure_factor.o: structure_factor.cpp structure_factor.h system.h analysis.h atom_trajectory.h coordinate.h molecule.h  trajectory.h trajectory_list.h wave_vectors.h version.h tokenize.h
	${CXX} $(CFLAGS) structure_factor.cpp

clustered_list.o: clustered_list.cpp clustered_list.h trajectory_list.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h trajectory.h
	${CXX} $(CFLAGS) clustered_list.cpp

vector_autocorrelation.o: vector_autocorrelation.cpp vector_autocorrelation.h analysis.h system.h atom_trajectory.h coordinate.h molecule.h analysis.h trajectory.h trajectory_list.h version.h tokenize.h
	${CXX} $(CFLAGS) vector_autocorrelation.cpp
	
trajectory_list_decay.o: trajectory_list_decay.cpp trajectory_list_decay.h analysis.h system.h trajectory_list.h version.h boolean_list.h trajectory.h molecule.h coordinate.h atom_trajectory.h
	${CXX} $(CFLAGS) trajectory_list_decay.cpp

error.o: error.cpp error.h control.cpp control.h
	${CXX} $(CFLAGS) error.cpp

multibody.o: multibody.cpp multibody.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h 
	${CXX} $(CFLAGS) multibody.cpp

multibody_set.o: multibody_set.cpp multibody_set.h trajectory.h coordinate.h system.h analysis.h atom_trajectory.h tokenize.h 
	${CXX} $(CFLAGS) multibody_set.cpp

clean: 
	rm -f amdat.o coordinate.o trajectory.o atom_trajectory.o molecule.o system.o analysis.o space-time_correlation_function.o van_hove_self.o progress.o mean_square_displacement.o van_hove_distinct.o control.o wave_vectors.o wave_density.o intermediate_scattering_function.o correlation_2d.o incoherent_scattering_function.o debyewaller_dist.o non_gaussian_parameter.o gaussian_comparison.o gaussian_comparison.o fast_particles.o radial_debye_waller.o tokenize.o mean_square_displacement_2d.o velocity_autocorrelation.o strings.o trajectory_list.o static_trajectory_list.o exptime_trajectory_list.o rgtensor.o trajmath.o rgtensor_stats.o displacement_distribution.o boolean_list.o value_list.o displacement_map.o trajectory_list_bins.o n_fold_order_parameter.o composition.o structure_factor.o clustered_list.o vector_autocorrelation.o trajectory_list_decay.o stiffness_dist.o error.o mean_displacement.o multibody.o multibody_set.o
