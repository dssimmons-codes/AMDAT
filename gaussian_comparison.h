/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Class to compare self Van Hove function to the so-called Gaussian approximation for the Van Hove, establishing boundaries and fractions for slow and fast particle categories*/
/*Written by David S. Simmons*/

#ifndef GAUSSIAN_COMPARISON
#define GAUSSIAN_COMPARISON

#include "system.h"
#include "non_gaussian_parameter.h"
#include "van_hove_self.h"

namespace std{
		
class Gaussian_Comparison
{
    Non_Gaussian_Parameter const * non_gaussian_parameter;
    Van_Hove_Self const * self_van_hove;
    System const * system;
    
    float mean_square_displacement;
    float * van_hove;				//stores van hove times shell volume
    float * gaussian_approx;			//stores gaussian approx for van hove times shell volume
    int n_bins;
    float bin_size;
    int time_index;				//stores t*, the time corresponding to the peak in the non-gaussian parameter
    int slowboundary_index;
    int fastboundary_index;
    float slowboundary;			//maximum mean-square-displacement of slow particles
    float fastboundary;			//minimum mean-square-displacement of fast particles
    float slowfraction;
    float fastfraction;
    float gaussian_slowfraction;
    float gaussian_fastfraction;
    
    bool errorstate;
    
    void scale_van_hoves();
    void find_crossovers();
    void determine_fractions();
    
    void clear();
  public:
    Gaussian_Comparison();
    Gaussian_Comparison(System* sys, const Non_Gaussian_Parameter* ngp, const Van_Hove_Self* vhs, const Mean_Square_Displacement* msd);
    Gaussian_Comparison(System* sys, int time_in, const Van_Hove_Self* vhs, const Mean_Square_Displacement* msd);
    void set(System* sys, const Non_Gaussian_Parameter* ngp, const Van_Hove_Self* vhs, const Mean_Square_Displacement* msd);
    void set(System* sys, int time_in, const Van_Hove_Self* vhs, const Mean_Square_Displacement* msd);
    void write(string)const;
    
    float show_slowboundary()const{return slowboundary;};
    float show_fastboundary()const{return fastboundary;};
    int show_time_index()const{return time_index;};
};

}

#endif