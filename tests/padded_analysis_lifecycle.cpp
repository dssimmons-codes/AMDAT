#include "bond_autocorrelation_function.h"
#include "mean_displacement.h"
#include "mean_square_displacement.h"
#include "mean_square_displacement_2d.h"
#include "space_time_correlation_function.h"

#include <cassert>
#include <cmath>

using namespace std;

void print_progress(int, int)
{
}

class Test_System : public System
{
  public:
    Test_System(int exponential_steps, int exponentials)
    {
      n_exponential_steps = exponential_steps;
      n_exponentials = exponentials;
      frt = false;
      n_timesteps = n_exponential_steps*n_exponentials+1;
      n_timegaps = n_exponential_steps+n_exponentials;
      displacement_limit = 0;
      rho = new float[1];
      rho[0] = 1.0f;

      timelist = new float[n_timesteps];
      for(int timeii=0; timeii<n_timesteps; ++timeii)
      {
        timelist[timeii] = static_cast<float>(timeii);
      }
    }

    ~Test_System()
    {
      delete [] timelist;
      delete [] rho;
    }
};

class Test_Space_Time_Correlation : public Space_Time_Correlation_Function
{
  public:
    explicit Test_Space_Time_Correlation(System * test_system)
    {
      system = test_system;
      n_bins = 4;
      n_times = 2;
      bin_size = 1.0f;
      max_value = 4.0f;

      correlation = new float*[n_times];
      weighting = new int*[n_times];
      timetable = new float[n_times];
      for(int timeii=0;timeii<n_times;++timeii)
      {
        correlation[timeii] = new float[n_bins];
        weighting[timeii] = new int[16]{};
        timetable[timeii] = static_cast<float>(timeii);
        for(int binii=0;binii<n_bins;++binii)
        {
          correlation[timeii][binii] = static_cast<float>(binii+1);
        }
      }
    }

    ~Test_Space_Time_Correlation()
    {
      clear_memory();
    }

    float inverse_value(int timeii, int wavenumberii) const
    {
      return spatial_inverse[timeii][wavenumberii];
    }
};

int main()
{
  Test_System small_system(1, 2);
  Test_System large_system(3, 4);
  Coordinate xyz(1, 1, 1);

  Test_Space_Time_Correlation space_time(&small_system);
  for(int iteration=0;iteration<100;++iteration)
  {
    const int n_wavenumbers = iteration%2 == 0 ? 8 : 13;
    space_time.calculate_spatial_inverse(n_wavenumbers);
    assert(std::isfinite(space_time.inverse_value(0, n_wavenumbers-1)));
  }

  for(int iteration=0; iteration<100; ++iteration)
  {
    Mean_Square_Displacement msd(&small_system);
    Mean_Square_Displacement msd_copy(msd);
    Mean_Square_Displacement msd_assigned;
    msd_assigned = msd;
    msd = msd;
    msd.set(&large_system);
    assert(msd.show(0) == 0.0f);
    assert(msd_copy.show(0) == 0.0f);
    assert(msd_assigned.show(0) == 0.0f);

    Mean_Displacement md(&small_system);
    Mean_Displacement md_copy(md);
    Mean_Displacement md_assigned;
    md_assigned = md;
    md = md;
    md.set(&large_system);
    assert(md.show(0) == Coordinate(0, 0, 0));
    assert(md_copy.show(0) == Coordinate(0, 0, 0));
    assert(md_assigned.show(0) == Coordinate(0, 0, 0));

    Mean_Square_Displacement_2D msd2d_small(&small_system, "xy");
    Mean_Square_Displacement_2D msd2d_large(&large_system, "xz");
    Mean_Square_Displacement_2D msd2d_copy(msd2d_small);
    msd2d_large = msd2d_small;
    msd2d_small = msd2d_small;
    assert(msd2d_large.show(0) == 0.0f);
    assert(msd2d_copy.show(0) == 0.0f);

    Bond_Autocorrelation_Function baf_system_only(&small_system);
    Bond_Autocorrelation_Function baf_large(&large_system, xyz);
    Bond_Autocorrelation_Function baf_copy(baf_system_only);
    baf_large = baf_system_only;
    baf_system_only = baf_system_only;
    baf_system_only.set(&large_system);
  }

  return 0;
}
