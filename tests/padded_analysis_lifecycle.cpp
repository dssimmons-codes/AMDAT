#include "bond_autocorrelation_function.h"
#include "mean_displacement.h"
#include "mean_square_displacement.h"
#include "mean_square_displacement_2d.h"

#include <cassert>

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

      timelist = new float[n_timesteps];
      for(int timeii=0; timeii<n_timesteps; ++timeii)
      {
        timelist[timeii] = static_cast<float>(timeii);
      }
    }

    ~Test_System()
    {
      delete [] timelist;
    }
};

int main()
{
  Test_System small_system(1, 2);
  Test_System large_system(3, 4);
  Coordinate xyz(1, 1, 1);

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
