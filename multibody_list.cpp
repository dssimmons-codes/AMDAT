/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Methods for Multibody_List class*/
/*Written by David S. Simmons*/

#include <fstream>
#include <stdlib.h>
#include "multibody_list.h"
#include "analysis.h"
#include "multibody_analysis.h"
#include "system.h"

using namespace std;

void Multibody_List::listloop(Two_Body_Mean_Orientational_Gyration_Tensor)
{
  int current_time;
  
  current_time = convert_time();
  
  for(int multibodyii=0;multibodyii<multibodies[current_time];multibodyii++)
	{
		analysis->listkernel(multibodies[current_time][multibodyii]);
	}
}