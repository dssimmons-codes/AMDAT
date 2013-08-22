/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Function to display percentage-progress on a task.*/
/*Written by David S. Simmons*/


#include <iostream>
#include <stdio.h>
//#include "version.h"

using namespace std;
#ifdef NO

#include "progress.h"

void print_progress(int done,int total)
{
	float percent;
	//if(!fileoutput)
	//{
	  percent = int(100*float(done)/float(total));
	  cout.width(3);
	  cout << "\b\b\b\b" << percent << "%";
	  cout.flush();
	//}
}

#endif