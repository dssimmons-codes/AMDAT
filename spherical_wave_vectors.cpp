/*Spherical_Wave_Vectors class methods.  This class bins a grind of wave vectors into spherical shells of wavenumber.*/
/*Written by David S. Simmons*/

using namespace std;

#include "spherical_wave_vectors.h"
#include <math.h>

void Spherical_Wave_Vectors::bin(int xii, int yii, int zii)
{
	Coordinate temp;
	int wavenumberii;
	int ii;
	
	temp.set(xii,yii,zii);
	wavenumberii=floor((temp.length()+delta_wavenumber)/delta_wavenumber);  //round to nearest wavenumber
	if(n_wavevectors[wavenumberii] >= binsize[wavenumberii])
	{
		Coordinate * transfer;
		transfer = new Coordinate [n_wavevectors[wavenumberii]];
		for(ii=0;ii<n_wavevectors[wavenumberii];ii++)
		{
			transfer[ii]=wavevector[wavenumberii][ii];
		}
		delete [] wavevector[wavenumberii];
		binsize[wavenumberii]*=2;
		wavevector[wavenumberii] = new Coordinate [binsize[wavenumberii]];
		for(ii=0;ii<n_wavevectors[wavenumberii];ii++)
		{
			wavevector[wavenumberii][ii]=transfer[ii];
		}
		delete [] transfer;
	}
	wavevector[wavenumberii][n_wavevectors[wavenumberii]]=temp;
	n_wavevectors[wavenumberii]++;
}