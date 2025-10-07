// Copyright (C) 2010-2025 David S. Simmons

/*_Wave_Vectors class methods.  This class bins a grind of wave vectors into spherical shells of wavenumber.*/
/*Written by David S. Simmons*/

/*Environment variables WV3D, WV2D, and WV1D must be defined in the makefile. These variables set the paths to the 3d, 2d, and 1d qvector files, respectively*/

#include "wave_vectors.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "progress.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

using namespace std;

Wave_Vectors::Wave_Vectors()
{
  system=0;
  wavegrid_spacing = 0;		//spacing of grid of wave-vectors; by default equal to 2pi/L, where L is the minimum system dimension
  maxrange = 0;			//distance wavegrid extends (symmetrically) from origin along each axis
  delta_wavenumber = 0;		//thickness of spherical wavenumber shells
  plane = "";
  n_atoms_looped = 0;
  n_wavenumbers=99;
  //n_wavevectors=new int[n_wavenumbers];
  //approx_wavenumber=new float [n_wavenumbers];
  //binsize=new int [n_wavenumbers];
  //wavevector = new Coordinate*[n_wavenumbers];
  wavevector.resize(n_wavenumbers);
  for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    //n_wavevectors[wavenumberii]=0;
    //binsize[wavenumberii]=0;
    //approx_wavenumber[wavenumberii]=0;
    //wavevector[wavenumberii] = new Coordinate[binsize[wavenumberii]];
    //wavevector[wavenumberii][n_wavevectors[wavenumberii]] = 0;
  }
}

Wave_Vectors::Wave_Vectors(System * sys, int shellcount)
{
  calculate(sys, shellcount);
}

Wave_Vectors::Wave_Vectors(System * sys)
{
  system = sys;
  int wavenumberii;
  n_wavenumbers = 99;
  wavegrid_spacing = 2.0*PI/(system->min_box_dimensions()).min();	//define spacing to be given by wavenumber corresponding to smallest system dimension
  maxrange = wavegrid_spacing * 100.0;		//wavevector grid 100 times this length in each direction from the origin
  delta_wavenumber = wavegrid_spacing/2.0;	//define the thickness of the wavenumber bins to be half a gridspacing

  //allocate memory
  //n_wavevectors = new int [n_wavenumbers];
  //wavevector = new Coordinate * [n_wavenumbers];
  wavevector.resize(n_wavenumbers);
  //binsize = new int [n_wavenumbers];
  //approx_wavenumber = new float [n_wavenumbers];
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    //n_wavevectors[wavenumberii]=0;  //zero out the number of wave vectors per shell
    //vectorcount_estimate = (4.0/3.0)*PI*pow(delta_wavenumber/wavegrid_spacing,3.0)*(3*wavenumberii*wavenumberii+3*wavenumberii+1);
    //binsize[wavenumberii]=int(vectorcount_estimate*1.5);	//estimate number of vectors per shell
    //wavevector[wavenumberii] = new Coordinate [binsize[wavenumberii]];
    //binsize[wavenumberii]=0;
    approx_wavenumber.push_back(wavegrid_spacing+delta_wavenumber*wavenumberii);
  }
  read_vectors();
}

Wave_Vectors::Wave_Vectors(System * sys, string pl, float max_length_scale)
{
    system = sys;
    int wavenumberii;
    n_wavenumbers = 99;
    plane = pl;
    float min_size;
    if (max_length_scale <= 0.0)
    {
      cout << "Maximum lengthscale for wavevector decomposition <= 0; defaulting to box size.\n";
      if (plane=="xy")
      {
	if((system->min_box_dimensions()).show_x()<(system->min_box_dimensions()).show_y())
	{
	  min_size=(system->min_box_dimensions()).show_x();
	}
	else
	{
	  min_size=(system->min_box_dimensions()).show_y();
	}
      }
      else if (plane == "xz")
      {
	if((system->min_box_dimensions()).show_x()<(system->min_box_dimensions()).show_z())
	{
	  min_size=(system->min_box_dimensions()).show_x();
	}
	else
	{
	  min_size=(system->min_box_dimensions()).show_z();
	}
      }
      else if (plane == "yz")
      {
	if((system->min_box_dimensions()).show_z()<(system->min_box_dimensions()).show_y())
	{
	  min_size=(system->min_box_dimensions()).show_z();
	}
	else
	{
	  min_size=(system->min_box_dimensions()).show_y();
	}
      }
      else if (plane == "xyz")
      {
	min_size = (system->min_box_dimensions()).min();
      }
      else if (plane == "x")
      {
	min_size = (system->min_box_dimensions()).show_x();
      }
      else if (plane == "y")
      {
	min_size = (system->min_box_dimensions()).show_y();
      }
      else if (plane == "z")
      {
	min_size = (system->min_box_dimensions()).show_z();
      }
      else
      {
	cout << "ERROR: plane command not recognized. valid inputs are: x, y, z, xy, xz, yz, xyz"<<endl; exit(1);
      }
    }
    else
    {
      if (max_length_scale > system->min_box_dimensions().min())
      {
	cout << endl<<"WARNING: Maximum length scale larger than smallest box dimension. Consider choosing a different maximum for wavevector decomposition.";
      }
      min_size = max_length_scale;
    }

  wavegrid_spacing = 2.0*PI/min_size;	//define spacing to be given by wavenumber corresponding to smallest system dimension
  maxrange = wavegrid_spacing * 100.0;		//wavevector grid 100 times this length in each direction from the origin
  delta_wavenumber = wavegrid_spacing/2.0;	//define the thickness of the wavenumber bins to be half a gridspacing

  //allocate memory
  //n_wavevectors = new int [n_wavenumbers];
  //wavevector = new Coordinate * [n_wavenumbers];
  wavevector.resize(n_wavenumbers);
  //binsize = new int [n_wavenumbers];
  //approx_wavenumber = new float [n_wavenumbers];
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    //n_wavevectors[wavenumberii]=0;  //zero out the number of wave vectors per shell
    //vectorcount_estimate = (4.0/3.0)*PI*pow(delta_wavenumber/wavegrid_spacing,3.0)*(3*wavenumberii*wavenumberii+3*wavenumberii+1);
    //binsize[wavenumberii]=int(vectorcount_estimate*1.5);	//estimate number of vectors per shell
    //wavevector[wavenumberii] = new Coordinate [binsize[wavenumberii]];
    //binsize[wavenumberii]=0;
    approx_wavenumber.push_back(wavegrid_spacing+delta_wavenumber*wavenumberii);
  }
  if (plane == "xyz")
  {
    read_vectors();
  }
  else if (plane == "xy" || plane == "xz" || plane == "yz")
  {
    read_vectors_2d(plane);
  }
  else
  {
    read_vectors_1d(plane);
  }
}

Wave_Vectors::~Wave_Vectors()
{
//   cout<<"Destructing Wave_Vectors\n";cout.flush();
//  for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
//  {
//    delete [] wavevector[wavenumberii];
//  }
//  delete [] wavevector;
//  delete [] binsize;
//  delete [] n_wavevectors;
//  delete [] approx_wavenumber;
}

Wave_Vectors::Wave_Vectors(const Wave_Vectors & copy)
{
  system = copy.system;
  wavegrid_spacing = copy.wavegrid_spacing;		//spacing of grid of wave-vectors; by default equal to 2pi/L, where L is the minimum system dimension
  maxrange = copy.maxrange;			//distance wavegrid extends (symmetrically) from origin along each axis
  delta_wavenumber = copy.delta_wavenumber;		//thickness of spherical wavenumber shells
  plane = copy.plane;
  n_atoms_looped = copy.n_atoms_looped;
  n_wavenumbers = copy.n_wavenumbers;

//  n_wavevectors = new int [n_wavenumbers];
//  approx_wavenumber = new float [n_wavenumbers];
//  binsize = new int [n_wavenumbers];
//  wavevector = new Coordinate*[n_wavenumbers];
  wavevector.resize(n_wavenumbers);
  approx_wavenumber=copy.approx_wavenumber;
  wavevector=copy.wavevector;
//  for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
//  {
//    n_wavevectors[wavenumberii]=copy.n_wavevectors[wavenumberii];
//    binsize[wavenumberii]=copy.binsize[wavenumberii];
//    approx_wavenumber[wavenumberii]=copy.approx_wavenumber[wavenumberii];
//    wavevector[wavenumberii]= new Coordinate[binsize[wavenumberii]];
//    wavevector[wavenumberii][n_wavevectors[wavenumberii]]=copy.wavevector[wavenumberii][n_wavevectors[wavenumberii]];
//  }
}

Wave_Vectors Wave_Vectors::operator = (const Wave_Vectors & copy)
{
  if(this!=&copy)
  {
    system = copy.system;
    wavegrid_spacing = copy.wavegrid_spacing;		//spacing of grid of wave-vectors; by default equal to 2pi/L, where L is the minimum system dimension
    maxrange = copy.maxrange;			//distance wavegrid extends (symmetrically) from origin along each axis
    delta_wavenumber = copy.delta_wavenumber;		//thickness of spherical wavenumber shells
    plane = copy.plane;
    n_atoms_looped = copy.n_atoms_looped;
    n_wavenumbers = copy.n_wavenumbers;
    
    approx_wavenumber=copy.approx_wavenumber;
    wavevector=copy.wavevector;

//    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
//    {
//      delete [] wavevector[wavenumberii];
//    }
//    delete [] wavevector;
//    delete [] binsize;
//    delete [] n_wavevectors;
//    delete [] approx_wavenumber;

//    n_wavevectors = new int [n_wavenumbers];
//    approx_wavenumber = new float [n_wavenumbers];
//    binsize = new int [n_wavenumbers];
//    wavevector = new Coordinate*[n_wavenumbers];

//    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
//    {
//      n_wavevectors[wavenumberii]=copy.n_wavevectors[wavenumberii];
//      binsize[wavenumberii]=copy.binsize[wavenumberii];
//      approx_wavenumber[wavenumberii]=copy.approx_wavenumber[wavenumberii];
//      wavevector[wavenumberii]= new Coordinate[binsize[wavenumberii]];
//      wavevector[wavenumberii][n_wavevectors[wavenumberii]]=copy.wavevector[wavenumberii][n_wavevectors[wavenumberii]];
//    }
  }
  return *this;
}


/*Method to read vectors in from file.*/
void Wave_Vectors::read_vectors_2d(string plane)
{
  string filename;
  char filestem[4];
  int wavenumberii;
  ifstream vectorfile;
  float tempx, tempy, tempz;
  Coordinate tempvec;

  cout << "\n";
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    print_progress(wavenumberii,n_wavenumbers-1);
    sprintf(filestem, ".%03d",wavenumberii+2);	//generate stem-number of vector file
    filename=WV2D;
    filename += filestem;

    vectorfile.open(filename.c_str());		//open file
    if(!vectorfile.is_open()){cout << "\nError opening vector file" << filename <<"\n.";exit(1);}

    /*first determine number of vectors*/
    while(!vectorfile.eof())			//loop until end of file
    {
      vectorfile >> tempx;
      if(vectorfile.eof()){tempz=1; break;}  //check whether end of file has been reached.  This protects against blank lines at the end of the file tricking the loop check for end of file after the last vector.
      vectorfile >> tempy >> tempz;	//temporarily copy coordinates of vector from file
      //binsize[wavenumberii]++;
    }
    //binsize[wavenumberii]++;
    //wavevector[wavenumberii] = new Coordinate [binsize[wavenumberii]];

    vectorfile.close();
    vectorfile.open(filename.c_str());		//open file
    if(!vectorfile.is_open()){cout << "\nError opening vector file"<< filename <<"\n.";exit(1);}

    while(!vectorfile.eof())			//loop until end of file
    {
      //if(n_wavevectors[wavenumberii]>=binsize[wavenumberii])
      //{cout<<"Error: vector memory allocation error."; exit(1);}
      vectorfile >> tempx;
      if(vectorfile.eof()){break;}  //check whether end of file has been reached.  This protects against blank lines at the end of the file tricking the loop check for end of file after the last vector.
      vectorfile >> tempy >> tempz;	//temporarily copy coordinates of vector from file

        if (plane=="xy")
        {
//          wavevector[wavenumberii][n_wavevectors[wavenumberii]].set(tempx,tempy,tempz);	//store vector in array
          tempvec.set(tempx,tempy,tempz);
	  tempvec*= wavegrid_spacing;
          wavevector[wavenumberii].push_back(tempvec);	//store vector in array
        }
        if (plane=="xz")
        {
//          wavevector[wavenumberii][n_wavevectors[wavenumberii]].set(tempx,tempz,tempy);	//store vector in array
	  tempvec.set(tempx,tempz,tempy);
	  tempvec*= wavegrid_spacing;
	  wavevector[wavenumberii].push_back(tempvec);	//store vector in array
        }
        if (plane=="yz")
        {
//          wavevector[wavenumberii][n_wavevectors[wavenumberii]].set(tempz,tempx,tempy);	//store vector in array
	  tempvec.set(tempz,tempx,tempy);
	  tempvec*= wavegrid_spacing;
	  wavevector[wavenumberii].push_back(tempvec);	//store vector in array
        }

//        wavevector[wavenumberii][n_wavevectors[wavenumberii]] *= wavegrid_spacing;		//scale vector appropriately
//      if(wavenumberii==1){cout << n_wavevectors[wavenumberii] << "\t" << wavevector[wavenumberii][n_wavevectors[wavenumberii]].show_x() << "\t" << wavevector[wavenumberii][n_wavevectors[wavenumberii]].show_y() << "\t" << wavevector[wavenumberii][n_wavevectors[wavenumberii]].show_z() << "\n";}
 //       n_wavevectors[wavenumberii]++;		//increment number of wavevectors stored.

    }
      vectorfile.close();
  }

  n_atoms_looped = 0;

}




/*Method to read vectors in from file.*/
/*could easily be improved - there is really no reason to store this data in a file.*/
void Wave_Vectors::read_vectors_1d(string plane)
{
  string filename;
  char filestem[4];
  int wavenumberii;
  ifstream vectorfile;
  float tempx, tempy, tempz;
  Coordinate tempvec;
  
  
  cout << "\n";
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    print_progress(wavenumberii,n_wavenumbers-1);
    sprintf(filestem, ".%03d",wavenumberii+2);	//generate stem-number of vector file
    filename=WV1D;
    filename += filestem;

    vectorfile.open(filename.c_str());		//open file
    if(!vectorfile.is_open()){cout << "\nError opening vector file" << filename <<"\n.";exit(1);}

    /*first determine number of vectors*/
    while(!vectorfile.eof())			//loop until end of file
    {
      vectorfile >> tempx;
      if(vectorfile.eof()){tempz=1; break;}  //check whether end of file has been reached.  This protects against blank lines at the end of the file tricking the loop check for end of file after the last vector.
      vectorfile >> tempy >> tempz;	//temporarily copy coordinates of vector from file
      //binsize[wavenumberii]++;
    }
    //binsize[wavenumberii]++;
    //wavevector[wavenumberii] = new Coordinate [binsize[wavenumberii]];

    vectorfile.close();
    vectorfile.open(filename.c_str());		//open file
    if(!vectorfile.is_open()){cout << "\nError opening vector file"<< filename <<"\n.";exit(1);}

    while(!vectorfile.eof())			//loop until end of file
    {
      //if(n_wavevectors[wavenumberii]>=binsize[wavenumberii])
      //{cout<<"Error: vector memory allocation error."; exit(1);}
      vectorfile >> tempx;
      if(vectorfile.eof()){break;}  //check whether end of file has been reached.  This protects against blank lines at the end of the file tricking the loop check for end of file after the last vector.
      vectorfile >> tempy >> tempz;	//temporarily copy coordinates of vector from file

        if (plane=="x")
        {
          //wavevector[wavenumberii][n_wavevectors[wavenumberii]].set(tempx,0,0);	//store vector in array
	  tempvec.set(tempx,0,0);
	  tempvec*= wavegrid_spacing;
	  wavevector[wavenumberii].push_back(tempvec);	//store vector in array
        }
        if (plane=="y")
        {
//          wavevector[wavenumberii][n_wavevectors[wavenumberii]].set(0,tempx,0);	//store vector in array
	  tempvec.set(0,tempx,0);
	  tempvec*= wavegrid_spacing;
	  wavevector[wavenumberii].push_back(tempvec);	//store vector in array
        }
        if (plane=="z")
        {
	  tempvec.set(0,0,tempx);
	  tempvec*= wavegrid_spacing;
	  wavevector[wavenumberii].push_back(tempvec);	//store vector in array
          //wavevector[wavenumberii][n_wavevectors[wavenumberii]].set(0,0,tempx);	//store vector in array
        }

//        wavevector[wavenumberii][n_wavevectors[wavenumberii]] *= wavegrid_spacing;		//scale vector appropriately
//      if(wavenumberii==1){cout << n_wavevectors[wavenumberii] << "\t" << wavevector[wavenumberii][n_wavevectors[wavenumberii]].show_x() << "\t" << wavevector[wavenumberii][n_wavevectors[wavenumberii]].show_y() << "\t" << wavevector[wavenumberii][n_wavevectors[wavenumberii]].show_z() << "\n";}
//        n_wavevectors[wavenumberii]++;		//increment number of wavevectors stored.

    }
      vectorfile.close();
  }

  n_atoms_looped = 0;

}




/*Method to read vectors in from file.*/
void Wave_Vectors::read_vectors()
{
  string filename;
  char filestem[4];
  int wavenumberii;
  ifstream vectorfile;
  float tempx, tempy, tempz;
  Coordinate tempvec;

  cout << "\n";
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    print_progress(wavenumberii,n_wavenumbers-1);
    sprintf(filestem, ".%03d",wavenumberii+2);	//generate stem-number of vector file
    filename=WV3D;
    filename += filestem;

    vectorfile.open(filename.c_str());		//open file
    if(!vectorfile.is_open()){cout << "\nError opening vector file" << filename <<"\n.";exit(1);}

    /*first determine number of vectors*/
    while(!vectorfile.eof())			//loop until end of file
    {
      vectorfile >> tempx >> tempy >> tempz;
      //binsize[wavenumberii]++;
    }
    //wavevector[wavenumberii] = new Coordinate [binsize[wavenumberii]];

    vectorfile.close();
    vectorfile.open(filename.c_str());		//open file
    if(!vectorfile.is_open()){cout << "\nError opening vector file"<< filename <<"\n.";exit(1);}

    while(!vectorfile.eof())			//loop until end of file
    {
      //if(n_wavevectors[wavenumberii]>=binsize[wavenumberii])
      //{cout<<"Error: vector memory allocation error."; exit(1);}
      vectorfile >> tempx;
      if(vectorfile.eof()){break;}  //check whether end of file has been reached.  This protects against blank lines at the end of the file tricking the loop check for end of file after the last vector.
      vectorfile >> tempy >> tempz;	//temporarily copy coordinates of vector from file
      
      tempvec.set(tempx,tempy,tempz);
      tempvec*= wavegrid_spacing;
      wavevector[wavenumberii].push_back(tempvec);	//store vector in array

//      wavevector[wavenumberii][n_wavevectors[wavenumberii]].set(tempx,tempy,tempz);	//store vector in array
//      wavevector[wavenumberii][n_wavevectors[wavenumberii]] *= wavegrid_spacing;		//scale vector appropriately
//      if(wavenumberii==1){cout << n_wavevectors[wavenumberii] << "\t" << wavevector[wavenumberii][n_wavevectors[wavenumberii]].show_x() << "\t" << wavevector[wavenumberii][n_wavevectors[wavenumberii]].show_y() << "\t" << wavevector[wavenumberii][n_wavevectors[wavenumberii]].show_z() << "\n";}
//      n_wavevectors[wavenumberii]++;		//increment number of wavevectors stored.
    }
        vectorfile.close();
  }

  n_atoms_looped = 0;

}

void Wave_Vectors::calculate(System * sys, int shellcount)
{
  int xii, yii, zii;			//indices over cartesian components of wavevectors
  int wavenumberii;
  int vectorcount_estimate;

  system = sys;

  n_wavenumbers = shellcount;

  wavegrid_spacing = 2.0*PI/(system->min_box_dimensions()).min();	//define spacing to be given by wavenumber corresponding to smallest system dimension
  maxrange = wavegrid_spacing * 100.0;		//for now, just use a wavevector grid 100 times this length in each direction from the origin
  delta_wavenumber = wavegrid_spacing/2.0;	//for now, just define the thickness of the wavenumber bins to be half a gridspacing

  //allocate memory
//  n_wavevectors = new int [n_wavenumbers];
//  wavevector = new Coordinate * [n_wavenumbers];
//  binsize = new int [n_wavenumbers];
//  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
//  {
//    n_wavevectors[wavenumberii]=0;  //zero out the number of wave vectors per shell
//    vectorcount_estimate = (4.0/3.0)*PI*pow(delta_wavenumber/wavegrid_spacing,3.0)*(3*wavenumberii*wavenumberii+3*waveumberii+1);
//    wavevector[wavenumberii] = new Coordinate [int(vectorcount_estimate*2)];
//    binsize[wavenumberii]=vectorcount_estimate*2;
//  }


  //loop over independent q vectors in grid, assigning each to the appropriate shell as determined for each geometry in bin().
  for(xii=-maxrange;xii<=maxrange;xii+=wavegrid_spacing)
  {
    for(yii=-maxrange;yii<=maxrange;yii+=wavegrid_spacing)
    {
      for(zii=0;zii<=maxrange;zii+=wavegrid_spacing)
      {
        bin(xii,yii,zii);
      }
    }
  }
}


float Wave_Vectors::show_mean_wavenumber(int index)const
{
  int vectorii;
  float mean_wavenumber = 0;
  for(vectorii=0;vectorii<wavevector[index].size();vectorii++)
  {
    mean_wavenumber+=wavevector[index][vectorii].length();
  }
  mean_wavenumber/=wavevector[index].size();
  return mean_wavenumber;
}



float Wave_Vectors::show_stdev_wavenumber(int index)const
{
  int vectorii;
  float variance=0;
  float stdev;
  float mean_wavenumber = show_mean_wavenumber(index);
  for(vectorii=0;vectorii<wavevector[index].size();vectorii++)
  {
    variance+=pow((wavevector[index][vectorii].length()-mean_wavenumber),2.0);
  }
  variance/=wavevector[index].size();
  stdev=pow(variance,0.5);
  return stdev;
}

Coordinate Wave_Vectors::show_mean_wavevector(int index)const
{
  int vectorii;
  Coordinate mean_wavevector(0,0,0);
  for(vectorii=0;vectorii<wavevector[index].size();vectorii++)
  {
    mean_wavevector+=wavevector[index][vectorii];
  }
  mean_wavevector/=wavevector[index].size();
  return mean_wavevector;
}




#ifdef NEVER
void Wave_Vectors::calculate(Coordinate boxsize, float deltak, float kmax, int maxvectors)
{
  float kx, ky, kz, k;
  float deltakx, deltaky, deltakz;
  int n_kx, n_ky, n_kz
  Coordinate boxsize;
  float kymax, kzmax;
  system=sys;
  
  deltakx=2*PI/boxsize.show_x();
  deltaky=2*PI/boxsize.show_y();
  deltakz=2*PI/boxsize.show_z();
  
  for(kx=deltakx;kx<kmax;kx+=deltakx)
  {
    kymax = pow(kmax*kmax-kx*kx,0.5);
    for(ky=deltaky;ky<kymax;ky+=deltaky)
    {
      kzmax = pow(kmax*kmax-kz*kz,0.5);
      for(kz=deltakyz;kz<kzmax;kz+=deltakz)
      {
	k = pow(kx*kx+ky*ky+kz*kz,0.5);
	
      }
    }
  }
  
  

  
}

#endif
