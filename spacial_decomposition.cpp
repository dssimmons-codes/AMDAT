/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for Spacial_Decomposition class*/
/*Written by David S. Simmons*/

/*Note that this code currently does not work properly with particle lists due to ambiguities in how lists that change in time should be handled.*/

#include <math.h>
#include "spacial_decomposition.h"
#include <iostream>
#include "progress.h"
#include <stdlib.h>
#include "version.h"

using namespace std;


Spacial_Decomposition::Spacial_Decomposition(System* sys, float min_cell_size, int capacity, int over)
{
  initialize(sys, min_cell_size, capacity, over);
}



/*-------------------------------------------------------------------------------------*/


void Spacial_Decomposition::initialize(System* sys, float min_box_size, int capacity, int over)
{
  int indexii;
  
  system = sys;
  overhang = over;
  
  if(min_box_size > system->size().min())
  {
    cout << "Error: minimum box size greater than smallest system dimension.";
    exit(1);
  }
  
  //calculate number of cells in each direction, rounding down so as to make cell equal to or larger than minimum specified size
  n_real_cells[0] = int(system->size().show_x()/min_box_size);
  n_real_cells[1] = int(system->size().show_y()/min_box_size);
  n_real_cells[2] = int(system->size().show_z()/min_box_size);

 //calculate cell size corresponding to number of cells above
  cell_size.set((system->size().show_x())/float(n_real_cells[0]),(system->size().show_y())/float(n_real_cells[1]),(system->size().show_z())/float(n_real_cells[2]));
  
  
  n_timesteps = system->show_n_timesteps();
  
  xmin = system->boundaries()[0].show_x();
  xmax = system->boundaries()[1].show_x();
  ymin = system->boundaries()[0].show_y();
  ymax = system->boundaries()[1].show_y();
  zmin = system->boundaries()[0].show_z();
  zmax = system->boundaries()[1].show_z();
  
  //add overhang cells to number of cells
  n_cells[0] = n_real_cells[0] + 2*overhang;
  n_cells[1] = n_real_cells[1] + 2*overhang;
  n_cells[2] = n_real_cells[2] + 2*overhang;
  
  cout << "\nAllocating memory to spacial decomposition array.";cout.flush();
  
  total_cells = n_timesteps*n_cells[0]*n_cells[1]*n_cells[2];
  spacial_cells = new Coordinate * [total_cells];
  atomID = new int * [total_cells];
  cell_capacity = new int [total_cells];
  atom_count = new int [total_cells];
  
 
  for(indexii=0;indexii<total_cells;indexii++)
  {
   cell_capacity[indexii] = capacity;
   atom_count[indexii] = 0;
   spacial_cells[indexii] = new Coordinate [cell_capacity[indexii]];
   atomID[indexii]=new int [cell_capacity[indexii]];
  }
  
  countmax=0;
}



/*-------------------------------------------------------------------------------------*/



Spacial_Decomposition::~Spacial_Decomposition()
{
  clear_memory();
}


/*-------------------------------------------------------------------------------------*/



void Spacial_Decomposition::clear_memory()
{
  int index;
  for(index=0;index<total_cells;index++)
  {
    delete [] spacial_cells[index];
    delete [] atomID [index];
  }
  delete [] spacial_cells;
  delete [] atomID;
  delete [] atom_count;
  delete [] cell_capacity;
}



/*-------------------------------------------------------------------------------------*/

int Spacial_Decomposition::convert_index(int timeii, int xii, int yii, int zii)const
{
  return timeii*n_cells[0]*n_cells[1]*n_cells[2] + (xii+overhang)*n_cells[1]*n_cells[2] + (yii+overhang)*n_cells[2]+(zii+overhang);
}




/*-------------------------------------------------------------------------------------*/


void Spacial_Decomposition::atomkernel(Trajectory * traj)
{
  
  int timeii, xii, yii, zii, index;
  Coordinate coordinate;
 

  for(timeii=0;timeii<n_timesteps;timeii++)	//set pointer to given atom
  {

    coordinate=(traj->show_coordinate(timeii));	//copy coordinate from present atom
    xii=int((coordinate.show_x()-xmin)/(cell_size.show_x()));	//calculate x cell index
    yii=int((coordinate.show_y()-ymin)/(cell_size.show_y()));	//calculate y cell index
    zii=int((coordinate.show_z()-zmin)/(cell_size.show_z()));	//calculate z cell index
   
    
    if(xii<0||xii>=n_real_cells[0]||yii<0||yii>=n_real_cells[1]||zii<0||zii>=n_real_cells[2])
    {
      cout << coordinate.show_x() <<"\t"<< coordinate.show_y() <<"\t"<< coordinate.show_z()<<"\n";
      cout << "Error: atom outside stated system boundary!\n";
      exit(1);
    }
    
    index = convert_index(timeii,xii,yii,zii);
    /*Check if there is room in cell; if not, increase memory allocation*/
    /*This process is very time-inefficient, should just increase base capacity if it happens much.*/

    if(atom_count[index] == cell_capacity[index])
      {change_capacity(index, 2*cell_capacity[index]);}
    
    (spacial_cells[index][atom_count[index]]) = coordinate;	//put coordinate in appropriate cell in array	
     atomID[index][atom_count[index]] = traj->show_trajectory_ID();
     
    (atom_count[index])++;		//increment number of atoms in that cell
    if((atom_count[index])>countmax){countmax = atom_count[index];}
    
  }
}


/*-------------------------------------------------------------------------------------*/


void Spacial_Decomposition::atomlist_kernel(int species_index, int molecule_index, int atom_type, int atom_index)
{
	int timeii, xii, yii, zii, index;
	Coordinate coordinate;
 
	cout<<"\n\natomlistkernel\n\n";
	 
	Atom_Trajectory * atompointer;		//pointer to present atom

	atompointer = ((system->show_molecule(species_index,molecule_index))->show_atom_trajectory(atom_type,atom_index));

	for(timeii=0;timeii<n_timesteps;timeii++)	//set pointer to given atom
	{

		coordinate=(atompointer->show_coordinate(timeii));	//copy coordinate from present atom
		xii=int((coordinate.show_x()-xmin)/(cell_size.show_x()));	//calculate x cell index
		yii=int((coordinate.show_y()-ymin)/(cell_size.show_y()));	//calculate y cell index
		zii=int((coordinate.show_z()-zmin)/(cell_size.show_z()));	//calculate z cell index
   		cout<<"\n"<<xii<<"\t"<<yii<<"\t"<<zii<<"\t";
    
		if(xii<0||xii>=n_real_cells[0]||yii<0||yii>=n_real_cells[1]||zii<0||zii>=n_real_cells[2])
		{
			cout << coordinate.show_x() <<"\t"<< coordinate.show_y() <<"\t"<< coordinate. show_z()<<"\n";
			cout << "Error: atom outside stated system boundary!\n";
			exit(1);
		}
    
		index = convert_index(timeii,xii,yii,zii);
		/*Check if there is room in cell; if not, increase memory allocation*/
		/*This process is very time-inefficient, should just increase base capacity if it happens much.*/

		if(atom_count[index] == cell_capacity[index])
		{change_capacity(index, 2*cell_capacity[index]);}
    
		(spacial_cells[index][atom_count[index]]) = coordinate;	//put coordinate in appropriate cell in array	
		atomID[index][atom_count[index]] = atompointer->show_atomID();
     
		(atom_count[index])++;		//increment number of atoms in that cell
		if((atom_count[index])>countmax){countmax = atom_count[index];}
    
	}
}



/*-------------------------------------------------------------------------------------*/


void Spacial_Decomposition::change_capacity(int index, int capacity)
{
      int atomii;
      Coordinate * temp;
      int * tempID;
      
      temp = new Coordinate [cell_capacity[index]];
      tempID = new int [cell_capacity[index]];
      for(atomii=0;atomii<cell_capacity[index];atomii++)
      {
        temp[atomii]=spacial_cells[index][atomii];	//copy cell to temp
        tempID[atomii] = atomID[index][atomii];
      }
      
      delete [] spacial_cells[index];	//deallocate cell memory
      delete [] atomID[index];
      (cell_capacity[index]) = capacity;	//double cell capactiy
      spacial_cells[index] = new Coordinate [cell_capacity[index]];	//reallocate twice the amount of cell memory
      atomID[index] = new int [cell_capacity[index]];
      
      for(atomii=0;atomii<cell_capacity[index];atomii++)
      {
        spacial_cells[index][atomii]=temp[atomii];	//copy back to cell
        atomID[index][atomii]=tempID[atomii];
      }

}



void Spacial_Decomposition::fill_overhang()
{
	
  int index;
  int copyindex;
  float xsize = xmax-xmin;
  float ysize = ymax-ymin;
  float zsize = zmax-zmin;
  int tii, xii, yii, zii, atomii;
  int xshift=0, yshift=0, zshift=0;
  int overhangcount=0;
  
  Coordinate xshifter(1,0,0);
  Coordinate yshifter(0,1,0);
  Coordinate zshifter(0,0,1);
  
  xshifter = xshifter * xsize;
  yshifter = yshifter * ysize;
  zshifter = zshifter * zsize;
  
  int n_overhang_cells=n_timesteps*(n_cells[0]*n_cells[1]*n_cells[2]-(n_cells[0]-2*overhang)*(n_cells[1]-2*overhang)*(n_cells[2]-2*overhang));
  
  for(tii=0;tii<n_timesteps;tii++)
  {
    for(xii=-overhang;xii<n_cells[0]-overhang;xii++)
    {
      for(yii=-overhang;yii<n_cells[1]-overhang;yii++)
      {
        for(zii=-overhang;zii<n_cells[2]-overhang;zii++)
        {
          if(zii==0 && (xii>=0 && xii<n_real_cells[0]) && (yii>=0 && yii<n_real_cells[1]))
          {zii+=(n_real_cells[2]-1);} //skip over existing cells
          //determine which cell to duplicate and which way to shift coordinates
          else
          {
            overhangcount++;

            if(xii>=n_real_cells[0])
            {
              //xshift = -1;
	      xshift = -int(xii/n_real_cells[0]);
            }
            else if (xii < 0)
            {
              //xshift = 1;
	      xshift = -(int(float(xii)/float(n_real_cells[0]))-1);
            }
            else xshift=0;
            
            if(yii>=n_real_cells[1])
            {
              //yshift = -1;
	      yshift = -int(yii/n_real_cells[1]);
            }
            else if (yii < 0)
            {
              //yshift = 1;
	      yshift = -(int(float(yii)/float(n_real_cells[1]))-1);
            }
            else yshift = 0;
            
            if(zii>=n_real_cells[2])
            {
              //zshift = -1;
	      zshift = -int(zii/n_real_cells[2]);
            }
            else if (zii < 0)
            {
              //zshift = 1;
	      zshift = -(int(float(zii)/float(n_real_cells[2]))-1);
            }
            else zshift=0;
            
            index = convert_index(tii,xii,yii,zii); 	//calculate index of present cell
            
            copyindex = convert_index(tii,xii+xshift*n_real_cells[0],yii+yshift*n_real_cells[1],zii+zshift*n_real_cells[2]);	//calculate index of cell to be copied
            
            atom_count[index]=atom_count[copyindex];	//copy atom count
            
            if(cell_capacity[index]!=cell_capacity[copyindex])		//copy cell_capacity if necessary
              {change_capacity(index,cell_capacity[copyindex]);}
            print_progress(overhangcount,n_overhang_cells);
            for(atomii=0;atomii<atom_count[index];atomii++)	//loop over atoms within cell
            {
              spacial_cells[index][atomii] = spacial_cells[copyindex][atomii];		//copy value of atom trajectory (not pointer)
              atomID[index][atomii] = atomID[copyindex][atomii];
              spacial_cells[index][atomii] -= (xshifter*xshift + yshifter*yshift + zshifter*zshift);
            }
          }
        } 
      }
    }
  }
}



/*loops over atoms in bin specified by an offset from bin of a specified atom; passes down timegap and coordinates of both atoms to cellkernel of analysis method*/
void Spacial_Decomposition::loop_cell(Analysis * analysis, int timegap, int nextii, Coordinate coordinate1, int atom1ID, int deltax, int deltay, int deltaz)const
{
  int xii, yii, zii;
  int atomii;
  
  //determine cell to be looped over
  xii=int((coordinate1.show_x()-xmin)/(cell_size.show_x()))+deltax;	//calculate x cell index
  yii=int((coordinate1.show_y()-ymin)/(cell_size.show_y()))+deltay;;	//calculate y cell index
  zii=int((coordinate1.show_z()-zmin)/(cell_size.show_z()))+deltaz;;	//calculate z cell index

  //determine index of cell
  int index = convert_index(nextii,xii,yii,zii);

  for(atomii=0; atomii < atom_count[index];atomii++)	//loop over atoms in cell
  {
    if(atom1ID != atomID[index][atomii])
    {
      analysis->cellkernel(timegap, coordinate1, spacial_cells[index][atomii]);
    }
  }
}


/*loops over atoms in bin specified by an offset from bin of a specified atom; passes down timegap and IDs of both atoms to cellkernelID of analysis method*/
void Spacial_Decomposition::loop_cell_ID(Analysis * analysis, int thisii, int nextii, int timegapii, Coordinate coordinate1, int atom1ID, int deltax, int deltay, int deltaz)const
{
	int xii, yii, zii;
	int atomii;
  
  //determine cell to be looped over
	xii=int((coordinate1.show_x()-xmin)/(cell_size.show_x()))+deltax;	//calculate x cell index
	yii=int((coordinate1.show_y()-ymin)/(cell_size.show_y()))+deltay;;	//calculate y cell index
	zii=int((coordinate1.show_z()-zmin)/(cell_size.show_z()))+deltaz;;	//calculate z cell index

  //determine index of cell
	int index = convert_index(nextii,xii,yii,zii);
	//cout << "\n" << atom_count[index];
	for(atomii=0; atomii < atom_count[index];atomii++)	//loop over atoms in cell
	{
		//cout<<"\n"<<atom1ID<<"\t"<<atomID[index][atomii];
		if(atom1ID != atomID[index][atomii])
		{
			analysis->cellkernelID(thisii, nextii, timegapii, atom1ID, atomID[index][atomii], coordinate1, spacial_cells[index][atomii]);
		}
	}
}