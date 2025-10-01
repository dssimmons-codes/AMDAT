
/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for class to hold a 3-D coordinate*/
/*Written by David S. Simmons*/

#include "coordinate.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>

using namespace std;



Coordinate::Coordinate(const Coordinate & copy)
{
  x=copy.x;
  y=copy.y;
  z=copy.z;
}



/*less-than operator; returns true if the first point is smaller in ANY dimension*/
bool Coordinate::operator < (const Coordinate & sample)const
{
  if(x < sample.x || y < sample.y || z < sample.z)
  {return 1;}
  else
  {return 0;}
}



/*greater-than operator; returns true if the first point is larger in ANY dimension*/
bool Coordinate::operator > (const Coordinate & sample)const
{
  if(x > sample.x || y > sample.y || z > sample.z)
  {return 1;}
  else
  {return 0;}
}


Coordinate Coordinate::integer()const
{
  Coordinate temp;
  temp.x = float(int(x));
  temp.y = float(int(y));
  temp.z = float(int(z));

  return temp;

}



Coordinate Coordinate::coord_floor()const
{
  Coordinate temp;

  temp.x = std::floor(x);
  temp.y = std::floor(y);
  temp.z = std::floor(z);

  return temp;
}


Coordinate Coordinate::coord_ceil()const
{
  Coordinate temp;

  temp.x = std::ceil(x);
  temp.y = std::ceil(y);
  temp.z = std::ceil(z);

  return temp;
}

Coordinate Coordinate::coord_round()const
{
  Coordinate temp;
  
  temp.x = std::round(x);
  temp.y = std::round(y);
  temp.z = std::round(z);
  
  return temp;
  
}


/*Methods to calculate shortest vector, considering box crossing*/
Coordinate Coordinate::vector_unwrapped(const Coordinate& boxsize)const
{
	float length; 
	float pbcx = x;
	float pbcy = y;
	float pbcz = z;

	if(x > boxsize.x/2)
		pbcx = x - boxsize.x;
	else if(x < -boxsize.x/2)
		pbcx = x + boxsize.x;

	if(y > boxsize.y/2)
		pbcy = y - boxsize.y;
	else if(y < -boxsize.y/2)
		pbcy = y + boxsize.y;

	if(z > boxsize.z/2)
		pbcz = z - boxsize.z;
	else if(z < -boxsize.z/2)
		pbcz = z + boxsize.z;

//	cout << " " << pbcx << " " << pbcy << " " << pbcz << "\n";
	Coordinate pbcr(pbcx, pbcy, pbcz);
	return pbcr;
}

/*Methods to calculate length of shortest vector, considering box crossing*/
float Coordinate::length_unwrapped(const Coordinate& boxsize)const
{
	float length;
	float minx, miny, minz;
	minx = min(std::abs(x),boxsize.x-std::abs(x));
	miny = min(std::abs(y),boxsize.y-std::abs(y));
	minz = min(std::abs(z),boxsize.z-std::abs(z));

	length = sqrtf(minx*minx+miny*miny+minz*minz);

	return length;
}


Coordinate Coordinate::closest_image(const Coordinate& other, const Coordinate& boxsize)const
{
  Coordinate imageflag(0,0,0);
  Coordinate diff=other-*this;
  if(boxsize.x-std::abs(diff.x)<std::abs(diff.x))
  {
    if(std::abs(x-(other.x+boxsize.x))<std::abs(x-(other.x-boxsize.x)))
    {
      imageflag.x=1;
    }
    else
    {
      imageflag.x=-1;
    }
  }
  if(boxsize.y-std::abs(diff.y)<std::abs(diff.y))
  {
    if(std::abs(y-(other.y+boxsize.y))<std::abs(y-(other.y-boxsize.y)))
    {
      imageflag.y=1;
    }
    else
    {
      imageflag.y=-1;
    }
  }
  if(boxsize.z-std::abs(diff.z)<std::abs(diff.z))
  {
    if(std::abs(z-(other.z+boxsize.z))<std::abs(z-(other.z-boxsize.z)))
    {
      imageflag.z=1;
    }
    else
    {
      imageflag.z=-1;
    }
  }
  return imageflag;
}

float Coordinate::min2D()const
{
	float minimum = x;
	if(y<minimum) minimum = y;

	return minimum;
}

float Coordinate::min()const
{
	float minimum = x;
	if(y<minimum) minimum = y;
	if(z<minimum) minimum = z;

	return minimum;
}

float Coordinate::max()const
{
	float maximum = x;
	if(y>maximum) maximum = y;
	if(z>maximum) maximum = z;

	return maximum;
}

float Coordinate::min(float a,float b)const
{
	if(a<b) return a;
	else return b;
}


Coordinate Coordinate::unit_vector()const
{
  Coordinate temp;
  float invl = 1/length();
  
  temp.x=x*invl;
  temp.y=y*invl;
  temp.z=z*invl;
  
  return temp;
}


float Coordinate::operator()(int index)const
{
  if(index == 0)
  {
    return x;
  }
  else if(index == 1)
  {
    return y;
  }
  else if(index == 2)
  {
    return z;
  }
  else
  {
    return 0;
    cout << "invalid index: coordinate index should be 0 for x, 1 for y, or 2 for z"<<endl;
    exit(1);
  }
}


void Coordinate::smallest(const Coordinate * coordlist, int listsize)
{
  x=coordlist[0].x;
  y=coordlist[0].y;
  z=coordlist[0].z;
  for(int coordii=1;coordii<listsize;coordii++)
  {
    if(coordlist[coordii].x<x) x=coordlist[coordii].x;
    if(coordlist[coordii].y<y) y=coordlist[coordii].y;
    if(coordlist[coordii].z<z) z=coordlist[coordii].z;
  }
}


bool Coordinate::within(const Coordinate& low, const Coordinate & high) const
{
  if(x<=high.x&&y<=high.y&&z<=high.z&&x>=low.x&&y>=low.y&&z>=low.z)
  {
    return true;
  }
  else
  {
    return false;
  }
  
}




bool Coordinate::operator!=(const Coordinate& comparator)
{
  return(x!=comparator.x||y!=comparator.y||z!=comparator.z);
}



bool Coordinate::operator==(const Coordinate& comparator)
{
  return(x==comparator.x&&y==comparator.y&&z==comparator.z);
}




float Coordinate::sum()const
{
  return x+y+z;
}
