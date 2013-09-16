
/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for class to hold a 3-D coordinate*/
/*Written by David S. Simmons*/

#include <math.h>
#include "coordinate.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

using namespace std;



Coordinate::Coordinate(const Coordinate & copy)
{
  x=copy.x;
  y=copy.y;
  z=copy.z;
}


Coordinate Coordinate::operator- (const Coordinate & decrement)const
{
	Coordinate temp;
	temp.x = x - decrement.x;
	temp.y = y - decrement.y;
	temp.z = z - decrement.z;

	return temp;
}



Coordinate Coordinate::operator+ (const Coordinate & increment)const
{
  Coordinate temp;
  temp.x = x + increment.x;
  temp.y = y + increment.y;
  temp.z = z + increment.z;

  return temp;
}



Coordinate Coordinate::operator* (const Coordinate & factor)const
{
  Coordinate temp;
  temp.x = x * factor.x;
  temp.y = y * factor.y;
  temp.z = z * factor.z;

  return temp;
}



Coordinate Coordinate::operator/ (const Coordinate & divisor)const
{
  Coordinate temp;
  temp.x = x / divisor.x;
  temp.y = y / divisor.y;
  temp.z = z / divisor.z;

  return temp;
}



/*Dot product method*/
float Coordinate::operator& (const Coordinate & multiplier)const
{
  return x*multiplier.x+y*multiplier.y+z*multiplier.z;
}



Coordinate Coordinate::operator/ (float divisor)const
{
	Coordinate temp;
	temp.x = x / divisor;
	temp.y = y / divisor;
	temp.z = z / divisor;

	return temp;
}



Coordinate Coordinate::operator* (float factor)const
{
	Coordinate temp;
	temp.x = x * factor;
	temp.y = y * factor;
	temp.z = z * factor;

	return temp;
}



void Coordinate::operator-= (const Coordinate & decrement)
{
  x -= decrement.x;
  y -= decrement.y;
  z -= decrement.z;
}



void Coordinate::operator+= (const Coordinate & increment)
{
  x += increment.x;
  y += increment.y;
  z += increment.z;
}



void Coordinate::operator*= (const Coordinate & multiplier)
{
  x*= multiplier.x;
  y*= multiplier.y;
  z*= multiplier.z;
}



void Coordinate::operator/= (const Coordinate & divisor)
{
  x/= divisor.x;
  y/= divisor.y;
  z/= divisor.z;
}



void Coordinate::operator*= (float multiplier)
{
  x*= multiplier;
  y*= multiplier;
  z*= multiplier;
}



void Coordinate::operator/= (float divisor)
{
  x/= divisor;
  y/= divisor;
  z/= divisor;
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

  temp.x = floor(x);
  temp.y = floor(y);
  temp.z = floor(z);

  return temp;
}

/*Methods to calculate in-box vector length*/
float Coordinate::length()const
{
	float length;
	//length = pow((pow(x,2)+pow(y,2)+pow(z,2)),.5);
	length = pow((x*x+y*y+z*z),.5);
	return length;
}

float Coordinate::length_xy()const
{
	float length;
	//length = pow((pow(x,2)+pow(y,2)+pow(z,2)),.5);
	length = pow((x*x+y*y),.5);
	return length;
}

float Coordinate::length_xz()const
{
	float length;
	//length = pow((pow(x,2)+pow(y,2)+pow(z,2)),.5);
	length = pow((x*x+z*z),.5);
	return length;
}

float Coordinate::length_yz()const
{
	float length;
	//length = pow((pow(x,2)+pow(y,2)+pow(z,2)),.5);
	length = pow((y*y+z*z),.5);
	return length;
}

/*Methods to calculate length of shortest vector, considering box crossing*/
float Coordinate::length_unwrapped(const Coordinate& boxsize)const
{
	float length;
	float minx, miny, minz;
	minx = min(x,boxsize.x-x);
	miny = min(y,boxsize.y-y);
	minz = min(z,boxsize.z-z);

	length = pow((minx*minx+miny*miny+minz*minz),.5);

	return length;
}


float Coordinate::min()const
{
	float minimum = x;
	if(y<minimum) minimum = y;
	if(z<minimum) minimum = z;

	return minimum;
}

float Coordinate::min(float a,float b)const
{
	if(a<b) return a;
	else return b;
}


Coordinate Coordinate::unit_vector()const
{
  Coordinate temp;
  float l = length();
  
  temp.x=x/l;
  temp.y=y/l;
  temp.z=z/l;
  
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
