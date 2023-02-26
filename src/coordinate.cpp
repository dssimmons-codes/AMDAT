
/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Methods for class to hold a 3-D coordinate*/
/*Written by David S. Simmons*/

#include <math.h>
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


Coordinate Coordinate::coord_round()const
{
  Coordinate temp;
  
  temp.x = round(x);
  temp.y = round(y);
  temp.z = round(z);
  
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


float Coordinate::length_sq()const
{
	float lengthsq;
	lengthsq = x*x+y*y+z*z;
	return lengthsq;
}


/*Methods to calculate length of shortest vector, considering box crossing*/
Coordinate Coordinate::vector_unwrapped(const Coordinate& boxsize)const
{
	float length; 
	float pbcx, pbcy, pbcz;

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

	Coordinate pbcr(pbcx, pbcy, pbcz);
	return pbcr;
}

/*Methods to calculate length of shortest vector, considering box crossing*/
float Coordinate::length_unwrapped(const Coordinate& boxsize)const
{
	float length;
	float minx, miny, minz;
	minx = min(abs(x),boxsize.x-abs(x));
	miny = min(abs(y),boxsize.y-abs(y));
	minz = min(abs(z),boxsize.z-abs(z));

	length = pow((minx*minx+miny*miny+minz*minz),.5);

	return length;
}


Coordinate Coordinate::closest_image(const Coordinate& other, const Coordinate& boxsize)const
{
  Coordinate imageflag(0,0,0);
  Coordinate diff=other-*this;
  if(boxsize.x-abs(diff.x)<abs(diff.x))
  {
    if(abs(x-(other.x+boxsize.x))<abs(x-(other.x-boxsize.x)))
    {
      imageflag.x=1;
    }
    else
    {
      imageflag.x=-1;
    }
  }
  if(boxsize.y-abs(diff.y)<abs(diff.y))
  {
    if(abs(y-(other.y+boxsize.y))<abs(y-(other.y-boxsize.y)))
    {
      imageflag.y=1;
    }
    else
    {
      imageflag.y=-1;
    }
  }
  if(boxsize.z-abs(diff.z)<abs(diff.z))
  {
    if(abs(z-(other.z+boxsize.z))<abs(z-(other.z-boxsize.z)))
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


bool Coordinate::within(const Coordinate& low, const Coordinate & high)
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
