/*Molecular Dynamics Analysis Toolkit (MDAT)*/
/*Header for class Coordinate, which holds 3-D spacial coordinate.*/
/*Written by David S. Simmons*/

#ifndef COORDINATE
#define COORDINATE

#include <math.h>

namespace std {

class Coordinate
{
    float x, y, z;
    
    float min(float,float)const;

  public:
    Coordinate(){x=0;y=0;z=0;};					//default constructor
    Coordinate(const Coordinate &);				//copy constructor
    Coordinate(float a, float b, float c){x=a;y=b;z=c;};	//Fully defined constructor
    
    void set(float a, float b, float c){x=a;y=b;z=c;};		//set coordinate value

    Coordinate operator- (const Coordinate&) const;
    Coordinate operator+ (const Coordinate&) const;
    Coordinate operator* (const Coordinate&) const;		//component product
    Coordinate operator/ (const Coordinate&) const;
    float operator& (const Coordinate&) const;			//dot product
    Coordinate operator/ (float) const;
    Coordinate operator* (float) const;		
    void operator-= (const Coordinate&);
    void operator+= (const Coordinate&);
    void operator/= (const Coordinate&);
    void operator*= (const Coordinate&);
    void operator/= (float);
    void operator*= (float);
    bool operator < (const Coordinate&)const;
    bool operator > (const Coordinate&)const;
    Coordinate integer()const;		//truncate decimals
    Coordinate coord_floor()const;			//floor command
    float length() const;			//returns vector length
    float length_xy() const;
    float length_yz() const;
    float length_xz() const;
    float length_unwrapped(const Coordinate&)const;	//returns length of shortest vector, allowing box wrapping
    float min()const;			//returns minimum of three coordinate values
    
    Coordinate unit_vector() const;

    float show_x()const{return x;};
    float show_y()const{return y;};
    float show_z()const{return z;};
 
    float operator () (int index)const;
    
    void smallest(const Coordinate *, int);		//from a list of coordinates, sets this coordinate to  the smallest values of x, y, and z in the list
};

}

#endif