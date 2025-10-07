// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
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

    /*Algebraic operators*/
    Coordinate& operator+=(const Coordinate& o) noexcept {x+=o.x; y+=o.y; z+=o.z; return *this;}
    Coordinate& operator-=(const Coordinate& o) noexcept {x-=o.x; y-=o.y; z-=o.z; return *this;}
    Coordinate& operator*=(const Coordinate& o) noexcept {x*=o.x; y*=o.y; z*=o.z; return *this;}
    Coordinate& operator/=(const Coordinate& o) noexcept {x/=o.x; y/=o.y; z/=o.z; return *this;}
    Coordinate& operator*=(float s) noexcept {x*=s; y*=s; z*=s; return *this;}
    Coordinate& operator/=(float s) noexcept {x/=s; y/=s; z/=s; return *this;}

    friend Coordinate operator+(Coordinate a, const Coordinate& b) noexcept {a += b; return a;}
    friend Coordinate operator-(Coordinate a, const Coordinate& b) noexcept {a -= b; return a;}
    friend Coordinate operator*(Coordinate a, const Coordinate& b) noexcept {a *= b; return a;}
    friend Coordinate operator/(Coordinate a, const Coordinate& b) noexcept {a /= b; return a;}
    friend Coordinate operator*(Coordinate a, float s) noexcept {a *= s; return a;}
    friend Coordinate operator*(float s, Coordinate a) noexcept {a *= s; return a;}
    friend Coordinate operator/(Coordinate a, float s) noexcept {a /= s; return a;}

    /*Vector operators*/
    float operator& (const Coordinate& multiplier) const {return x*multiplier.x+y*multiplier.y+z*multiplier.z;}			//dot product

    /*Boolean operators*/
    bool operator < (const Coordinate&)const;
    bool operator > (const Coordinate&)const;
    Coordinate integer()const;		//truncate decimals
    Coordinate coord_floor()const;			//floor command
    Coordinate coord_ceil()const;    //ceil command
    Coordinate coord_round()const;		//round;

    /*Methods to calculate in-box vector length*/
    float length() const noexcept {return std::sqrt(x*x + y*y + z*z);}  //returns vector length
    float length_xy() const noexcept {return std::sqrt(x*x+y*y);}       //returns vector length in xy plane
    float length_yz() const noexcept {return std::sqrt(y*y+z*z);}       //returns vector length in yz plane
    float length_xz() const noexcept {return std::sqrt(x*x+z*z);}       //returns vector length in xz plane
    float length_sq() const noexcept {return x*x+y*y+z*z;}			    //returns vector length squared
    
    float length_unwrapped(const Coordinate&)const;	//returns length of shortest vector, allowing box wrapping
    Coordinate vector_unwrapped(const Coordinate&)const;	//returns shortest vector, allowing box wrapping
    Coordinate closest_image(const Coordinate&, const Coordinate&)const;	//returns index (integer) coordinate providing image correction that must be applied to get the shortest distance
    float min()const;			//returns minimum of three coordinate values
    float max()const;			//returns maximum of three coordinate values
    
    Coordinate unit_vector() const;

    float show_x()const noexcept{return x;};
    float show_y()const noexcept{return y;};
    float show_z()const noexcept{return z;};


    float operator () (int index)const;
    
    void smallest(const Coordinate *, int);		//from a list of coordinates, sets this coordinate to  the smallest values of x, y, and z in the list
    bool within(const Coordinate& low, const Coordinate & high) const;	//returns true if this coordinate is within the rectangular prism with low and high as lower and upper corners; returns false otherwise
    
    bool operator!=(const Coordinate&);
    bool operator==(const Coordinate&);
    
    float sum()const;
};

}

#endif
