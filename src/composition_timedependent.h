// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Class outputting information about system composition as a function of time*/
/*Is exactly the same as parent class Composition, but changes standard output to include timedependence*/
/*Written by David Simmons*/


#ifndef COMPOSITION_TIMEDEPENDENT
#define COMPOSITION

#include "analysis.h"
#include "system.h"
#include <string>
#include "analysis_onetime.h"

namespace std{


class Composition_TimeDependent:public Composition
{
    public:

    Composition_TimeDependent():Composition(){;};

    Composition_TimeDependent(System * sys, int n_xbins, int n_ybins, int n_zbins,float lx, float ly, float lz, int timescheme=-1):Composition(sys,n_xbins,n_ybins,n_zbins,lx,ly,lz,timescheme){;};

    Composition_TimeDependent(const Composition & copy):Composition(copy){;};

    ~Composition_TimeDependent(){;};

    void write(string);
    void write(ofstream&);
};
}

#endif // COMPOSITION
