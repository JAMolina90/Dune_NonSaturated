#include <iostream>
#include <time.h>
#include "Globals.h"
#include "geometry.h"
#include "parameter.h"
#include "matrix.h"
#include "Getparam.h"
#include "Getgeometry.h"
#include "CalculateSFD.h"
#include "CBSM.h"

int main()
{
    /*VARIABLES DEFINITION*/
    //Geometry variables definition
    geometry geo;

    //Unknowns definition
    matrix<nprec> unknowns;

    //Parameters
    parameter par;

    /*GEOMETRY PROCESSING*/
    //Parameters Collection
    Getparam(par);

    //Geometry Collection
    Getgeometry(geo,par,unknowns);

    //Matrix of shape functions derivatives
    matrix<nprec> SFD;
    CalculateSFD(SFD,geo);

    /*CBS Method*/
    clock_t t;
    t=clock();

    CBSM(geo,par,unknowns,SFD);

    t=clock()-t;
    std::cout << "Elapsed Time: " << ((float)t)/CLOCKS_PER_SEC << " seconds" << std::endl;
    return 0;
}
