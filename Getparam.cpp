#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <string.h>
#include <stdlib.h>
#include "Globals.h"
#include "parameter.h"

void Getparam(parameter& par)
{
    std::ifstream readfile("parameters.dat");
    assert(readfile.is_open());
    readfile.ignore(1000,'\n');
    readfile >> par.g >> par.dvisco0 >> par.rho0 >> par.c0 >> par.dvisco1 >> par.rho1 >> par.c1 >> par.surft >> par.Uref >> par.Lref;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.t0 >> par.ntime >> par.niteru >> par.nitert >> par.nplt;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.pcgtol >> par.itpcgmax;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.dyntimestep >> par.delt >> par.maxCou >> par.maxCouH;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.ppre >> par.upre >> par.vpre;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.theta1 >> par.theta2 >> par.theta3 >> par.theta4;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.restart >> par.nplotrestart;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.masscheck >> par.massstep;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.resCheck >> par.tolerance;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.GInput;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.spx >> par.spy;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.slipnx >> par.slipny;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.nsp;
    if(par.nsp>0)
    {
        par.SPc.SetDimensions(par.nsp,2);
        par.SPn.SetDimensions(par.nsp,2);
        for(unsigned int j=0; j<par.nsp; j++)
        {
            readfile >> par.SPc(j,0) >> par.SPc(j,1);
            readfile >> par.SPn(j,0) >> par.SPn(j,1);
        };
    };
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.method;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.hydroP;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.indexfctr >> par.deps >> par.pardts >> par.tolfctt;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.zoneworeinit >> par.x1zwor >> par.y1zwor >> par.x2zwor >> par.y2zwor;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.indexfctu;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.SCM >> par.z00 >> par.angfri >> par.veltresh >> par.bdens >> par.aval >> par.graind >> par.tfactor;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.nz0;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.nsmoothuast;
    for(int i=0; i<2; i++){readfile.ignore(1000,'\n');};
    readfile >> par.SLalpha >> par.SLgamma >> par.SLdtfactor;
    readfile.close();
    par.DELTP=par.delt*par.Uref/par.Lref;
    par.rhoref=par.rho0;
    if (par.dvisco0>1e-9)
    {
        par.dviscoref=par.dvisco0;
    } else if(par.dvisco1>1e-9)
    {
        par.dviscoref=par.dvisco1;
    } else
    {
        par.dviscoref=(nprec)1.0;
    };
    par.Re=par.dviscoref/(par.rhoref*par.Uref*par.Lref);
    par.Fr2=(par.g*par.Lref)/(par.Uref*par.Uref);
    par.We=(par.surft)/(par.rhoref*par.Uref*par.Uref*par.Lref);
    par.c0=par.c0/(par.Uref*par.Uref);
    par.c1=par.c1/(par.Uref*par.Uref);

    if(par.nplt>0){par.jw=par.ntime/par.nplt;} else {par.jw=0;};
    if(par.jw==1){std::cout << "WARNING: code writes a plot file every timestep!" << std::endl;};
    if(par.nplt>100){std::cout << "WARNING: code writes " << par.nplt << " plot files" << std::endl;};
    if(par.pcgtol<1e-20){par.pcgtol=1e-20; std::cout << "WARNING: PCGSolver tolerance has been modified, new value is " << par.pcgtol << std::endl;};
    if(std::abs(par.Uref)<1e-9){std::cout << "Uref cannot be 0" << std::endl;exit(1);};
    if(std::abs(par.Lref)<1e-9){std::cout << "Lref cannot be 0" << std::endl;exit(1);};

    if(par.restart)
    {
        std::string tmp1;
        std::ifstream readfile("nsvtk2D_final.vtk");
        if (!readfile.is_open()){std::cout << "Final vtk file not found" << std::endl; exit(1);};
        readfile.ignore(1000,'\n');
        readfile >> tmp1 >> tmp1 >> par.t0;
        readfile.close();
    };

    if (par.SCM!=1&&par.SCM!=2)
    {
        std::cout << "Stress Calculation Method not well defined" << std::endl;
        exit(1);
    };
    if (par.nz0<0.0)
    {
        std::cout << "Height where velocity is taken can not be negative" << std::endl;
        exit(1);
    };

    if (par.z00<(nprec)1e-15)
    {
        par.z00=(nprec)1e-15;
    };

    par.tfinal=par.t0+(nprec)par.ntime*par.delt*par.Uref/par.Lref;
    par.interplot=((nprec)par.ntime)/((nprec)par.nplt)*par.delt*par.Uref/par.Lref;
}
