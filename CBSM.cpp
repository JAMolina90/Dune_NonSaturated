#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "Globals.h"
#include "myvector.h"
#include "matrix.h"
#include "geometry.h"
#include "parameter.h"
#include "AuxVar.h"
#include "Neighbour.h"
#include "UpwindVar.h"
#include "SStepMatrix.h"
#include "SetBCI.h"
#include "UpdateRhoVisco.h"
#include "H.h"
#include "FirstStep.h"
#include "PreSecondStep.h"
#include "SecondStep.h"
#include "Get3old.h"
#include "Getend.h"
#include "FirstStepUpwind.h"
#include "ThirdStepUpwind.h"
#include "FCTU.h"
#include "VTKoutput.h"

void CBSM(const geometry& geo0, parameter& par, matrix<nprec>& unknowns0, const matrix<nprec>& SFD0)
{
    geometry geo1;
    matrix<nprec> unknowns1(unknowns0.GetRows(),unknowns0.GetColumns());
    matrix<nprec> SFD1;

    //Initialize vectors and matrix for the System of Linear Equations
    AuxVar av0,av1;
    AuxVarGeo avG;
    AuxVarIni(av0,geo0,SFD0,par);
    AuxVarGeoIni(avG,geo0);

    //Variables for Upwind Schemes
    UpwindVar upvel0;
    UpwindVar upvel1;
    matrix<nprec> RHSUW1;

    if (par.indexfctu!=0)
    {
        Upwindgeo(geo0,upvel0);
        //RHSUW.SetDimensions(geo.nd.GetNnodes(),3);//La columna 0 se queda libre por si se quiere meter ahí algo relativo a las presiones en el futuro
    };

    //Elementos vecinos
    Neig NB0;
    Neighbour(geo0,NB0,avG);
    avG.EdgeCase.SetDimensions(avG.Edges.GetRows());
    avG.EdgeCase.Initialize(1);
    avG.PrevEdgecutted.SetDimensions(avG.Edges.GetRows());
    avG.PrevEdgecutted.Initialize(false);
    avG.PrevEdgeValue.SetDimensions(avG.Edges.GetRows());
    avG.PrevEdgeValue.Initialize((nprec)0.0);
    avG.NodeNumberEdge.SetDimensions(avG.Edges.GetRows());

    //Indexes for writings
    int kw=0;
    int kf=0;
    if(par.restart)
    {
        if(par.nplotrestart==0)
        {
            kf=par.nplt;
        }
        else
        {
            kf=par.nplotrestart;
        };
        std::cout << "VTK files will have index from " << kf << std::endl;
    };

    //Set boundary an Initial conditions
    SetBCI(geo0, par, unknowns0, av0);
    unknowns1=unknowns0;

    //Initial Diffusion of interface
    //Reinitialization(geo,par,av,SFD,unknowns,upreinit);
    //UpdateRhoVisco(geo0,unknowns0,par,av0);

    /*Hvalues*/
    HnodesCalculationIni(avG,geo0);

    //New Mesh
    NewMesh(geo0,geo1,par,avG,SFD0,SFD1,upvel0,upvel1,RHSUW1,unknowns0,unknowns1);
    AuxVarIni(av1,geo1,SFD1,par);
    SetBCI(geo1,par,unknowns1,av1);

    /*Initial VTK output*/
    if (!par.restart)
    {
        VTKoutput(0,0,geo1,par,unknowns1,false,av1,avG);
    };

    //Loop over timestep
    nprec time;
    if (!par.dyntimestep)
    {
    time=par.t0;
    for (unsigned int itime=0; itime<par.ntime; itime++)
    {
        time+=par.DELTP;
        std::cout << "itime = " << itime+1 << std::endl;

        //Pressure interpolation for non-activated nodes
        PressureInterpolation(geo0,geo1,unknowns0,unknowns1,avG);

        UastCalculation(avG,geo1,unknowns1,par,SFD1,NB0);
        Sediment(av1,avG,geo1,unknowns1,par,SFD1);
        
	if (time>=60.0)
        {
        HnodesCalculation(avG,geo1,unknowns1,par,SFD1,0,0);

        //New Mesh
        if ((itime+1)%10==0)
	{
	   NewMesh(geo0,geo1,par,avG,SFD0,SFD1,upvel0,upvel1,RHSUW1,unknowns0,unknowns1);
           AuxVarIni(av1,geo1,SFD1,par);
	   std::cout << "New Mesh Calculated" << std::endl;
	};
        };
        SetBCI(geo1,par,unknowns1,av1);

        av1.RHS.Initialize((nprec)0.0);
        //Velocity Predictor
        //InterfNormalsCalculation(geo,av,SFD,unknowns,NB);
        //STCalculation(geo,par,SFD,av,unknowns);

        //Activar la linea de abajo si se quiere compatibilizar tensiones en el fondo en FirstStep y las calculadas para sedimento
        //UastCalculation(avG,geo1,unknowns1,par,SFD1,NB0);

        FirstStep(geo1,par,av1,SFD1,unknowns1,avG,0);

        //Second Step Matrix
        SStepMatrix(geo1, par, SFD1, av1, unknowns1);

        //Surf Tension Term
        //InterfNormalsCalculation(geo,av,SFD,unknowns,NB);
        //STCalculation(geo,par,SFD,av,unknowns);

        //First Step (Navier-Stokes)
        if (par.indexfctu!=2)
        {
            ConvTerm(geo1,par,av1,SFD1,unknowns1,avG);
        };

        if (par.indexfctu!=0)
        {
            RHSUW1.Initialize((nprec)0.0);

            //First Step Upwind
            FirstStepUpwind(geo1,par,av1,unknowns1,RHSUW1,upvel1);

            //FCT
            if(par.indexfctu==1)
            {
                matrix<nprec> FS(geo1.nd.GetNnodes(),2);
                FCTU(par,geo1,SFD1,av1.DMAT,vectorize(unknowns1,1,"column"),vectorize(unknowns1,2,"column"),\
                     vectorize(RHSUW1,1,"column")+vectorize(unknowns1,1,"column"),\
                     vectorize(RHSUW1,2,"column")+vectorize(unknowns1,2,"column"),\
                     vectorize(av1.ConvTerm,0,"column")+vectorize(unknowns1,1,"column"),\
                     vectorize(av1.ConvTerm,1,"column")+vectorize(unknowns1,2,"column"),FS,vectorize(unknowns1,3,"column"));
                for(unsigned int i=0; i<geo1.nd.GetNnodes(); i++)
                {
                    av1.ConvTerm(i,0)=FS(i,0)-unknowns1.Read(i,1);
                    av1.ConvTerm(i,1)=FS(i,1)-unknowns1.Read(i,2);
                };
            } else
            {
                for(unsigned int i=0; i<geo1.nd.GetNnodes(); i++)
                {
                    av1.ConvTerm(i,0)=RHSUW1(i,1);
                    av1.ConvTerm(i,1)=RHSUW1(i,2);
                };
            };
        };

        FirstStep(geo1,par,av1,SFD1,unknowns1,avG,1);

        //Second Step
        PreSecondStep(geo1,par,SFD1,unknowns1,av1);
        SecondStep(geo1,par,av1);

        //if (par.indexfctu!=2)
        {
            //Perform pressure correction for velocities
            Get3old(geo1,par,av1,SFD1,unknowns1);
        };


        //Final Computation
        Getend(geo1,par,av1,SFD1,unknowns1);

        //Final Result
        for(unsigned int i=0; i<geo1.nd.GetNnodes(); i++)
        {
            unknowns1(i,0)+=av1.RHS.Read(i,0);
            unknowns1(i,1)+=av1.RHS.Read(i,1);
            unknowns1(i,2)+=av1.RHS.Read(i,2);
        };

        avG.PrevEdgecutted.Initialize(false);
        avG.PrevEdgeValue.Initialize((nprec)0.0);
        for (unsigned int i=0; i<avG.Edges.GetRows(); i++)
        {
            if (avG.EdgeCase.Read(i)==2)
            {
                avG.PrevEdgecutted(i)=true;
                avG.PrevEdgeValue(i)=unknowns1(avG.NodeNumberEdge.Read(i),0);
            };
        };

        //Writes files for movie: vtk-i
        //output for ParaView
        if (par.nplt>1)
        {
            kw++;
            if (kw==par.jw)
            {
                kf++;
                kw=0;
                VTKoutput(kf,time,geo1,par,unknowns1,false,av1,avG);
            };
        };

        for (unsigned int i=0; i<geo0.nd.GetNnodes(); i++)
        {
            if (avG.ActivatedNodes.Read(i)==0)
            {
                unknowns0(i,1)=(nprec)0.0;
                unknowns0(i,2)=(nprec)0.0;
                unknowns0(i,3)=(nprec)0.0;
            };
        };
    };
    }
    else
    {
    time=par.t0;
    unsigned int itime=0;
    nprec timepar=0.0;
    while (time<par.tfinal)
    {
        itime++;
        std::cout << "itime = " << itime << std::endl;

        TimeStep(geo1,SFD1,unknowns1,par,timepar,time);

        //Pressure interpolation for non-activated nodes
        PressureInterpolation(geo0,geo1,unknowns0,unknowns1,avG);

        //Desactivar la linea de abajo si se quiere compatibilizar tensiones en el fondo en FirstStep y las calculadas para sedimento
        UastCalculation(avG,geo1,unknowns1,par,SFD1,NB0);

        Sediment(av1,avG,geo1,unknowns1,par,SFD1);

        HnodesCalculation(avG,geo1,unknowns1,par,SFD1,timepar,time);

        time+=par.DELTP;
        timepar+=par.DELTP;
        std::cout << "Time-step = " << par.DELTP << std::endl;
        std::cout << "Total time = " << time << std::endl;

        //New Mesh
        NewMesh(geo0,geo1,par,avG,SFD0,SFD1,upvel0,upvel1,RHSUW1,unknowns0,unknowns1);
        AuxVarIni(av1,geo1,SFD1,par);
        SetBCI(geo1,par,unknowns1,av1);

        av1.RHS.Initialize((nprec)0.0);
        //Velocity Predictor
        //InterfNormalsCalculation(geo,av,SFD,unknowns,NB);
        //STCalculation(geo,par,SFD,av,unknowns);

        //Activar la linea de abajo si se quiere compatibilizar tensiones en el fondo en FirstStep y las calculadas para sedimento
        //UastCalculation(avG,geo1,unknowns1,par,SFD1,NB0);

        FirstStep(geo1,par,av1,SFD1,unknowns1,avG,0);

        //Second Step Matrix
        SStepMatrix(geo1, par, SFD1, av1, unknowns1);

        //Surf Tension Term
        //InterfNormalsCalculation(geo,av,SFD,unknowns,NB);
        //STCalculation(geo,par,SFD,av,unknowns);

        //First Step (Navier-Stokes)
        if (par.indexfctu!=2)
        {
            ConvTerm(geo1,par,av1,SFD1,unknowns1,avG);
        };

        if (par.indexfctu!=0)
        {
            RHSUW1.Initialize((nprec)0.0);

            //First Step Upwind
            FirstStepUpwind(geo1,par,av1,unknowns1,RHSUW1,upvel1);

            //FCT
            if(par.indexfctu==1)
            {
                matrix<nprec> FS(geo1.nd.GetNnodes(),2);
                FCTU(par,geo1,SFD1,av1.DMAT,vectorize(unknowns1,1,"column"),vectorize(unknowns1,2,"column"),\
                     vectorize(RHSUW1,1,"column")+vectorize(unknowns1,1,"column"),\
                     vectorize(RHSUW1,2,"column")+vectorize(unknowns1,2,"column"),\
                     vectorize(av1.ConvTerm,0,"column")+vectorize(unknowns1,1,"column"),\
                     vectorize(av1.ConvTerm,1,"column")+vectorize(unknowns1,2,"column"),FS,vectorize(unknowns1,3,"column"));
                for(unsigned int i=0; i<geo1.nd.GetNnodes(); i++)
                {
                    av1.ConvTerm(i,0)=FS(i,0)-unknowns1.Read(i,1);
                    av1.ConvTerm(i,1)=FS(i,1)-unknowns1.Read(i,2);
                };
            } else
            {
                for(unsigned int i=0; i<geo1.nd.GetNnodes(); i++)
                {
                    av1.ConvTerm(i,0)=RHSUW1(i,1);
                    av1.ConvTerm(i,1)=RHSUW1(i,2);
                };
            };
        };

        FirstStep(geo1,par,av1,SFD1,unknowns1,avG,1);

        //Second Step
        PreSecondStep(geo1,par,SFD1,unknowns1,av1);
        SecondStep(geo1,par,av1);

        //if (par.indexfctu!=2)
        {
            //Perform pressure correction for velocities
            Get3old(geo1,par,av1,SFD1,unknowns1);
        };


        //Final Computation
        Getend(geo1,par,av1,SFD1,unknowns1);

        //Final Result
        for(unsigned int i=0; i<geo1.nd.GetNnodes(); i++)
        {
            unknowns1(i,0)+=av1.RHS.Read(i,0);
            unknowns1(i,1)+=av1.RHS.Read(i,1);
            unknowns1(i,2)+=av1.RHS.Read(i,2);
        };

        avG.PrevEdgecutted.Initialize(false);
        avG.PrevEdgeValue.Initialize((nprec)0.0);
        for (unsigned int i=0; i<avG.Edges.GetRows(); i++)
        {
            if (avG.EdgeCase.Read(i)==2)
            {
                avG.PrevEdgecutted(i)=true;
                avG.PrevEdgeValue(i)=unknowns1(avG.NodeNumberEdge.Read(i),0);
            };
        };

        //Writes files for movie: vtk-i
        //output for ParaView
        if (par.nplt>1)
        {
            if ((timepar-1e-15<par.interplot)&&(timepar+1e-15>par.interplot))
            {
                kf++;
                timepar=(nprec)0.0;
                VTKoutput(kf,time,geo1,par,unknowns1,false,av1,avG);
            };
        };

        for (unsigned int i=0; i<geo0.nd.GetNnodes(); i++)
        {
            if (avG.ActivatedNodes.Read(i)==0)
            {
                unknowns0(i,1)=(nprec)0.0;
                unknowns0(i,2)=(nprec)0.0;
                unknowns0(i,3)=(nprec)0.0;
            };
        };
    };
    }

    /*Final VTK output*/
    VTKoutput(0,time,geo1,par,unknowns1,true,av1,avG);
    //VTKoutput(1,par.DELTP*par.ntime,geoRef,par,unknownsRef,false,avRef);
}
