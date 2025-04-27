#include <cmath>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include "Globals.h"
#include "myvector.h"
#include "matrix.h"
#include "geometry.h"
#include "parameter.h"
#include "AuxVar.h"
#include "UpwindVar.h"

void FirstStepUpwind(const geometry& geo, const parameter& par, AuxVar& av, const matrix<nprec>& unknowns, matrix<nprec>& RHSUW, UpwindVar& upvel)
{
    Vector<nprec> ux(geo.nd.GetNnodes());
    ux.Initialize((nprec)0.0);
    Vector<nprec> uy(geo.nd.GetNnodes());
    uy.Initialize((nprec)0.0);
    Vector<nprec> vx(geo.nd.GetNnodes());
    vx.Initialize((nprec)0.0);
    Vector<nprec> vy(geo.nd.GetNnodes());
    vy.Initialize((nprec)0.0);
    unsigned int i1,i2,i3;
    nprec auxu,auxv,fu,fv,u12,v12,un;

    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);

        auxu=(unknowns.Read(i1,1)+unknowns.Read(i2,1))/((nprec)2.0);
        auxv=(unknowns.Read(i1,2)+unknowns.Read(i2,2))/((nprec)2.0);
        ux(i1)+=auxu*upvel.GEOEL.Read(i,0)*upvel.GEOEL.Read(i,2);
        uy(i1)+=auxu*upvel.GEOEL.Read(i,1)*upvel.GEOEL.Read(i,2);
        ux(i2)-=auxu*upvel.GEOEL.Read(i,0)*upvel.GEOEL.Read(i,2);
        uy(i2)-=auxu*upvel.GEOEL.Read(i,1)*upvel.GEOEL.Read(i,2);
        vx(i1)+=auxv*upvel.GEOEL.Read(i,0)*upvel.GEOEL.Read(i,2);
        vy(i1)+=auxv*upvel.GEOEL.Read(i,1)*upvel.GEOEL.Read(i,2);
        vx(i2)-=auxv*upvel.GEOEL.Read(i,0)*upvel.GEOEL.Read(i,2);
        vy(i2)-=auxv*upvel.GEOEL.Read(i,1)*upvel.GEOEL.Read(i,2);

        auxu=(unknowns.Read(i2,1)+unknowns.Read(i3,1))/((nprec)2.0);
        auxv=(unknowns.Read(i2,2)+unknowns.Read(i3,2))/((nprec)2.0);
        ux(i2)+=auxu*upvel.GEOEL.Read(i,3)*upvel.GEOEL.Read(i,5);
        uy(i2)+=auxu*upvel.GEOEL.Read(i,4)*upvel.GEOEL.Read(i,5);
        ux(i3)-=auxu*upvel.GEOEL.Read(i,3)*upvel.GEOEL.Read(i,5);
        uy(i3)-=auxu*upvel.GEOEL.Read(i,4)*upvel.GEOEL.Read(i,5);
        vx(i2)+=auxv*upvel.GEOEL.Read(i,3)*upvel.GEOEL.Read(i,5);
        vy(i2)+=auxv*upvel.GEOEL.Read(i,4)*upvel.GEOEL.Read(i,5);
        vx(i3)-=auxv*upvel.GEOEL.Read(i,3)*upvel.GEOEL.Read(i,5);
        vy(i3)-=auxv*upvel.GEOEL.Read(i,4)*upvel.GEOEL.Read(i,5);

        auxu=(unknowns.Read(i3,1)+unknowns.Read(i1,1))/((nprec)2.0);
        auxv=(unknowns.Read(i3,2)+unknowns.Read(i1,2))/((nprec)2.0);
        ux(i3)+=auxu*upvel.GEOEL.Read(i,6)*upvel.GEOEL.Read(i,8);
        uy(i3)+=auxu*upvel.GEOEL.Read(i,7)*upvel.GEOEL.Read(i,8);
        ux(i1)-=auxu*upvel.GEOEL.Read(i,6)*upvel.GEOEL.Read(i,8);
        uy(i1)-=auxu*upvel.GEOEL.Read(i,7)*upvel.GEOEL.Read(i,8);
        vx(i3)+=auxv*upvel.GEOEL.Read(i,6)*upvel.GEOEL.Read(i,8);
        vy(i3)+=auxv*upvel.GEOEL.Read(i,7)*upvel.GEOEL.Read(i,8);
        vx(i1)-=auxv*upvel.GEOEL.Read(i,6)*upvel.GEOEL.Read(i,8);
        vy(i1)-=auxv*upvel.GEOEL.Read(i,7)*upvel.GEOEL.Read(i,8);
    };

    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        ux(i)/=upvel.VOL.Read(i);
        uy(i)/=upvel.VOL.Read(i);
        vx(i)/=upvel.VOL.Read(i);
        vy(i)/=upvel.VOL.Read(i);
    };

    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);
        if (av.dudn0(i))
        {
            if(std::abs(av.normalel.Read(i,0))<1e-9)
            {
                uy(i1)=(nprec)0.0;
                vy(i1)=(nprec)0.0;
                uy(i2)=(nprec)0.0;
                vy(i2)=(nprec)0.0;
                uy(i3)=(nprec)0.0;
                vy(i3)=(nprec)0.0;
            }else if (std::abs(av.normalel.Read(i,1))<1e-9)
            {
                ux(i1)=(nprec)0.0;
                vx(i1)=(nprec)0.0;
                ux(i2)=(nprec)0.0;
                vx(i2)=(nprec)0.0;
                ux(i3)=(nprec)0.0;
                vx(i3)=(nprec)0.0;
            }else
            {
                uy(i1)=-ux(i1)*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                vy(i1)=-vx(i1)*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                uy(i2)=-ux(i2)*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                vy(i2)=-vx(i2)*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                uy(i3)=-ux(i3)*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                vy(i3)=-vx(i3)*av.normalel.Read(i,0)/av.normalel.Read(i,1);
            };
        };
    };

    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        i1=geo.elem.ReadNode(i,0);i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);

        //Convective term (as flux)
        u12=(unknowns.Read(i1,1)+unknowns.Read(i2,1))/((nprec)2.0);
        v12=(unknowns.Read(i1,2)+unknowns.Read(i2,2))/((nprec)2.0);
        un=(u12*upvel.GEOEL.Read(i,0)+v12*upvel.GEOEL.Read(i,1))*upvel.GEOEL.Read(i,2);
        if (un>(nprec)0.0)
        {
            fu=un*unknowns.Read(i1,1);
            fv=un*unknowns.Read(i1,2);
        }else
        {
            fu=un*unknowns.Read(i2,1);
            fv=un*unknowns.Read(i2,2);
        };
        RHSUW(i1,1)-=fu;
        RHSUW(i2,1)+=fu;
        RHSUW(i1,2)-=fv;
        RHSUW(i2,2)+=fv;

        u12=(unknowns.Read(i2,1)+unknowns.Read(i3,1))/((nprec)2.0);
        v12=(unknowns.Read(i2,2)+unknowns.Read(i3,2))/((nprec)2.0);
        un=(u12*upvel.GEOEL.Read(i,3)+v12*upvel.GEOEL.Read(i,4))*upvel.GEOEL.Read(i,5);
        if (un>(nprec)0.0)
        {
            fu=un*unknowns.Read(i2,1);
            fv=un*unknowns.Read(i2,2);
        }else
        {
            fu=un*unknowns.Read(i3,1);
            fv=un*unknowns.Read(i3,2);
        };
        RHSUW(i2,1)-=fu;
        RHSUW(i3,1)+=fu;
        RHSUW(i2,2)-=fv;
        RHSUW(i3,2)+=fv;

        u12=(unknowns.Read(i3,1)+unknowns.Read(i1,1))/((nprec)2.0);
        v12=(unknowns.Read(i3,2)+unknowns.Read(i1,2))/((nprec)2.0);
        un=(u12*upvel.GEOEL.Read(i,6)+v12*upvel.GEOEL.Read(i,7))*upvel.GEOEL.Read(i,8);
        if (un>(nprec)0.0)
        {
            fu=un*unknowns.Read(i3,1);
            fv=un*unknowns.Read(i3,2);
        }else
        {
            fu=un*unknowns.Read(i1,1);
            fv=un*unknowns.Read(i1,2);
        };
        RHSUW(i3,1)-=fu;
        RHSUW(i1,1)+=fu;
        RHSUW(i3,2)-=fv;
        RHSUW(i1,2)+=fv;
    };

    nprec anx,any,aleng;
    for(unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        i1=geo.bd.ReadNode(i,0); i2=geo.bd.ReadNode(i,1);
        aleng=geo.bd.ReadLen(i)/((nprec)2.0);
        anx=geo.bd.ReadNormal(i,0)*aleng;
        any=geo.bd.ReadNormal(i,1)*aleng;

        if(unknowns.Read(i1,1)*anx+unknowns.Read(i1,2)*any>(nprec)0.0){RHSUW(i1,1)-=unknowns.Read(i1,1)*(unknowns.Read(i1,1)*anx+unknowns.Read(i1,2)*any);RHSUW(i1,2)-=unknowns.Read(i1,2)*(unknowns.Read(i1,1)*anx+unknowns.Read(i1,2)*any);};
        if(unknowns.Read(i2,1)*anx+unknowns.Read(i2,2)*any>(nprec)0.0){RHSUW(i2,1)-=unknowns.Read(i2,1)*(unknowns.Read(i2,1)*anx+unknowns.Read(i2,2)*any);RHSUW(i2,2)-=unknowns.Read(i2,2)*(unknowns.Read(i2,1)*anx+unknowns.Read(i2,2)*any);};
    };

    //Final Calculation
    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        RHSUW(i,1)*=par.DELTP/upvel.VOL.Read(i);
        RHSUW(i,2)*=par.DELTP/upvel.VOL.Read(i);
    };

    /*
    //BOUNDARY CONDITIONS
    //u prescribed boundary condition correction
    if (par.ksi)
    {
        unsigned int in1,in2;
        for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
        {
            if(geo.bd.ReadBC(i)==13)
            {
                in1=geo.bd.ReadNode(i,0); in2=geo.bd.ReadNode(i,1);
                RHSUW(in1,1)=(nprec)0.0;
                RHSUW(in1,2)=(nprec)0.0;
                RHSUW(in2,1)=(nprec)0.0;
                RHSUW(in2,2)=(nprec)0.0;
            };
        };
    };

    //Wall and slip condition correction
    if (geo.wnodes.GetRows()>0)
    {
        unsigned int ip;
        for (unsigned int i=0; i<geo.wnodes.GetRows(); i++)
        {
            ip=geo.wnodes.Read(i,0);
            RHSUW(ip,1)=(nprec)0.0;
            RHSUW(ip,2)=(nprec)0.0;
        };
    };

    if (geo.slipnodes.GetRows()>0)
    {
        nprec aux1, aux2, nx, ny;
        unsigned int ip;
        for (unsigned int i=0; i<geo.slipnodes.GetRows(); i++)
        {
            ip=geo.slipnodes.Read(i,0);
            nx=geo.slipnormal.Read(i,0);
            ny=geo.slipnormal.Read(i,1);
            if (!((std::abs(nx)<1e-9)&&(std::abs(ny)<1e-9)))
            {
                aux1=unknowns.Read(ip,1)*nx+unknowns.Read(ip,2)*ny;
                aux2=RHSUW(ip,1)*nx+RHSUW(ip,2)*ny;
                RHSUW(ip,1)-=(aux2+aux1)*nx;
                RHSUW(ip,2)-=(aux2+aux1)*ny;
            } else
            {
                RHSUW(ip,1)=(nprec)0.0;
                RHSUW(ip,2)=(nprec)0.0;
            };

        };
    };

    if (geo.opennodes.GetRows()>0)
    {
        nprec aux1, aux2, nx, ny;
        unsigned int ip;
        for (unsigned int i=0; i<geo.opennodes.GetRows(); i++)
        {
            ip=geo.opennodes.Read(i,0);
            nx=-geo.opennormal.Read(i,1);
            ny=geo.opennormal.Read(i,0);
            if (!((std::abs(nx)<1e-9)&&(std::abs(ny)<1e-9)))
            {
                aux1=unknowns.Read(ip,1)*nx+unknowns.Read(ip,2)*ny;
                aux2=RHSUW(ip,1)*nx+RHSUW(ip,2)*ny;
                RHSUW(ip,1)-=(aux2+aux1)*nx;
                RHSUW(ip,2)-=(aux2+aux1)*ny;
            } else
            {
                RHSUW(ip,1)=(nprec)0.0;
                RHSUW(ip,2)=(nprec)0.0;
            };
        };
    };
    if (geo.nodesu0.GetSize()>0)
    {
        unsigned int ip;
        for (unsigned int i=0; i<geo.nodesu0.GetSize(); i++)
        {
            ip=geo.nodesu0.Read(i);
            RHSUW(ip,1)=(nprec)0.0;
            RHSUW(ip,2)=(nprec)0.0;
        };
    };
    */
};
