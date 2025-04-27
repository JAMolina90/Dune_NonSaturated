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

void ConvTerm(const geometry& geo, const parameter& par, AuxVar& av, const matrix<nprec>& SFD, const matrix<nprec>& unknowns, AuxVarGeo& avG)
{
    nprec area;
    unsigned int nelem = geo.elem.GetNelem();
    Vector<nprec> fxsec(nelem);
    Vector<nprec> fysec(nelem);
    Vector<nprec> uxnvec(nelem);
    Vector<nprec> uynvec(nelem);
    Vector<nprec> vxnvec(nelem);
    Vector<nprec> vynvec(nelem);
    matrix<nprec> fxel(3,2);
    matrix<nprec> fyel(3,2);
    matrix<nprec> velocn(3,2);
    matrix<nprec> velocn12(3,2);
    matrix<nprec> unkel(3,3);
    Vector<nprec> geome(6);
    unsigned int ip,n1,n2,n3;
    nprec ct2,ct3,uxn,uyn,vxn,vyn,uxn12,vxn12,vyn12,uen12,ven12,uen,ven,soc2,soc3,vau;
    av.ConvTerm.Initialize((nprec)0.0);
    for (unsigned int i=0; i<nelem; i++)
    {
        area = SFD.Read(i,6)/((nprec)2.0);
        fxsec(i) = (nprec)0.0;
        fysec(i) = (nprec)0.0;
        n1=geo.elem.ReadNode(i,0); n2=geo.elem.ReadNode(i,1); n3=geo.elem.ReadNode(i,2);
        for (unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            geome(j) = SFD.Read(i,j);
            geome(j+3) = SFD.Read(i,j+3);
            velocn(j,0) = unknowns.Read(ip,1);
            velocn(j,1) = unknowns.Read(ip,2);

            velocn12(j,0) = (unknowns.Read(ip,1)+av.VelPred.Read(ip,0))/((nprec)2.0);
            velocn12(j,1) = (unknowns.Read(ip,2)+av.VelPred.Read(ip,1))/((nprec)2.0);

            unkel(j,0) = unknowns.Read(ip,0);
            unkel(j,1) = unknowns.Read(ip,1);
            unkel(j,2) = unknowns.Read(ip,2);
            fxel(j,0) = velocn12(j,0)*unkel(j,1);
            fxel(j,1) = velocn12(j,0)*unkel(j,2);
            fyel(j,0) = velocn12(j,1)*unkel(j,1);
            fyel(j,1) = velocn12(j,1)*unkel(j,2);
        };
        if (!av.dudn0(i))
        {
            ct2=geome(0)*fxel(0,0)+geome(1)*fxel(1,0)+geome(2)*fxel(2,0)+\
                geome(3)*fxel(0,1)+geome(4)*fxel(1,1)+geome(5)*fxel(2,1);
            ct3=geome(0)*fyel(0,0)+geome(1)*fyel(1,0)+geome(2)*fyel(2,0)+\
                geome(3)*fyel(0,1)+geome(4)*fyel(1,1)+geome(5)*fyel(2,1);

            uxn=geome(0)*velocn(0,0)+geome(1)*velocn(1,0)+geome(2)*velocn(2,0);
            uyn=geome(3)*velocn(0,0)+geome(4)*velocn(1,0)+geome(5)*velocn(2,0);
            vxn=geome(0)*velocn(0,1)+geome(1)*velocn(1,1)+geome(2)*velocn(2,1);
            vyn=geome(3)*velocn(0,1)+geome(4)*velocn(1,1)+geome(5)*velocn(2,1);
            uxn12=geome(0)*velocn12(0,0)+geome(1)*velocn12(1,0)+geome(2)*velocn12(2,0);
            //uyn12=geome(3)*velocn12(0,0)+geome(4)*velocn12(1,0)+geome(5)*velocn12(2,0);
            vxn12=geome(0)*velocn12(0,1)+geome(1)*velocn12(1,1)+geome(2)*velocn12(2,1);
            vyn12=geome(3)*velocn12(0,1)+geome(4)*velocn12(1,1)+geome(5)*velocn12(2,1);
        }else
        {
            if(std::abs(av.normalel.Read(i,0))<1e-9)
            {
                uxn12=geome(0)*velocn12(0,0)+geome(1)*velocn12(1,0)+geome(2)*velocn12(2,0);
                //uyn12=(nprec)0.0;
                vxn12=geome(0)*velocn12(0,1)+geome(1)*velocn12(1,1)+geome(2)*velocn12(2,1);
                vyn12=(nprec)0.0;
                uxn=geome(0)*velocn(0,0)+geome(1)*velocn(1,0)+geome(2)*velocn(2,0);
                uyn=(nprec)0.0;
                vxn=geome(0)*velocn(0,1)+geome(1)*velocn(1,1)+geome(2)*velocn(2,1);
                vyn=(nprec)0.0;
            }else if (std::abs(av.normalel.Read(i,1))<1e-9)
            {
                uxn12=(nprec)0.0;
                //uyn12=geome(3)*velocn12(0,0)+geome(4)*velocn12(1,0)+geome(5)*velocn12(2,0);
                vxn12=(nprec)0.0;
                vyn12=geome(3)*velocn12(0,1)+geome(4)*velocn12(1,1)+geome(5)*velocn12(2,1);
                uxn=(nprec)0.0;
                uyn=geome(3)*velocn(0,0)+geome(4)*velocn(1,0)+geome(5)*velocn(2,0);
                vxn=(nprec)0.0;
                vyn=geome(3)*velocn(0,1)+geome(4)*velocn(1,1)+geome(5)*velocn(2,1);
            }else
            {
                uxn12=geome(0)*velocn12(0,0)+geome(1)*velocn12(1,0)+geome(2)*velocn12(2,0);
                //uyn12=-uxn12*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                vxn12=geome(0)*velocn12(0,1)+geome(1)*velocn12(1,1)+geome(2)*velocn12(2,1);
                vyn12=-vxn12*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                uxn=geome(0)*velocn(0,0)+geome(1)*velocn(1,0)+geome(2)*velocn(2,0);
                uyn=-uxn*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                vxn=geome(0)*velocn(0,1)+geome(1)*velocn(1,1)+geome(2)*velocn(2,1);
                vyn=-vxn*av.normalel.Read(i,0)/av.normalel.Read(i,1);
            };

            ct2=uxn*(velocn12(0,0)+velocn12(1,0)+velocn12(2,0))/((nprec)3.0)+\
                uyn*(velocn12(0,1)+velocn12(1,1)+velocn12(2,1))/((nprec)3.0)+\
                (uxn12+vyn12)*(velocn(0,0)+velocn(1,0)+velocn(2,0))/((nprec)3.0);
            ct3=vxn*(velocn12(0,0)+velocn12(1,0)+velocn12(2,0))/((nprec)3.0)+\
                vyn*(velocn12(0,1)+velocn12(1,1)+velocn12(2,1))/((nprec)3.0)+\
                (uxn12+vyn12)*(velocn(0,1)+velocn(1,1)+velocn(2,1))/((nprec)3.0);
        };

        uxnvec(i)=uxn;
        uynvec(i)=uyn;
        vxnvec(i)=vxn;
        vynvec(i)=vyn;

        for (unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            av.RHS(ip,1)+=-ct2*area/((nprec)3.0);
            av.RHS(ip,2)+=-ct3*area/((nprec)3.0);
        };

        //Divergence term if compressibility is activated
        if (par.c0!=(nprec)0.0&&par.c1!=(nprec)0.0)
        {
            av.RHS(n1,1)+=(uxn12+vyn12)*((nprec)2.0*velocn12(0,0)+velocn12(1,0)+velocn12(2,0))*area/((nprec)12.0);
            av.RHS(n2,1)+=(uxn12+vyn12)*(velocn12(0,0)+(nprec)2.0*velocn12(1,0)+velocn12(2,0))*area/((nprec)12.0);
            av.RHS(n3,1)+=(uxn12+vyn12)*(velocn12(0,0)+velocn12(1,0)+(nprec)2.0*velocn12(2,0))*area/((nprec)12.0);
            av.RHS(n1,2)+=(uxn12+vyn12)*((nprec)2.0*velocn12(0,1)+velocn12(1,1)+velocn12(2,1))*area/((nprec)12.0);
            av.RHS(n2,2)+=(uxn12+vyn12)*(velocn12(0,1)+(nprec)2.0*velocn12(1,1)+velocn12(2,1))*area/((nprec)12.0);
            av.RHS(n3,2)+=(uxn12+vyn12)*(velocn12(0,1)+velocn12(1,1)+(nprec)2.0*velocn12(2,1))*area/((nprec)12.0);
        };

        //Second order terms
        uen12=(velocn12(0,0)+velocn12(1,0)+velocn12(2,0))/((nprec)3.0);
        ven12=(velocn12(0,1)+velocn12(1,1)+velocn12(2,1))/((nprec)3.0);
        uen=(velocn(0,0)+velocn(1,0)+velocn(2,0))/((nprec)3.0);
        ven=(velocn(0,1)+velocn(1,1)+velocn(2,1))/((nprec)3.0);

        for(unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            av.RHS(ip,1)+=par.DELTP/((nprec)2.0)*(uen12*uxn*uxn+uen12*vxn*uyn+ven12*uyn*uxn+ven12*vyn*uyn)*area/((nprec)3.0);
            av.RHS(ip,2)+=par.DELTP/((nprec)2.0)*(uen12*uxn*vxn+uen12*vxn*vyn+ven12*uyn*vxn+ven12*vyn*vyn)*area/((nprec)3.0);
        };

        soc2=(nprec)0.0;
        soc3=(nprec)0.0;
        if (!(par.c0!=(nprec)0.0&&par.c1!=(nprec)0.0))
        {
            soc2+=par.DELTP/((nprec)2.0)*uen*(uxn+vyn)*area;
            soc3+=par.DELTP/((nprec)2.0)*ven*(uxn+vyn)*area;
        };
        for(unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            vau=geome(j)*uen12+geome(j+3)*ven12;
            av.RHS(ip,1)+=-soc2*vau;
            av.RHS(ip,2)+=-soc3*vau;
        };
        fxsec(i) = soc2*((nprec)2.0)/area;
        fysec(i) = soc3*((nprec)2.0)/area;

        for (unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            av.RHS(ip,1)-=par.DELTP/((nprec)2.0)*area*(geome(j)*uen12*(uen12*uxn+ven12*uyn)+geome(j+3)*ven12*(uen12*uxn+ven12*uyn));
            av.RHS(ip,2)-=par.DELTP/((nprec)2.0)*area*(geome(j)*uen12*(uen12*vxn+ven12*vyn)+geome(j+3)*ven12*(uen12*vxn+ven12*vyn));
        };
    };

    //Second order flux boundary terms
    int i0,i1,ie;
    velocn12.Initialize((nprec)0.0);
    velocn.Initialize((nprec)0.0);
    nprec aleng,anx,any,uma,vma;
    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        aleng=geo.bd.ReadLen(i)/((nprec)12.0);
        anx=geo.bd.ReadNormal(i,0)*aleng;
        any=geo.bd.ReadNormal(i,1)*aleng;
        i0=geo.bd.ReadNode(i,0);
        i1=geo.bd.ReadNode(i,1);
        velocn(0,0)=unknowns.Read(i0,1);
        velocn(1,0)=unknowns.Read(i1,1);
        velocn(0,1)=unknowns.Read(i0,2);
        velocn(1,1)=unknowns.Read(i1,2);
        velocn12(0,0) = (unknowns.Read(i0,1)+av.VelPred.Read(i0,0))/((nprec)2.0);
        velocn12(0,1) = (unknowns.Read(i0,2)+av.VelPred.Read(i0,1))/((nprec)2.0);
        velocn12(1,0) = (unknowns.Read(i1,1)+av.VelPred.Read(i1,0))/((nprec)2.0);
        velocn12(1,1) = (unknowns.Read(i1,2)+av.VelPred.Read(i1,1))/((nprec)2.0);
        ie=geo.bd.ReadEO(i);
        uma=velocn12(0,0)+velocn12(1,0);
        vma=velocn12(0,1)+velocn12(1,1);

        av.RHS(i0,1)+=fxsec(ie)*(anx*(uma+velocn12(0,0))+any*(vma+velocn12(0,1)));
        av.RHS(i1,1)+=fxsec(ie)*(anx*(uma+velocn12(1,0))+any*(vma+velocn12(1,1)));
        av.RHS(i0,2)+=fysec(ie)*(anx*(uma+velocn12(0,0))+any*(vma+velocn12(0,1)));
        av.RHS(i1,2)+=fysec(ie)*(anx*(uma+velocn12(1,0))+any*(vma+velocn12(1,1)));

        av.RHS(i0,1)+=par.DELTP/((nprec)2.0)*(uma/((nprec)2.0)*geo.bd.ReadNormal(i,0)+vma/((nprec)2.0)*geo.bd.ReadNormal(i,1))*(uma/((nprec)2.0)*uxnvec(ie)+vma/((nprec)2.0)*uynvec(ie))*geo.bd.ReadLen(i)/((nprec)2.0);
        av.RHS(i1,1)+=par.DELTP/((nprec)2.0)*(uma/((nprec)2.0)*geo.bd.ReadNormal(i,0)+vma/((nprec)2.0)*geo.bd.ReadNormal(i,1))*(uma/((nprec)2.0)*uxnvec(ie)+vma/((nprec)2.0)*uynvec(ie))*geo.bd.ReadLen(i)/((nprec)2.0);
        av.RHS(i0,2)+=par.DELTP/((nprec)2.0)*(uma/((nprec)2.0)*geo.bd.ReadNormal(i,0)+vma/((nprec)2.0)*geo.bd.ReadNormal(i,1))*(uma/((nprec)2.0)*vxnvec(ie)+vma/((nprec)2.0)*vynvec(ie))*geo.bd.ReadLen(i)/((nprec)2.0);
        av.RHS(i1,2)+=par.DELTP/((nprec)2.0)*(uma/((nprec)2.0)*geo.bd.ReadNormal(i,0)+vma/((nprec)2.0)*geo.bd.ReadNormal(i,1))*(uma/((nprec)2.0)*vxnvec(ie)+vma/((nprec)2.0)*vynvec(ie))*geo.bd.ReadLen(i)/((nprec)2.0);
    };

    for(unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        av.RHS(i,1)*=par.DELTP;
        av.RHS(i,2)*=par.DELTP;
    };

    //Solver
    Eigen::Matrix<nprec,Eigen::Dynamic,1> b(geo.nd.GetNnodes());
    Eigen::Matrix<nprec,Eigen::Dynamic,1> x(geo.nd.GetNnodes());

    for (unsigned int i=0; i<av.FRHS.GetSize(); i++)
    {
        b(i)=av.RHS.Read(i,1);
    };

    x=av.solverCMM.solve(b);

    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        av.ConvTerm(i,0)=x(i);
    };

    for (unsigned int i=0; i<av.FRHS.GetSize(); i++)
    {
        b(i)=av.RHS.Read(i,2);
    };

    x=av.solverCMM.solve(b);

    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        av.ConvTerm(i,1)=x(i);
    };

    av.RHS.Initialize((nprec)0.0);

};

void FirstStep(const geometry& geo, const parameter& par, AuxVar& av, const matrix<nprec>& SFD, const matrix<nprec>& unknowns, AuxVarGeo& avG, const int& ind)
{

    nprec area;
    unsigned int nelem = geo.elem.GetNelem();
    Vector<nprec> fxsec(nelem);
    Vector<nprec> fysec(nelem);
    Vector<nprec> tauxx(nelem);
    Vector<nprec> tauxy(nelem);
    Vector<nprec> tauyy(nelem);
    Vector<nprec> uxnvec(nelem);
    Vector<nprec> uynvec(nelem);
    Vector<nprec> vxnvec(nelem);
    Vector<nprec> vynvec(nelem);
    matrix<nprec> fxel(3,2);
    matrix<nprec> fyel(3,2);
    matrix<nprec> velocn(3,2);
    matrix<nprec> velocn12(3,2);
    matrix<nprec> unkel(3,3);
    Vector<nprec> geome(6);
    Vector<nprec> dviscoel(3);
    unsigned int ip,n1,n2,n3;
    nprec uast,tauc,ct2,ct3,uxn,uyn,vxn,vyn,uxn12,uyn12,vxn12,vyn12,uen12,ven12,uen,ven,soc2,soc3,vau,debx,deby;

    for (unsigned int i=0; i<nelem; i++)
    {
        area = SFD.Read(i,6)/((nprec)2.0);
        fxsec(i) = (nprec)0.0;
        fysec(i) = (nprec)0.0;
        n1=geo.elem.ReadNode(i,0); n2=geo.elem.ReadNode(i,1); n3=geo.elem.ReadNode(i,2);
        for (unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            geome(j) = SFD.Read(i,j);
            geome(j+3) = SFD.Read(i,j+3);
            dviscoel(j) = av.dviscon.Read(ip)+(av.dviscon1.Read(ip)-av.dviscon.Read(ip))*par.theta4;
            velocn(j,0) = unknowns.Read(ip,1);
            velocn(j,1) = unknowns.Read(ip,2);
            if (ind==0)
            {
                velocn12(j,0) = velocn(j,0);
                velocn12(j,1) = velocn(j,1);
            }
            else
            {
                velocn12(j,0) = (unknowns.Read(ip,1)+av.VelPred.Read(ip,0))/((nprec)2.0);
                velocn12(j,1) = (unknowns.Read(ip,2)+av.VelPred.Read(ip,1))/((nprec)2.0);
            };
            unkel(j,0) = unknowns.Read(ip,0);
            unkel(j,1) = unknowns.Read(ip,1);
            unkel(j,2) = unknowns.Read(ip,2);
            fxel(j,0) = velocn12(j,0)*unkel(j,1);
            fxel(j,1) = velocn12(j,0)*unkel(j,2);
            fyel(j,0) = velocn12(j,1)*unkel(j,1);
            fyel(j,1) = velocn12(j,1)*unkel(j,2);
        };
        if (!av.dudn0(i))
        {
            ct2=geome(0)*fxel(0,0)+geome(1)*fxel(1,0)+geome(2)*fxel(2,0)+\
                geome(3)*fxel(0,1)+geome(4)*fxel(1,1)+geome(5)*fxel(2,1);
            ct3=geome(0)*fyel(0,0)+geome(1)*fyel(1,0)+geome(2)*fyel(2,0)+\
                geome(3)*fyel(0,1)+geome(4)*fyel(1,1)+geome(5)*fyel(2,1);

            uxn=geome(0)*velocn(0,0)+geome(1)*velocn(1,0)+geome(2)*velocn(2,0);
            uyn=geome(3)*velocn(0,0)+geome(4)*velocn(1,0)+geome(5)*velocn(2,0);
            vxn=geome(0)*velocn(0,1)+geome(1)*velocn(1,1)+geome(2)*velocn(2,1);
            vyn=geome(3)*velocn(0,1)+geome(4)*velocn(1,1)+geome(5)*velocn(2,1);
            uxn12=geome(0)*velocn12(0,0)+geome(1)*velocn12(1,0)+geome(2)*velocn12(2,0);
            uyn12=geome(3)*velocn12(0,0)+geome(4)*velocn12(1,0)+geome(5)*velocn12(2,0);
            vxn12=geome(0)*velocn12(0,1)+geome(1)*velocn12(1,1)+geome(2)*velocn12(2,1);
            vyn12=geome(3)*velocn12(0,1)+geome(4)*velocn12(1,1)+geome(5)*velocn12(2,1);
        }
        else
        {
            if(std::abs(av.normalel.Read(i,0))<1e-9)
            {
                uxn12=geome(0)*velocn12(0,0)+geome(1)*velocn12(1,0)+geome(2)*velocn12(2,0);
                uyn12=(nprec)0.0;
                vxn12=geome(0)*velocn12(0,1)+geome(1)*velocn12(1,1)+geome(2)*velocn12(2,1);
                vyn12=(nprec)0.0;
                uxn=geome(0)*velocn(0,0)+geome(1)*velocn(1,0)+geome(2)*velocn(2,0);
                uyn=(nprec)0.0;
                vxn=geome(0)*velocn(0,1)+geome(1)*velocn(1,1)+geome(2)*velocn(2,1);
                vyn=(nprec)0.0;
            }else if (std::abs(av.normalel.Read(i,1))<1e-9)
            {
                uxn12=(nprec)0.0;
                uyn12=geome(3)*velocn12(0,0)+geome(4)*velocn12(1,0)+geome(5)*velocn12(2,0);
                vxn12=(nprec)0.0;
                vyn12=geome(3)*velocn12(0,1)+geome(4)*velocn12(1,1)+geome(5)*velocn12(2,1);
                uxn=(nprec)0.0;
                uyn=geome(3)*velocn(0,0)+geome(4)*velocn(1,0)+geome(5)*velocn(2,0);
                vxn=(nprec)0.0;
                vyn=geome(3)*velocn(0,1)+geome(4)*velocn(1,1)+geome(5)*velocn(2,1);
            }else
            {
                uxn12=geome(0)*velocn12(0,0)+geome(1)*velocn12(1,0)+geome(2)*velocn12(2,0);
                uyn12=-uxn12*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                vxn12=geome(0)*velocn12(0,1)+geome(1)*velocn12(1,1)+geome(2)*velocn12(2,1);
                vyn12=-vxn12*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                uxn=geome(0)*velocn(0,0)+geome(1)*velocn(1,0)+geome(2)*velocn(2,0);
                uyn=-uxn*av.normalel.Read(i,0)/av.normalel.Read(i,1);
                vxn=geome(0)*velocn(0,1)+geome(1)*velocn(1,1)+geome(2)*velocn(2,1);
                vyn=-vxn*av.normalel.Read(i,0)/av.normalel.Read(i,1);
            };

            ct2=uxn*(velocn12(0,0)+velocn12(1,0)+velocn12(2,0))/((nprec)3.0)+\
                uyn*(velocn12(0,1)+velocn12(1,1)+velocn12(2,1))/((nprec)3.0)+\
                (uxn12+vyn12)*(velocn(0,0)+velocn(1,0)+velocn(2,0))/((nprec)3.0);
            ct3=vxn*(velocn12(0,0)+velocn12(1,0)+velocn12(2,0))/((nprec)3.0)+\
                vyn*(velocn12(0,1)+velocn12(1,1)+velocn12(2,1))/((nprec)3.0)+\
                (uxn12+vyn12)*(velocn(0,1)+velocn(1,1)+velocn(2,1))/((nprec)3.0);
        };

        uxnvec(i)=uxn;
        uynvec(i)=uyn;
        vxnvec(i)=vxn;
        vynvec(i)=vyn;

        if (!av.dpdn0(i))
        {
            debx=geome(0)*unkel(0,0)+geome(1)*unkel(1,0)+geome(2)*unkel(2,0);
            deby=geome(3)*unkel(0,0)+geome(4)*unkel(1,0)+geome(5)*unkel(2,0);
        }else
        {
            if(std::abs(av.normalel.Read(i,0))<1e-9)
            {
                debx=geome(0)*unkel(0,0)+geome(1)*unkel(1,0)+geome(2)*unkel(2,0);
                deby=(nprec)0.0;
            }else if (std::abs(av.normalel.Read(i,1))<1e-9)
            {
                debx=(nprec)0.0;
                deby=geome(3)*unkel(0,0)+geome(4)*unkel(1,0)+geome(5)*unkel(2,0);
            }else
            {
                debx=geome(0)*unkel(0,0)+geome(1)*unkel(1,0)+geome(2)*unkel(2,0);
                deby=-debx*av.normalel.Read(i,0)/av.normalel.Read(i,1);
            };
        };

        tauc=(dviscoel(0)+dviscoel(1)+dviscoel(2))/((nprec)3.0);
        tauc*=par.Re;
        tauxx(i)=uxn12*(nprec)2.0;
        tauxy(i)=(uyn12+vxn12);
        tauyy(i)=vyn12*(nprec)2.0;
        if (av.walllawind.Read(i))
        {
            if(std::abs(av.normalel.Read(i,0))<1e-9)
            {
                if (av.walllawY.Read(i)>par.z00)
                {
                    uast=std::abs(unknowns.Read(av.walllawnode.Read(i),1))*0.41/std::log(av.walllawY.Read(i)/par.z00);
                    uast=uast*uast*av.rhomel.Read(i);
                    if (uyn>(nprec)0.0)
                    {
                        tauxy(i)=uast;
                    }
                    else
                    {
                        tauxy(i)=-uast;
                    };
                };
            }
            else if (std::abs(av.normalel.Read(i,1))<1e-9)
            {
                if (av.walllawY.Read(i)>par.z00)
                {
                    uast=std::abs(unknowns.Read(av.walllawnode.Read(i),2))*0.41/std::log(av.walllawY.Read(i)/par.z00);
                    uast=uast*uast*av.rhomel.Read(i);
                    if (vxn>(nprec)0.0)
                    {
                        tauxy(i)=uast;
                    }
                    else
                    {
                        tauxy(i)=-uast;
                    };
                };
            }
            else
            {
                std::cout << "Not Contemplated wall log law boundary condition with inclined boundaries" << std::endl;
                exit(1);
            };
        };

        //First order terms contributions to RHS
        if (ind==0)
        {
            ct3-=par.Fr2;
        }
        else
        {
            ct2=(nprec)0.0;
            ct3=-par.Fr2;
        };

        //Surface tension  //Surface tension has + sign because in curvature calculation we must include a - sign
        //ct2+=av.surftenterm.Read(i,0);
        //ct3+=av.surftenterm.Read(i,1);

        //Viscous and Pressure terms
        if (!av.dtaudn0(i))
        {
            for (unsigned int j=0; j<3; j++)
            {
                ip=geo.elem.ReadNode(i,j);
                av.RHS(ip,1)+=-ct2*area/((nprec)3.0)-(geome(j)*tauxx(i)*tauc+geome(j+3)*tauxy(i)*tauc)*area/av.rhomel(i);
                av.RHS(ip,2)+=-ct3*area/((nprec)3.0)-(geome(j)*tauxy(i)*tauc+geome(j+3)*tauyy(i)*tauc)*area/av.rhomel(i);
            };
        } else
        {
            if(std::abs(av.normalel.Read(i,0))<1e-9)
            {
                for (unsigned int j=0; j<3; j++)
                {
                    ip=geo.elem.ReadNode(i,j);
                    av.RHS(ip,1)+=-ct2*area/((nprec)3.0)-(geome(j)*tauxx(i)*tauc)*area/av.rhomel(i);
                    av.RHS(ip,2)+=-ct3*area/((nprec)3.0)-(geome(j)*tauxy(i)*tauc)*area/av.rhomel(i);
                };
            }else if (std::abs(av.normalel.Read(i,1))<1e-9)
            {
                for (unsigned int j=0; j<3; j++)
                {
                    ip=geo.elem.ReadNode(i,j);
                    av.RHS(ip,1)+=-ct2*area/((nprec)3.0)-(geome(j+3)*tauxy(i)*tauc)*area/av.rhomel(i);
                    av.RHS(ip,2)+=-ct3*area/((nprec)3.0)-(geome(j+3)*tauyy(i)*tauc)*area/av.rhomel(i);
                };
            }else
            {
                std::cout << "Not Contemplated dtaudn=0 condition with inclined boundaries" << std::endl;
                exit(1);
            };
        };


        if (ind==1)
        {
            for (unsigned int j=0; j<3; j++)
            {
                ip=geo.elem.ReadNode(i,j);
                av.RHS(ip,1)+=-par.theta3*debx*(area/((nprec)3.0))/av.rhomel(i);
                av.RHS(ip,2)+=-par.theta3*deby*(area/((nprec)3.0))/av.rhomel(i);
            };
        }
        else
        {
            for (unsigned int j=0; j<3; j++)
            {
                ip=geo.elem.ReadNode(i,j);
                av.RHS(ip,1)+=-debx*(area/((nprec)3.0))/av.rhomel(i);
                av.RHS(ip,2)+=-deby*(area/((nprec)3.0))/av.rhomel(i);
            };
        };

        //Divergence term if compressibility is activated
        if (ind==0)
        {
        if (par.c0!=(nprec)0.0&&par.c1!=(nprec)0.0)
        {
            av.RHS(n1,1)+=(uxn12+vyn12)*((nprec)2.0*velocn12(0,0)+velocn12(1,0)+velocn12(2,0))*area/((nprec)12.0);
            av.RHS(n2,1)+=(uxn12+vyn12)*(velocn12(0,0)+(nprec)2.0*velocn12(1,0)+velocn12(2,0))*area/((nprec)12.0);
            av.RHS(n3,1)+=(uxn12+vyn12)*(velocn12(0,0)+velocn12(1,0)+(nprec)2.0*velocn12(2,0))*area/((nprec)12.0);
            av.RHS(n1,2)+=(uxn12+vyn12)*((nprec)2.0*velocn12(0,1)+velocn12(1,1)+velocn12(2,1))*area/((nprec)12.0);
            av.RHS(n2,2)+=(uxn12+vyn12)*(velocn12(0,1)+(nprec)2.0*velocn12(1,1)+velocn12(2,1))*area/((nprec)12.0);
            av.RHS(n3,2)+=(uxn12+vyn12)*(velocn12(0,1)+velocn12(1,1)+(nprec)2.0*velocn12(2,1))*area/((nprec)12.0);
        };
        };

        //Second order terms
        uen12=(velocn12(0,0)+velocn12(1,0)+velocn12(2,0))/((nprec)3.0);
        ven12=(velocn12(0,1)+velocn12(1,1)+velocn12(2,1))/((nprec)3.0);
        uen=(velocn(0,0)+velocn(1,0)+velocn(2,0))/((nprec)3.0);
        ven=(velocn(0,1)+velocn(1,1)+velocn(2,1))/((nprec)3.0);

        if (ind==0)
        {
        for(unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            av.RHS(ip,1)+=par.DELTP/((nprec)2.0)*(uen12*uxn*uxn+uen12*vxn*uyn+ven12*uyn*uxn+ven12*vyn*uyn)*area/((nprec)3.0);
            av.RHS(ip,2)+=par.DELTP/((nprec)2.0)*(uen12*uxn*vxn+uen12*vxn*vyn+ven12*uyn*vxn+ven12*vyn*vyn)*area/((nprec)3.0);
        };
        };

        soc2=(nprec)0.0;
        soc3=(nprec)0.0;
        if (ind==1)
        {
            soc2+=par.DELTP/((nprec)2.0)*(par.theta3*debx/av.rhomel(i))*area;
            soc3+=par.DELTP/((nprec)2.0)*(par.theta3*deby/av.rhomel(i))*area;
        }
        else
        {
            soc2+=par.DELTP/((nprec)2.0)*(debx/av.rhomel(i))*area;
            soc3+=par.DELTP/((nprec)2.0)*(deby/av.rhomel(i))*area;
        };

        if (ind==0)
        {
        if (!(par.c0!=(nprec)0.0&&par.c1!=(nprec)0.0))
        {
            soc2+=par.DELTP/((nprec)2.0)*uen*(uxn+vyn)*area;
            soc3+=par.DELTP/((nprec)2.0)*ven*(uxn+vyn)*area;
        };
        };
        soc3-=par.DELTP/((nprec)2.0)*par.Fr2*area;
        //soc2+=par.DELTP/((nprec)2.0)*av.surftenterm.Read(i,0)*area;
        //soc3+=par.DELTP/((nprec)2.0)*av.surftenterm.Read(i,1)*area;
        for(unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            vau=geome(j)*uen12+geome(j+3)*ven12;
            av.RHS(ip,1)+=-soc2*vau;
            av.RHS(ip,2)+=-soc3*vau;
        };
        fxsec(i) = soc2*((nprec)2.0)/area;
        fysec(i) = soc3*((nprec)2.0)/area;

        if (ind==0)
        {
        for (unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            av.RHS(ip,1)-=par.DELTP/((nprec)2.0)*area*(geome(j)*uen12*(uen12*uxn+ven12*uyn)+geome(j+3)*ven12*(uen12*uxn+ven12*uyn));
            av.RHS(ip,2)-=par.DELTP/((nprec)2.0)*area*(geome(j)*uen12*(uen12*vxn+ven12*vyn)+geome(j+3)*ven12*(uen12*vxn+ven12*vyn));
        };
        };
    };

    //Second order flux boundary terms
    int i0,i1,ie;
    velocn12.Initialize((nprec)0.0);
    velocn.Initialize((nprec)0.0);
    nprec aleng,anx,any,uma,vma,asid,bsid;
    nprec tauxxmod,tauyymod,tauxymod,aux;
    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        aleng=geo.bd.ReadLen(i)/((nprec)12.0);
        anx=geo.bd.ReadNormal(i,0)*aleng;
        any=geo.bd.ReadNormal(i,1)*aleng;
        i0=geo.bd.ReadNode(i,0);
        i1=geo.bd.ReadNode(i,1);
        dviscoel(0)=(av.dviscon.Read(i0)+(av.dviscon1.Read(i0)-av.dviscon.Read(i0))*par.theta4);
        dviscoel(1)=(av.dviscon.Read(i1)+(av.dviscon1.Read(i1)-av.dviscon.Read(i1))*par.theta4);
        velocn(0,0)=unknowns.Read(i0,1);
        velocn(1,0)=unknowns.Read(i1,1);
        velocn(0,1)=unknowns.Read(i0,2);
        velocn(1,1)=unknowns.Read(i1,2);
        if (ind==0)
        {
            velocn12(0,0) = velocn(0,0);
            velocn12(0,1) = velocn(0,1);
            velocn12(1,0) = velocn(1,0);
            velocn12(1,1) = velocn(1,1);
        }
        else
        {
            velocn12(0,0) = (unknowns.Read(i0,1)+av.VelPred.Read(i0,0))/((nprec)2.0);
            velocn12(0,1) = (unknowns.Read(i0,2)+av.VelPred.Read(i0,1))/((nprec)2.0);
            velocn12(1,0) = (unknowns.Read(i1,1)+av.VelPred.Read(i1,0))/((nprec)2.0);
            velocn12(1,1) = (unknowns.Read(i1,2)+av.VelPred.Read(i1,1))/((nprec)2.0);
        };
        ie=geo.bd.ReadEO(i);
        uma=velocn12(0,0)+velocn12(1,0);
        vma=velocn12(0,1)+velocn12(1,1);
        tauc=(dviscoel(0)+dviscoel(1))/((nprec)2.0);
        tauc*=par.Re;
        tauc/=av.rhomel(ie);

        if (ind==0)
        {
        av.RHS(i0,1)+=fxsec(ie)*(anx*(uma+velocn12(0,0))+any*(vma+velocn12(0,1)));
        av.RHS(i1,1)+=fxsec(ie)*(anx*(uma+velocn12(1,0))+any*(vma+velocn12(1,1)));
        av.RHS(i0,2)+=fysec(ie)*(anx*(uma+velocn12(0,0))+any*(vma+velocn12(0,1)));
        av.RHS(i1,2)+=fysec(ie)*(anx*(uma+velocn12(1,0))+any*(vma+velocn12(1,1)));
        };

        //At interface surface tension forces are zero because it was the boundary condition at reinitialization

        if((geo.bd.ReadBC(i)!=1)&&(geo.bd.ReadBC(i)!=3)&&(geo.bd.ReadBC(i)!=5))//&&(geo.bd.ReadBC(i)!=20))
        {
            asid=(nprec)6.0*(anx*tauxx(ie)*tauc+any*tauxy(ie)*tauc);
            bsid=(nprec)6.0*(anx*tauxy(ie)*tauc+any*tauyy(ie)*tauc);
            av.RHS(i0,1)+=asid;
            av.RHS(i0,2)+=bsid;
            av.RHS(i1,1)+=asid;
            av.RHS(i1,2)+=bsid;
        }
        else if (geo.bd.ReadBC(i)==3)
        {
            if (!((std::abs(anx)<std::abs(any)+1e-9)&&(std::abs(anx)>std::abs(any)-1e-9)))
            {
                tauxxmod=tauxx(ie);
                tauyymod=tauyy(ie);
                tauxymod=(anx*any*(tauxxmod-tauyymod))/(anx*anx-any*any);
                asid=(nprec)6.0*(anx*tauxxmod*tauc+any*tauxymod*tauc);
                bsid=(nprec)6.0*(anx*tauxymod*tauc+any*tauyymod*tauc);
                av.RHS(i0,1)+=asid;
                av.RHS(i0,2)+=bsid;
                av.RHS(i1,1)+=asid;
                av.RHS(i1,2)+=bsid;

            } else
            {
                tauxxmod=(tauxx(ie)+tauyy(ie))/((nprec)2.0);
                tauyymod=tauxxmod;
                tauxymod=tauxy(ie);
                asid=(nprec)6.0*(anx*tauxxmod*tauc+any*tauxymod*tauc);
                bsid=(nprec)6.0*(anx*tauxymod*tauc+any*tauyymod*tauc);
                av.RHS(i0,1)+=asid;
                av.RHS(i0,2)+=bsid;
                av.RHS(i1,1)+=asid;
                av.RHS(i1,2)+=bsid;
            };
        }else if (geo.bd.ReadBC(i)==1)
        {
            if (std::abs(geo.bd.ReadNormal(i,0))<1e-9)
            {
                tauxxmod=tauxx(ie);
                tauyymod=0.0;//(unknowns.Read(i0,0)+unknowns.Read(i1,0))/(dviscoel(0)+dviscoel(1));
                tauxymod=tauxy(ie);
                asid=(nprec)6.0*(anx*tauxxmod*tauc+any*tauxymod*tauc);
                bsid=(nprec)6.0*(anx*tauxymod*tauc+any*tauyymod*tauc);
                av.RHS(i0,1)+=asid;
                av.RHS(i0,2)+=bsid;
                av.RHS(i1,1)+=asid;
                av.RHS(i1,2)+=bsid;
            }else if (std::abs(geo.bd.ReadNormal(i,1))<1e-9)
            {
                tauxxmod=0.0;//(unknowns.Read(i0,0)+unknowns.Read(i1,0))/(dviscoel(0)+dviscoel(1));
                tauyymod=tauyy(ie);
                tauxymod=tauxy(ie);
                asid=(nprec)6.0*(anx*tauxxmod*tauc+any*tauxymod*tauc);
                bsid=(nprec)6.0*(anx*tauxymod*tauc+any*tauyymod*tauc);
                av.RHS(i0,1)+=asid;
                av.RHS(i0,2)+=bsid;
                av.RHS(i1,1)+=asid;
                av.RHS(i1,2)+=bsid;
            }else
            {
                tauxxmod=tauxx(ie);
                tauyymod=tauyy(ie);
                tauxymod=(-geo.bd.ReadNormal(i,0)*geo.bd.ReadNormal(i,0)*tauxxmod-geo.bd.ReadNormal(i,1)*geo.bd.ReadNormal(i,1)*tauyymod)/((nprec)2.0*geo.bd.ReadNormal(i,0)*geo.bd.ReadNormal(i,1));
                //Esto es para cuando la presion no la imponia como 0 ((unknowns.Read(i0,0)+unknowns.Read(i1,0))/(dviscoel(0)+dviscoel(1))-geo.bd.ReadNormal(i,0)*geo.bd.ReadNormal(i,0)*tauxxmod-geo.bd.ReadNormal(i,1)*geo.bd.ReadNormal(i,1)*tauyymod)/((nprec)2.0*geo.bd.ReadNormal(i,0)*geo.bd.ReadNormal(i,1));
                asid=(nprec)6.0*(anx*tauxxmod*tauc+any*tauxymod*tauc);
                bsid=(nprec)6.0*(anx*tauxymod*tauc+any*tauyymod*tauc);
                av.RHS(i0,1)+=asid;
                av.RHS(i0,2)+=bsid;
                av.RHS(i1,1)+=asid;
                av.RHS(i1,2)+=bsid;

            };
        } else if (geo.bd.ReadBC(i)==5)
        {
            if (std::abs(geo.bd.ReadNormal(i,0))<1e-9)
            {
                asid=(nprec)6.0*(anx*tauxx(ie)*tauc);
                bsid=(nprec)6.0*(anx*tauxy(ie)*tauc);
                av.RHS(i0,1)+=asid;
                av.RHS(i0,2)+=bsid;
                av.RHS(i1,1)+=asid;
                av.RHS(i1,2)+=bsid;
            }else if (std::abs(geo.bd.ReadNormal(i,1))<1e-9)
            {
                asid=(nprec)6.0*(any*tauxy(ie)*tauc);
                bsid=(nprec)6.0*(any*tauyy(ie)*tauc);
                av.RHS(i0,1)+=asid;
                av.RHS(i0,2)+=bsid;
                av.RHS(i1,1)+=asid;
                av.RHS(i1,2)+=bsid;
            }else
            {
                std::cout << "Not Contemplated dtaudn=0 condition with inclined boundaries" << std::endl;
                exit(1);
            };
        };
        //Extra Forces due to the sediment grains
        if (geo.bd.ReadBC(i)==20)
        {
            aux=-av.TwoWayForce.Read(i0)*par.Re/par.dviscoref*(1.0/av.rhomel(ie));
            asid=(nprec)6.0*(-any*aux);
            bsid=(nprec)6.0*(anx*aux);
            av.RHS(i0,1)+=asid;
            av.RHS(i0,2)+=bsid;
            av.RHS(i1,1)+=asid;
            av.RHS(i1,2)+=bsid;
        };
        /*
        else if (geo.bd.ReadBC(i)==20)
        {
            if (avG.Uast(i)>(nprec)0.0)
            {
                aux=-avG.Uast(i)*avG.Uast(i)*(par.rho0);   //Tension real
            }
            else
            {
                aux=avG.Uast(i)*avG.Uast(i)*(par.rho0);  //Tension real
            };
            aux*=par.Re/par.dviscoref;                    //Dejarlo asi para que las unidades cuadren
            aux/=av.rhomel(ie);
            asid=(nprec)6.0*(-any*aux);
            bsid=(nprec)6.0*(anx*aux);
            std::cout << asid << " " << (nprec)6.0*(anx*tauxx(ie)*tauc+any*tauxy(ie)*tauc) << " ";
            std::cout << bsid << " " << (nprec)6.0*(anx*tauxy(ie)*tauc+any*tauyy(ie)*tauc) << std::endl;
            av.RHS(i0,1)+=asid;
            av.RHS(i0,2)+=bsid;
            av.RHS(i1,1)+=asid;
            av.RHS(i1,2)+=bsid;
        }
        */
        if (ind==0)
        {
        av.RHS(i0,1)+=par.DELTP/((nprec)2.0)*(uma/((nprec)2.0)*geo.bd.ReadNormal(i,0)+vma/((nprec)2.0)*geo.bd.ReadNormal(i,1))*(uma/((nprec)2.0)*uxnvec(ie)+vma/((nprec)2.0)*uynvec(ie))*geo.bd.ReadLen(i)/((nprec)2.0);
        av.RHS(i1,1)+=par.DELTP/((nprec)2.0)*(uma/((nprec)2.0)*geo.bd.ReadNormal(i,0)+vma/((nprec)2.0)*geo.bd.ReadNormal(i,1))*(uma/((nprec)2.0)*uxnvec(ie)+vma/((nprec)2.0)*uynvec(ie))*geo.bd.ReadLen(i)/((nprec)2.0);
        av.RHS(i0,2)+=par.DELTP/((nprec)2.0)*(uma/((nprec)2.0)*geo.bd.ReadNormal(i,0)+vma/((nprec)2.0)*geo.bd.ReadNormal(i,1))*(uma/((nprec)2.0)*vxnvec(ie)+vma/((nprec)2.0)*vynvec(ie))*geo.bd.ReadLen(i)/((nprec)2.0);
        av.RHS(i1,2)+=par.DELTP/((nprec)2.0)*(uma/((nprec)2.0)*geo.bd.ReadNormal(i,0)+vma/((nprec)2.0)*geo.bd.ReadNormal(i,1))*(uma/((nprec)2.0)*vxnvec(ie)+vma/((nprec)2.0)*vynvec(ie))*geo.bd.ReadLen(i)/((nprec)2.0);
        };
    };

    for(unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        av.RHS(i,1)*=par.DELTP;
        av.RHS(i,2)*=par.DELTP;
    };

    //Solver
    Eigen::Matrix<nprec,Eigen::Dynamic,1> b(geo.nd.GetNnodes());
    Eigen::Matrix<nprec,Eigen::Dynamic,1> x(geo.nd.GetNnodes());

    for (unsigned int i=0; i<av.FRHS.GetSize(); i++)
    {
        b(i)=av.RHS.Read(i,1);
    };

    x=av.solverCMM.solve(b);

    if (ind==1)
    {
        for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
        {
            av.RHS(i,1)=x(i)+av.ConvTerm(i,0);
            av.ConvTerm(i,0)=(nprec)0.0;
        };
    }
    else
    {
        for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
        {
            av.VelPred(i,0)=x(i);
        };
    };

    for (unsigned int i=0; i<av.FRHS.GetSize(); i++)
    {
        b(i)=av.RHS.Read(i,2);
    };

    x=av.solverCMM.solve(b);

    if (ind==1)
    {
        for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
        {
            av.RHS(i,2)=x(i)+av.ConvTerm(i,1);
            av.ConvTerm(i,1)=(nprec)0.0;
        };
    }
    else
    {
        for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
        {
            av.VelPred(i,1)=x(i);
        };
    };

    if (ind==0)
    {
        av.RHS.Initialize((nprec)0.0);

        //u prescribed boundary condition correction
        if (par.ksi)
        {
            unsigned int in1,in2;
            for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
            {
                if(geo.bd.ReadBC(i)==13)
                {
                    in1=geo.bd.ReadNode(i,0); in2=geo.bd.ReadNode(i,1);
                    av.VelPred(in1,0)=(nprec)0.0;
                    av.VelPred(in1,1)=(nprec)0.0;
                    av.VelPred(in2,0)=(nprec)0.0;
                    av.VelPred(in2,1)=(nprec)0.0;
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
                av.VelPred(ip,0)=(nprec)0.0;
                av.VelPred(ip,1)=(nprec)0.0;
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
                    aux2=av.VelPred(ip,0)*nx+av.VelPred(ip,1)*ny;
                    av.VelPred(ip,0)-=(aux2+aux1)*nx;
                    av.VelPred(ip,1)-=(aux2+aux1)*ny;
                } else
                {
                    av.VelPred(ip,0)=(nprec)0.0;
                    av.VelPred(ip,1)=(nprec)0.0;
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
                    aux2=av.VelPred(ip,0)*nx+av.VelPred(ip,1)*ny;
                    av.VelPred(ip,0)-=(aux2+aux1)*nx;
                    av.VelPred(ip,1)-=(aux2+aux1)*ny;
                } else
                {
                    av.VelPred(ip,0)=(nprec)0.0;
                    av.VelPred(ip,1)=(nprec)0.0;
                };
            };
        };

        if (geo.nodesu0.GetSize()>0)
        {
            unsigned int ip;
            for (unsigned int i=0; i<geo.nodesu0.GetSize(); i++)
            {
                ip=geo.nodesu0.Read(i);
                av.VelPred(ip,0)=(nprec)0.0;
                av.VelPred(ip,1)=(nprec)0.0;
            };
        };

        for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
        {
            av.VelPred(i,0)+=unknowns.Read(i,1);
            av.VelPred(i,1)+=unknowns.Read(i,2);
        };
    };


    //IMPOSE BOUNDARY CONDITIONS WHEN IND==0 AND ERASE AV.RHS
    /*
    //Velocity correction for wall, slip condition and u prescribed
    if (geo.wnodes.GetRows()>0)
    {
        for (unsigned int i=0; i<geo.wnodes.GetRows(); i++)
        {
            ip=geo.wnodes.Read(i,0);
            av.RHS(ip,1)=(nprec)0.0;
            av.RHS(ip,2)=(nprec)0.0;
        };
    };
    if (geo.slipnodes.GetRows()>0)
    {
        unsigned int ip;
        nprec aux1,aux2, nx, ny;
        for (unsigned int i=0; i<geo.slipnodes.GetRows(); i++)
        {
            ip=geo.slipnodes.Read(i,0);
            nx=geo.slipnormal.Read(i,0);
            ny=geo.slipnormal.Read(i,1);
            if (!((std::abs(nx)<1e-9)&&(std::abs(ny)<1e-9)))
            {
                aux1=unknowns.Read(ip,1)*nx+unknowns.Read(ip,2)*ny;
                aux2=av.RHS(ip,1)*nx+av.RHS(ip,2)*ny;
                av.RHS(ip,1)-=(aux2+aux1)*nx;
                av.RHS(ip,2)-=(aux2+aux1)*ny;
            } else
            {
                av.RHS(ip,1)=(nprec)0.0;
                av.RHS(ip,2)=(nprec)0.0;
            };
        };
    };
    if (geo.opennodes.GetRows()>0)
    {
        unsigned int ip;
        nprec aux1, aux2, nx, ny;
        for (unsigned int i=0; i<geo.opennodes.GetRows(); i++)
        {
            ip=geo.opennodes.Read(i,0);
            nx=-geo.opennormal.Read(i,1);
            ny=geo.opennormal.Read(i,0);
            if (!((std::abs(nx)<1e-9)&&(std::abs(ny)<1e-9)))
            {
                aux1=unknowns.Read(ip,1)*nx+unknowns.Read(ip,2)*ny;
                aux2=av.RHS(ip,1)*nx+av.RHS(ip,2)*ny;
                av.RHS(ip,1)-=(aux2+aux1)*nx;
                av.RHS(ip,2)-=(aux2+aux1)*ny;
            } else
            {
                av.RHS(ip,1)=(nprec)0.0;
                av.RHS(ip,2)=(nprec)0.0;
            };
        };
    };
    if (par.ksi)
    {
        unsigned int n1,n2;
        for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
        {
            if(geo.bd.ReadBC(i)==13)
            {
                n1=geo.bd.ReadNode(i,0); n2=geo.bd.ReadNode(i,1);
                av.RHS(n1,1)=(nprec)0.0;
                av.RHS(n1,2)=(nprec)0.0;
                av.RHS(n2,1)=(nprec)0.0;
                av.RHS(n2,2)=(nprec)0.0;
            };
        };
    };

    if (geo.nodesu0.GetSize()>0)
    {
        unsigned int ip;
        for (unsigned int i=0; i<geo.nodesu0.GetSize(); i++)
        {
            ip=geo.nodesu0.Read(i);
            av.RHS(ip,1)=(nprec)0.0;
            av.RHS(ip,2)=(nprec)0.0;
        };
    };
    */
};
