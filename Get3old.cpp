#include "Globals.h"
#include "geometry.h"
#include "parameter.h"
#include "AuxVar.h"
#include "myvector.h"
#include "matrix.h"

void Get3old(const geometry& geo, const parameter& par, AuxVar& av, const matrix<nprec>& SFD, matrix<nprec>& unknowns)
{
    av.u2a.Initialize((nprec)0.0);
    av.v2a.Initialize((nprec)0.0);
    nprec zero=(nprec)0.0;
    nprec two=(nprec)2.0;
    nprec three=(nprec)3.0;
    Vector<nprec> RG(geo.nd.GetNnodes());
    RG=(vectorize(unknowns,0,"column")*((nprec)1.0-par.theta3))+((vectorize(av.RHS,0,"column"))*par.theta2);
    Vector<nprec> debx(geo.elem.GetNelem());
    debx.Initialize(zero);
    Vector<nprec> deby(geo.elem.GetNelem());
    deby.Initialize(zero);
    nprec um,vm;
    unsigned int ip,n1,n2,n3,ie,nx,ny;
    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        n1=geo.elem.ReadNode(i,0); n2=geo.elem.ReadNode(i,1); n3=geo.elem.ReadNode(i,2);
        um=zero; vm=zero;
        for (unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            um+=unknowns.Read(ip,1);
            vm+=unknowns.Read(ip,2);
        };
        um/=three;
        vm/=three;

        if (!av.dpdn0(i))
        {
            debx(i)=SFD.Read(i,0)*RG(n1)+SFD.Read(i,1)*RG(n2)+SFD.Read(i,2)*RG(n3);
            deby(i)=SFD.Read(i,3)*RG(n1)+SFD.Read(i,4)*RG(n2)+SFD.Read(i,5)*RG(n3);
        }else
        {
            if(std::abs(av.normalel.Read(i,0))<1e-9)
            {
                debx(i)=SFD.Read(i,0)*RG(n1)+SFD.Read(i,1)*RG(n2)+SFD.Read(i,2)*RG(n3);
                deby(i)=(nprec)0.0;
            }else if (std::abs(av.normalel.Read(i,1))<1e-9)
            {
                debx(i)=(nprec)0.0;
                deby(i)=SFD.Read(i,3)*RG(n1)+SFD.Read(i,4)*RG(n2)+SFD.Read(i,5)*RG(n3);
            }else
            {
                debx(i)=SFD.Read(i,0)*RG(n1)+SFD.Read(i,1)*RG(n2)+SFD.Read(i,2)*RG(n3);
                deby(i)=-debx(i)*av.normalel.Read(i,0)/av.normalel.Read(i,1);
            };
        };

        //Perform RHS assembling
        for (unsigned int j=0; j<3; j++)
        {
            ip=geo.elem.ReadNode(i,j);
            av.u2a(ip)+=(-(debx(i))*SFD.Read(i,6)/((nprec)6.0))/av.rhomel(i);
            av.v2a(ip)+=(-(deby(i))*SFD.Read(i,6)/((nprec)6.0))/av.rhomel(i);
            av.u2a(ip)+=(-(par.DELTP/two)*(um*SFD.Read(i,j)+vm*SFD.Read(i,j+3))*debx(i)*SFD.Read(i,6)/two)/av.rhomel(i);
            av.v2a(ip)+=(-(par.DELTP/two)*(um*SFD.Read(i,j)+vm*SFD.Read(i,j+3))*deby(i)*SFD.Read(i,6)/two)/av.rhomel(i);
        };
    };

    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
	if(geo.bd.ReadBC(i)!=13&&geo.bd.ReadBC(i)!=3&&geo.bd.ReadBC(i)!=20&&geo.bd.ReadBC(i)!=2)       //It seem that dpdn must be 0 when u is imposed^M
        {
        n1=geo.bd.ReadNode(i,0); n2=geo.bd.ReadNode(i,1);
        nx=geo.bd.ReadNormal(i,0); ny=geo.bd.ReadNormal(i,1);
        ie=geo.bd.ReadEO(i);
        um=unknowns.Read(n1,1)+unknowns.Read(n2,1);
        vm=unknowns.Read(n1,2)+unknowns.Read(n2,2);
        av.u2a(n1)+=((par.DELTP/two)*(nx*(um+unknowns.Read(n1,1))+ny*(vm+unknowns.Read(n1,2)))*debx(ie)*geo.bd.ReadLen(i)/((nprec)6.0))/av.rhomel(ie);
        av.v2a(n1)+=((par.DELTP/two)*(nx*(um+unknowns.Read(n1,1))+ny*(vm+unknowns.Read(n1,2)))*deby(ie)*geo.bd.ReadLen(i)/((nprec)6.0))/av.rhomel(ie);
        av.u2a(n2)+=((par.DELTP/two)*(nx*(um+unknowns.Read(n2,1))+ny*(vm+unknowns.Read(n2,2)))*debx(ie)*geo.bd.ReadLen(i)/((nprec)6.0))/av.rhomel(ie);
        av.v2a(n2)+=((par.DELTP/two)*(nx*(um+unknowns.Read(n2,1))+ny*(vm+unknowns.Read(n2,2)))*deby(ie)*geo.bd.ReadLen(i)/((nprec)6.0))/av.rhomel(ie);
	};
    };

    av.u2a*=par.DELTP;
    av.v2a*=par.DELTP;

    //Solver
    Eigen::Matrix<nprec,Eigen::Dynamic,1> b(geo.nd.GetNnodes());
    Eigen::Matrix<nprec,Eigen::Dynamic,1> x(geo.nd.GetNnodes());

    for (unsigned int i=0; i<av.FRHS.GetSize(); i++)
    {
        b(i)=av.u2a.Read(i);
    };

    x=av.solverCMM.solve(b);
    for (unsigned int i=0; i<av.FRHS.GetSize(); i++)
    {
        av.u2a(i)=x(i);
    };

    for (unsigned int i=0; i<av.FRHS.GetSize(); i++)
    {
        b(i)=av.v2a.Read(i);
    };

    x=av.solverCMM.solve(b);
    for (unsigned int i=0; i<av.FRHS.GetSize(); i++)
    {
        av.v2a(i)=x(i);
    };
}
