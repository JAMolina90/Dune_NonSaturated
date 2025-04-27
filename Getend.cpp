#include "Globals.h"
#include "geometry.h"
#include "parameter.h"
#include "AuxVar.h"
#include "myvector.h"
#include "matrix.h"

void Getend(const geometry& geo, const parameter& par, AuxVar& av, const matrix<nprec>& SFD, const matrix<nprec>& unknowns)
{
    //u prescribed boundary condition correction
    if (par.ksi)
    {
        unsigned int in1,in2;
        for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
        {
            if(geo.bd.ReadBC(i)==13)
            {
                in1=geo.bd.ReadNode(i,0); in2=geo.bd.ReadNode(i,1);
                av.u2a(in1)=(nprec)0.0;
                av.u2a(in2)=(nprec)0.0;
                av.v2a(in1)=(nprec)0.0;
                av.v2a(in2)=(nprec)0.0;
                av.RHS(in1,1)=(nprec)0.0;
                av.RHS(in1,2)=(nprec)0.0;
                av.RHS(in2,1)=(nprec)0.0;
                av.RHS(in2,2)=(nprec)0.0;
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
            av.u2a(ip)=(nprec)0.0;
            av.v2a(ip)=(nprec)0.0;
            av.RHS(ip,1)=(nprec)0.0;
            av.RHS(ip,2)=(nprec)0.0;
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
                aux2=av.RHS(ip,1)*nx+av.RHS(ip,2)*ny;
                av.RHS(ip,1)-=(aux2+aux1)*nx;
                av.RHS(ip,2)-=(aux2+aux1)*ny;
                aux2=av.u2a(ip)*nx+av.v2a(ip)*ny;
                av.u2a(ip)-=aux2*nx;
                av.v2a(ip)-=aux2*ny;
            } else
            {
                av.RHS(ip,1)=(nprec)0.0;
                av.RHS(ip,2)=(nprec)0.0;
                av.u2a(ip)=(nprec)0.0;
                av.v2a(ip)=(nprec)0.0;
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
                aux2=av.RHS(ip,1)*nx+av.RHS(ip,2)*ny;
                av.RHS(ip,1)-=(aux2+aux1)*nx;
                av.RHS(ip,2)-=(aux2+aux1)*ny;
                aux2=av.u2a(ip)*nx+av.v2a(ip)*ny;
                av.u2a(ip)-=aux2*nx;
                av.v2a(ip)-=aux2*ny;
            } else
            {
                av.RHS(ip,1)=(nprec)0.0;
                av.RHS(ip,2)=(nprec)0.0;
                av.u2a(ip)=(nprec)0.0;
                av.v2a(ip)=(nprec)0.0;
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
            av.u2a(ip)=(nprec)0.0;
            av.v2a(ip)=(nprec)0.0;
        };
    };

    //sum dU** to RHS
    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        av.RHS(i,1)+=av.u2a(i);
        av.RHS(i,2)+=av.v2a(i);
    };

}
