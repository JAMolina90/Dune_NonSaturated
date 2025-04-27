#include <cmath>
#include <limits>
#include "Globals.h"
#include "geometry.h"
#include "myvector.h"
#include "matrix.h"
#include "AuxVar.h"
#include "parameter.h"
#include "UpwindVar.h"
#include "Neighbour.h"

const double eps=1e-9;

void HvaluesIni(AuxVarGeo& s, const geometry& geo)
{
    nprec x0=75.0;
    nprec a=50.0;
    nprec r;
    nprec hb=0.0;
    nprec h0=7.5; //Altura Cresta
    s.Hmax=std::numeric_limits<nprec>::min();
    nprec X,h;
    //unsigned int nod;
    for (unsigned int i=0; i<s.nodesh.GetSize(); i++)
    {
        X=s.nodeshX(i);
        r=std::sqrt((X-x0)*(X-x0));
        if (r/a>(nprec)1.0)
        {
            h=hb;
            s.Hmax=std::max(h,s.Hmax);
            s.h(i)=h;
        }
        else
        {
            h=hb+h0*std::cos(M_PI*r/(2.0*a))*std::cos(M_PI*r/(2.0*a));
            s.Hmax=std::max(h,s.Hmax);
            s.h(i)=h;
        };
    };
}

void HnodesCalculationIni(AuxVarGeo& s, const geometry& geo)
{
    Vector<bool> check(geo.nd.GetNnodes());
    check.Initialize(false);
    unsigned int cont=0;
    unsigned int cont2=0;

    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        if (geo.bd.ReadBC(i)==20)
        {
            cont2++;
            if (!check(geo.bd.ReadNode(i,0)))
            {
                cont++;
                check(geo.bd.ReadNode(i,0))=true;
            };
            if (!check(geo.bd.ReadNode(i,1)))
            {
                cont++;
                check(geo.bd.ReadNode(i,1))=true;
            };
        };
    };

    s.nodesh.SetDimensions(cont);
    s.nodeshX.SetDimensions(cont);

    s.nodeshqs.SetDimensions(cont);
    s.nodeshMeanqs.SetDimensions(cont);
    s.nodeshveff.SetDimensions(cont);
    s.nodeshuast.SetDimensions(cont);
    s.nodeshsedrho.SetDimensions(cont);
    s.nodeshqs.Initialize((nprec)0.0);
    s.nodeshuast.Initialize((nprec)0.0);
    s.nodeshsedrho.Initialize((nprec)0.0);
    s.nodeshveff.Initialize((nprec)0.0);

    s.h.SetDimensions(cont);
    s.linesh.SetDimensions(cont2,2);
    s.lineslen.SetDimensions(cont2);
    s.linklist.SetDimensions(geo.nd.GetNnodes());
    s.linklist.Initialize(-1);
    s.conectivity.SetDimensions(cont,3);
    s.conectivity.Initialize(0);
    s.nodesbd.SetDimensions(2);

    cont=0;
    cont2=0;
    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        if (check(i))
        {
            s.nodesh(cont)=i;
            s.linklist(i)=cont;
            s.nodeshX(cont)=geo.nd.ReadC(i,0);
            cont++;
        };
    };

    //Sorting
    {
    unsigned int nodeprov,linkprov;
    nprec Xprov;

    for (unsigned int i=1; i<s.nodesh.GetSize(); i++)
    {
        for (unsigned int j=0; j<s.nodesh.GetSize()-i; j++)
        {
            if (s.nodeshX(j)>s.nodeshX(j+1))
            {
                linkprov=s.linklist(s.nodesh(j));
                s.linklist(s.nodesh(j))=s.linklist(s.nodesh(j+1));
                s.linklist(s.nodesh(j+1))=linkprov;
                Xprov = s.nodeshX(j+1);
                nodeprov = s.nodesh(j+1);
                s.nodeshX(j+1) = s.nodeshX(j);
                s.nodesh(j+1) = s.nodesh(j);
                s.nodeshX(j) = Xprov;
                s.nodesh(j) = nodeprov;
            };
        };
    };
    }

    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        if (geo.bd.ReadBC(i)==20)
        {
            if (s.linklist(geo.bd.ReadNode(i,0))>=0)
            {
                s.linesh(cont2,0)=(unsigned int)s.linklist(geo.bd.ReadNode(i,0));
                s.conectivity(s.linesh(cont2,0),0)++;
                if (s.conectivity(s.linesh(cont2,0),0)==1)
                {
                    s.conectivity(s.linesh(cont2,0),1)=cont2;
                }
                else
                {
                    s.conectivity(s.linesh(cont2,0),2)=cont2;
                };
            }
            else
            {
                std::cout << "Error in nodesh numbering" << std::endl;
                exit(1);
            };

            if (s.linklist(geo.bd.ReadNode(i,1))>=0)
            {
                s.linesh(cont2,1)=(unsigned int)s.linklist(geo.bd.ReadNode(i,1));
                s.conectivity(s.linesh(cont2,1),0)++;
                if (s.conectivity(s.linesh(cont2,1),0)==1)
                {
                    s.conectivity(s.linesh(cont2,1),1)=cont2;
                }
                else
                {
                    s.conectivity(s.linesh(cont2,1),2)=cont2;
                };
            }
            else
            {
                std::cout << "Error in nodesh numbering" << std::endl;
                exit(1);
            };

            s.lineslen(cont2)=geo.bd.ReadLen(i);

            cont2++;
        };
    };

    cont=0;
    for (unsigned int i=0; i<s.conectivity.GetRows(); i++)
    {
        if (s.conectivity(i,0)==1)
        {
            s.nodesbd(cont)=i;
            cont++;
        };
    };

    if (s.nodeshX(s.nodesbd(0))>s.nodeshX(s.nodesbd(1)))
    {
        cont=s.nodesbd(1);
        s.nodesbd(1)=s.nodesbd(0);
        s.nodesbd(0)=cont;
    };

    //H values
    HvaluesIni(s,geo);
    //s.nodeshqs.Initialize(1e-5);
    //s.nodeshsedrho.Initialize(1e-5);
};


void UastCalculation(AuxVarGeo& s, const geometry& geo, const matrix<nprec>& unknowns1, parameter& par, const matrix<nprec>& SFD, Neig& NB)
{
    Vector<unsigned int> uastnumber(geo.nd.GetNnodes()); //Pos0 --> number of boundaries for each node
    uastnumber.Initialize(0);
    Vector<nprec> uastboundary(geo.nd.GetNnodes()); //mean qs for each node belonging to a boundary with BC=20
    uastboundary.Initialize((nprec)0.0);

    if (par.SCM==1)
    {
        {
        unsigned int el,nb1,nb2,nfree,el0,elselect,n1,n2,n3;
        int el1;
        nprec nx,ny,ux,uy,uast,aux,z,SFx,SFy,aux1,aux2,aux3;
        nprec HZ=par.z00*par.nz0;
        nprec X,Y;
        bool neigfound;
        for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
        {
            if (geo.bd.ReadBC(i)==20)
            {
                nb1=geo.bd.ReadNode(i,0);
                nb2=geo.bd.ReadNode(i,1);
                nx=-geo.bd.ReadNormal(i,0);
                ny=-geo.bd.ReadNormal(i,1);
                el=geo.bd.ReadEO(i);
                uast=(nprec)0.0;
                SFx=(nprec)0.0;
                SFy=(nprec)0.0;
                for (unsigned int j=0; j<3; j++)
                {
                    nfree=geo.elem.ReadNode(el,j);
                    if ((nfree!=nb1)&&(nfree!=nb2))
                    {
                        SFx=SFD.Read(el,j);
                        SFy=SFD.Read(el,j+3);
                        break;
                    };
                };
                z=(nprec)1.0/(std::sqrt(SFx*SFx+SFy*SFy));
                if (HZ>z)
                {
                    X=(geo.nd.ReadC(nb1,0)+geo.nd.ReadC(nb2,0))/((nprec)2.0)+HZ*nx;
                    Y=(geo.nd.ReadC(nb1,1)+geo.nd.ReadC(nb2,1))/((nprec)2.0)+HZ*ny;
                    el0=s.NeigCorresp(el);
                    neigfound=false;
                    for (unsigned int j=NB.L3_lstE.Read(el0); j<NB.L3_lstE.Read(el0+1); j++)
                    {
                        el1=s.OldElemIndex.Read(NB.L3_elem.Read(j));
                        if (el1>=0)
                        {
                            aux1=SFD.Read((unsigned int)el1,0)*X+SFD.Read((unsigned int)el1,3)*Y+SFD.Read((unsigned int)el1,7);
                            aux2=SFD.Read((unsigned int)el1,1)*X+SFD.Read((unsigned int)el1,4)*Y+SFD.Read((unsigned int)el1,8);
                            aux3=SFD.Read((unsigned int)el1,2)*X+SFD.Read((unsigned int)el1,5)*Y+SFD.Read((unsigned int)el1,9);

                            if ((-eps<aux1)&&(aux1<(nprec)1.0+eps)&&(-eps<aux2)&&(aux2<(nprec)1.0+eps)&&(-eps<aux3)&&(aux3<(nprec)1.0+eps))
                            {
                                elselect=(unsigned int)el1;
                                neigfound=true;
                                break;
                            };
                        };
                    };

                    if (!neigfound)
                    {
                        std::cout << "Neighbor not found" << std::endl;
                        for (unsigned int j=0; j<geo.elem.GetNelem(); j++)
                        {
                            aux1=SFD.Read(j,0)*X+SFD.Read(j,3)*Y+SFD.Read(j,7);
                            aux2=SFD.Read(j,1)*X+SFD.Read(j,4)*Y+SFD.Read(j,8);
                            aux3=SFD.Read(j,2)*X+SFD.Read(j,5)*Y+SFD.Read(j,9);

                            if ((-eps<aux1)&&(aux1<(nprec)1.0+eps)&&(-eps<aux2)&&(aux2<(nprec)1.0+eps)&&(-eps<aux3)&&(aux3<(nprec)1.0+eps))
                            {
                                elselect=j;
                                neigfound=true;
                                break;
                            };
                        };
                    };

                    if (neigfound)
                    {
                        n1=geo.elem.ReadNode(elselect,0);
                        n2=geo.elem.ReadNode(elselect,1);
                        n3=geo.elem.ReadNode(elselect,2);
                        ux=unknowns1.Read(n1,1)*aux1+unknowns1.Read(n2,1)*aux2+unknowns1.Read(n3,1)*aux3;
                        uy=unknowns1.Read(n1,2)*aux1+unknowns1.Read(n2,2)*aux2+unknowns1.Read(n3,2)*aux3;
                        aux=ux*nx+uy*ny;
                        if ((ux-aux*nx)>(nprec)0.0)
                        {
                            uast=std::sqrt((ux-aux*nx)*(ux-aux*nx)+(uy-aux*ny)*(uy-aux*ny));
                        }
                        else
                        {
                            uast=-std::sqrt((ux-aux*nx)*(ux-aux*nx)+(uy-aux*ny)*(uy-aux*ny));
                        };
                        uast*=0.41/std::log(HZ/par.z00);
                    }
                    else
                    {
                        if (z>par.z00)
                        {
                        ux=unknowns1.Read(nfree,1);
                        uy=unknowns1.Read(nfree,2);
                        aux=ux*nx+uy*ny;
                        if ((ux-aux*nx)>(nprec)0.0)
                        {
                            uast=std::sqrt((ux-aux*nx)*(ux-aux*nx)+(uy-aux*ny)*(uy-aux*ny));
                        }
                        else
                        {
                            uast=-std::sqrt((ux-aux*nx)*(ux-aux*nx)+(uy-aux*ny)*(uy-aux*ny));
                        };
                        uast*=0.41/std::log(z/par.z00);
                        };
                    };
                }
                else
                {
                    if (z>par.z00)
                    {
                    ux=unknowns1.Read(nfree,1);
                    uy=unknowns1.Read(nfree,2);
                    aux=ux*nx+uy*ny;
                    if ((ux-aux*nx)>(nprec)0.0)
                    {
                        uast=std::sqrt((ux-aux*nx)*(ux-aux*nx)+(uy-aux*ny)*(uy-aux*ny));
                    }
                    else
                    {
                        uast=-std::sqrt((ux-aux*nx)*(ux-aux*nx)+(uy-aux*ny)*(uy-aux*ny));
                    };
                    uast*=0.41/std::log(z/par.z00);
                    };
                };
                uastboundary(nb1)+=uast;
                uastboundary(nb2)+=uast;
                uastnumber(nb1)++;
                uastnumber(nb2)++;
            };
        };
        }
    }
    else
    {
        unsigned int el,n1,n2,n3,nb1,nb2;
        nprec nx,ny,ux,uy,vx,vy,tauxx,tauyy,tauxy,tw,aux1,aux2,uast;
        for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
        {
            if (geo.bd.ReadBC(i)==20)
            {
                nb1=geo.bd.ReadNode(i,0);
                nb2=geo.bd.ReadNode(i,1);

                nx=geo.bd.ReadNormal(i,0);
                ny=geo.bd.ReadNormal(i,1);
                el=geo.bd.ReadEO(i);

                n1=geo.elem.ReadNode(el,0);
                n2=geo.elem.ReadNode(el,1);
                n3=geo.elem.ReadNode(el,2);

                ux=SFD.Read(el,0)*unknowns1.Read(n1,1)+SFD.Read(el,1)*unknowns1.Read(n2,1)+SFD.Read(el,2)*unknowns1.Read(n3,1);
                uy=SFD.Read(el,3)*unknowns1.Read(n1,1)+SFD.Read(el,4)*unknowns1.Read(n2,1)+SFD.Read(el,5)*unknowns1.Read(n3,1);
                vx=SFD.Read(el,0)*unknowns1.Read(n1,2)+SFD.Read(el,1)*unknowns1.Read(n2,2)+SFD.Read(el,2)*unknowns1.Read(n3,2);
                vy=SFD.Read(el,3)*unknowns1.Read(n1,2)+SFD.Read(el,4)*unknowns1.Read(n2,2)+SFD.Read(el,5)*unknowns1.Read(n3,2);

                tauxx=(nprec)2.0*ux*par.dvisco0;
                tauyy=(nprec)2.0*vy*par.dvisco0;
                tauxy=(uy+vx)*par.dvisco0;

                aux1=tauxx*nx+tauxy*ny;
                aux2=tauxy*nx+tauyy*ny;

                tw=-(-ny*aux1+nx*aux2);
                uast=(nprec)0.0;
                if (tw>=(nprec)0.0)
                {
                    uast=std::sqrt(tw/par.rho0);
                }
                else
                {
                    uast=-std::sqrt(-tw/par.rho0);
                };
                uastboundary(nb1)+=uast;
                uastboundary(nb2)+=uast;
                uastnumber(nb1)++;
                uastnumber(nb2)++;
            };
        };
    };

    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        if (uastnumber(i)>0)
        {
            uastboundary(i)/=(nprec)uastnumber(i);
        };
    };

    //Pass to BD nodes
    s.nodeshuast.Initialize((nprec)0.0);
    {
    unsigned int i,j,nb1,nb2;
    nprec X,x1,x2;
    bool found;
    //#pragma omp parallel for default(none) shared(std::cout,s,geo,unknowns1,par,SFD,qsboundary,qsnumber) private(i,j,X,nb1,nb2,x1,x2,found)
    for (i=0; i<s.nodesh.GetSize(); i++)
    {
        X=s.nodeshX(i);
        found=false;
        for (j=0; j<geo.bd.GetNbounds(); j++)
        {
            if (geo.bd.ReadBC(j)==20)
            {
                nb1=geo.bd.ReadNode(j,0);
                nb2=geo.bd.ReadNode(j,1);
                x1=geo.nd.ReadC(nb1,0);
                x2=geo.nd.ReadC(nb2,0);
                if (x1>x2)
                {
                    nb1=geo.bd.ReadNode(j,1);
                    nb2=geo.bd.ReadNode(j,0);
                    x1=geo.nd.ReadC(nb1,0);
                    x2=geo.nd.ReadC(nb2,0);
                };
                if ((X>x1-eps)&&(X<x2+eps))
                {
                    s.nodeshuast(i)=uastboundary(nb1)+((X-x1)/(x2-x1))*(uastboundary(nb2)-uastboundary(nb1));

                    found=true;
                    break;
                };
            };
        };
        if (!found)
        {
            std::cout << "Correspondence for Node " << i << " not found" << std::endl;
            exit(1);
        };
    };
    }

    //Smoother
    if (par.nsmoothuast!=0)
    {
        Vector<nprec> newsol(s.nodeshuast.GetSize());
        for (unsigned int j=0; j<par.nsmoothuast; j++)
        {
            newsol.Initialize((nprec)0.0);
            for (unsigned int i=0; i<s.nodeshuast.GetSize()-1; i++)
            {
                newsol(i)+=(s.nodeshuast(i)+s.nodeshuast(i+1))/4.0;
                newsol(i+1)+=(s.nodeshuast(i)+s.nodeshuast(i+1))/4.0;
            };
            s.nodeshuast=newsol;
        };
    };
};

void Sediment(AuxVar& av, AuxVarGeo& s, const geometry& geo, const matrix<nprec>& unknowns1, parameter& par, const matrix<nprec>& SFD)
{
    Vector<nprec> Force (s.nodesh.GetSize());
    Force.Initialize((nprec)0.0);
    s.nodeshMeanqs.Initialize((nprec)0.0);

    matrix<nprec> LOsol(s.nodesh.GetSize(),2);
    Vector<nprec> LMM(s.nodesh.GetSize());
    LMM.Initialize((nprec)0.0);
    Vector<nprec> veff(s.nodesh.GetSize());
    for (unsigned int i=0; i<s.nodesh.GetSize()-1; i++)
    {
        LMM(i)+=s.lineslen(i)/((nprec)2.0);
        LMM(i+1)+=s.lineslen(i)/((nprec)2.0);
    };

    nprec rhoq=par.bdens*2.0;
    nprec Cd=3.0;
    nprec zm=0.04; //meters
    nprec z0=2.5e-5;//par.z00;
    nprec z1=0.005; //meters
    nprec gamma=par.SLgamma;
    nprec alpha=par.SLalpha;
    nprec g=9.81;
    nprec Phia=5.7e-4;
    nprec k=0.41;
    //nprec err=std::numeric_limits<nprec>::max();
    //nprec erra=std::numeric_limits<nprec>::max();
    {
    nprec eps=1e-6;
    unsigned int n1,n2;
    nprec len,ue,flr,flu;
    matrix<nprec> F(s.nodesh.GetSize(),2); //Source Terms
    Vector<nprec> U(s.nodesh.GetSize());
    nprec tau,taug0,tauta,taut,A1,A0,aux,aux2,aux3,veffm,SF1,SF2;
    Vector<nprec> uastt(s.nodesh.GetSize());
    Vector<nprec> gradh(s.nodesh.GetSize());
    uastt.Initialize((nprec)0.0);
    gradh.Initialize((nprec)0.0);
    {
        Vector<unsigned int> cont(s.nodesh.GetSize());
        cont.Initialize(0);
        nprec auxp,auxp2;
        for (unsigned int i=0; i<s.nodesh.GetSize()-1; i++)
        {
            len = s.lineslen.Read(i);
            SF1 = -(nprec)1.0/len;
            SF2 = -SF1;
            auxp=s.h(i)*SF1+s.h(i+1)*SF2;
            gradh(i)+=auxp;
            gradh(i+1)+=auxp;
            auxp2=par.veltresh;//*std::sqrt(std::cos(std::atan(auxp))+std::sin(std::atan(auxp))/std::sin(par.angfri*M_PI/180.0));
            uastt(i)+=auxp2;
            uastt(i+1)+=auxp2;
            cont(i)++;
            cont(i+1)++;
        };
        for (unsigned int i=0; i<s.nodesh.GetSize(); i++)
        {
            gradh(i)/=(nprec)cont(i);
            uastt(i)/=(nprec)cont(i);
        };
    }
    unsigned int niter=0;
    nprec timef=par.DELTP*par.tfactor;
    nprec totalt=(nprec)0.0;
    nprec newdt;
    nprec dt=par.DELTP*par.SLdtfactor;
    while (totalt<timef)
    {
        niter++;
        LOsol.Initialize((nprec)0.0);
        F.Initialize((nprec)0.0);
        U.Initialize((nprec)0.0);
        newdt=dt;

        for (unsigned int i=0; i<s.nodesh.GetSize(); i++)
        {
            if(s.nodeshsedrho.Read(i)>eps)
            {
                U(i)=s.nodeshqs.Read(i)/s.nodeshsedrho.Read(i);
            }
            else
            {
                U(i)=(nprec)0.0;
            };
            aux=std::abs(s.nodeshuast.Read(i));
            aux2=std::abs(U(i));
            tau=aux*aux*par.rho0;
            taut=uastt(i)*uastt(i)*par.rho0;
            tauta=taut/0.8;
            taug0=s.nodeshsedrho.Read(i)*g/(2.0*alpha);
            if (aux2>eps)
            {
                F(i,0)=gamma*g/(2.0*alpha)*(tau-taut)/taut*s.nodeshsedrho.Read(i)/aux2;
                F(i,0)-=gamma*std::pow(g/(2.0*alpha),2.0)*s.nodeshsedrho.Read(i)*s.nodeshsedrho.Read(i)/(aux2*taut);
            }
            else
            {
                F(i,0)=0.0;
            };
            if (s.h(i)>1e-3){
            F(i,0)+=Phia*std::max((tau-taug0)/tauta-1.0,0.0);};

            veffm=0.0;
            if (tau>taug0)
            {
                A1=std::sqrt(1.0+z1/zm*taug0/(tau-taug0));
                A0=std::sqrt(1.0+z0/zm*taug0/(tau-taug0));
                if (taug0>eps)
                {
                    veffm=aux/k*std::sqrt(1.0-taug0/tau)*(2.0*A1-2.0*A0+std::log((A1-1.0)*(A0+1.0)/((A1+1.0)*(A0-1.0))));
                }
                else
                {
                    veffm=aux/k*std::log(z1/z0);
                };
                //veffm=aux/k*std::sqrt(1.0-taug0/tau)*(2.0*A1-2.0+std::log(z1/z0));
                //veffm=std::min(veffm,veffmax);
            }
            else
            {
                veffm=2.0*aux/k*(std::sqrt(z1/zm)-std::sqrt(z0/zm));
            };

            veff(i)=0.0;
            if (aux>eps)
            {
                veff(i)=s.nodeshuast.Read(i)/aux*veffm;
            };
            s.nodeshveff(i)=veff(i);
            aux3=std::abs(veff(i)-U(i));

            F(i,1)=3.0/4.0*Cd*par.rho0/rhoq*1.0/par.graind*(veff(i)-U(i))*aux3;
            if (aux>eps){F(i,1)-=g/(2.0*alpha)*s.nodeshuast.Read(i)/aux;};
            if (aux>uastt(i)){F(i,1)-=g*std::sin(std::atan(gradh(i)));};
            if(s.nodeshsedrho.Read(i)>eps)
            {
                F(i,1)-=U(i)*F(i,0)/s.nodeshsedrho.Read(i);
            };
        };

        for (unsigned int i=0; i<s.nodesh.GetSize()-1; i++)
        {
            n1=i; n2=i+1;
            len = s.lineslen.Read(i);
            ue=(U(n1)+U(n2))/((nprec)2.0);

            //Low Order Solution
            flr=(nprec)0.0;
            flu=(nprec)0.0;
            if (ue>(nprec)0.0)
            {
                flr=ue*(s.nodeshsedrho.Read(n1));
                flu=ue*(U.Read(n1));
            }else
            {
                flr=ue*(s.nodeshsedrho.Read(n2));
                flu=ue*(U.Read(n2));
            };
            LOsol(n1,0)-=flr; LOsol(n1,1)-=flu;
            LOsol(n2,0)+=flr; LOsol(n2,1)+=flu;

            //Low Order Source Term
            LOsol(n1,0)+=F(n1,0)*len/2.0;
            LOsol(n1,1)+=F(n1,1)*len/2.0;
            LOsol(n2,0)+=F(n2,0)*len/2.0;
            LOsol(n2,1)+=F(n2,1)*len/2.0;

            //Boundary Integrals
            if (i==0)
            {
                //LOsol(i,0)+=U(i)*s.nodeshsedrho.Read(i);
                //LOsol(i,1)+=U(i)*U.Read(i);
            };
            if (i==(s.nodesh.GetSize()-2))
            {
                LOsol(i+1,0)-=U(i+1)*s.nodeshsedrho(i+1);
                LOsol(i+1,1)-=U(i+1)*U.Read(i+1);
            };
        };

        for (unsigned int i=0; i<s.nodesh.GetSize(); i++)
        {
            if (LOsol(i,0)*dt/LMM(i)+s.nodeshsedrho(i)<(nprec)0.0)
            {
                newdt=std::min(newdt,-s.nodeshsedrho(i)/LOsol(i,0)*LMM(i));
            };
        };
        if (totalt+newdt>=timef){newdt=timef-totalt+1e-15;};
        //err=(nprec)0.0;
        totalt+=newdt;
        for (unsigned int i=0; i<s.nodesh.GetSize(); i++)
        {
            Force(i)+=newdt*3.0/4.0*Cd*par.rho0/rhoq*1.0/par.graind*(veff(i)-U(i))*std::abs(veff(i)-U(i))*s.nodeshsedrho.Read(i);
            LOsol(i,0)*=newdt/LMM(i);
            LOsol(i,1)*=newdt/LMM(i);
            //err+=LOsol(i,0)*LOsol(i,0)+LOsol(i,1)*LOsol(i,1);
            s.nodeshsedrho(i)+=LOsol(i,0);
            if (s.nodeshsedrho(i)<0.0){s.nodeshsedrho(i)=0.0;}; //To avoid very small negative values. Negative values are prevented with adaptive timestep
            s.nodeshMeanqs(i)+=(s.nodeshqs(i)+(LOsol(i,1)+U(i))*s.nodeshsedrho(i))/((nprec)2.0)*dt;
            s.nodeshqs(i)=(LOsol(i,1)+U(i))*s.nodeshsedrho(i);
        };
        //err=std::sqrt(err/((nprec)s.nodesh.GetSize()*2.0));
        s.nodeshqs(0)=0.0;
        s.nodeshsedrho(0)=0.0;
        //if(err>erra){err=0.0;};
        //erra=err;
    };

    nprec meanF=0.0;
    for (unsigned int i=0; i<s.nodesh.GetSize(); i++)
    {
        Force(i)/=totalt;
	meanF+=Force(i);
        s.nodeshMeanqs(i)/=totalt;
    };
    meanF/=(nprec)s.nodesh.GetSize();
    s.nodeshMeanqs(0)=0.0;

    std::cout << "Number of iterations: " << niter << ", Mean Force: " << meanF << std::endl;
    }

    //Pass From Bed Nodes to Mesh nodes
    av.TwoWayForce.Initialize((nprec)0.0);
    {
    Vector<bool> Check(geo.nd.GetNnodes());
    Check.Initialize(false);
    unsigned int nb1,nb2;
    nprec X,x1,x2;
    bool found;
    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        if (geo.bd.ReadBC(i)==20)
        {
            nb1=geo.bd.ReadNode(i,0); nb2=geo.bd.ReadNode(i,1);

            if (!Check(nb1))
            {
                Check(nb1)=true;
                X=geo.nd.ReadC(nb1,0);
                found=false;
                for (unsigned int j=0; j<s.linesh.GetRows(); j++)
                {
                    x1=s.nodeshX.Read(s.linesh.Read(j,0));
                    x2=s.nodeshX.Read(s.linesh.Read(j,1));
                    if (X>x1-eps&&X<x2+eps)
                    {
                        found=true;
                        av.TwoWayForce(nb1)=Force(s.linesh.Read(j,0))+((X-x1)/(x2-x1))*(Force(s.linesh.Read(j,1))-Force(s.linesh.Read(j,0)));
                        break;
                    };

                };
                if (!found)
                {
                    std::cout << "Correspondence for Node " << nb1 << " in fluid mesh not found" << std::endl;
                    exit(1);
                };
            };

            if (!Check(nb2))
            {
                Check(nb2)=true;
                X=geo.nd.ReadC(nb2,0);
                found=false;
                for (unsigned int j=0; j<s.linesh.GetRows(); j++)
                {
                    x1=s.nodeshX.Read(s.linesh.Read(j,0));
                    x2=s.nodeshX.Read(s.linesh.Read(j,1));
                    if (X>x1-eps&&X<x2+eps)
                    {
                        found=true;
                        av.TwoWayForce(nb2)=Force(s.linesh.Read(j,0))+((X-x1)/(x2-x1))*(Force(s.linesh.Read(j,1))-Force(s.linesh.Read(j,0)));
                        break;
                    };

                };
                if (!found)
                {
                    std::cout << "Correspondence for Node " << nb2 << " in fluid mesh not found" << std::endl;
                    exit(1);
                };
            };
        };
    };
    }
}


void HnodesCalculation(AuxVarGeo& s, const geometry& geo, const matrix<nprec>& unknowns1, parameter& par, const matrix<nprec>& SFD, nprec timepar, nprec currenttime)
{
    Vector<nprec> LOsol(s.h.GetSize());
    Vector<nprec> HOsol(s.h.GetSize());
    LOsol.Initialize((nprec)0.0);
    HOsol.Initialize((nprec)0.0);
    matrix<nprec> CMM(s.nodesh.GetSize(),3);    //pos 0 ---> current node, pos 1 ---> previous node, pos 2 ---> following node
    CMM.Initialize((nprec)0.0);
    Vector<nprec> LMM(s.nodesh.GetSize());
    LMM.Initialize((nprec)0.0);

    //Advection Equation
    {
    //Parameters for diffusion K according to P. Ortiz, 2009. Coupling the dynamics of boundary layers and evolutionary dunes.
    nprec beta=par.aval;     //beta
    nprec sc=std::tan(M_PI/180.0*par.angfri);       //critical slope
    Vector<nprec> rhs(s.nodesh.GetSize());
    rhs.Initialize((nprec)0.0);
    Vector<nprec> rhsup(s.nodesh.GetSize());
    rhsup.Initialize((nprec)0.0);
    Vector<nprec> U(s.nodesh.GetSize());
    nprec SFD1,SFD2,dhdx,K,Um,flux,dudx,hm;
    for (unsigned int i=0; i<s.nodesh.GetSize(); i++)
    {
        if (s.h.Read(i)>1e-2)   //No poner valores menores que este porque si no en el Upwind salen valores negativos (por temas de redondeo chungos)!!!!
        {
            U(i)=s.nodeshMeanqs.Read(i)*par.tfactor/(par.bdens*s.h.Read(i));
        }
        else
        {
            U(i)=(nprec)0.0;
        };
    };

    if (par.dyntimestep)
    {
    nprec um;
    nprec timestep=par.DELTP;
    nprec maxCo=(nprec)0.0;
    nprec minlen=s.lineslen(0);
    for (unsigned int i=0; i<s.nodesh.GetSize()-1; i++)
    {
        um=U(i)+U(i+1);
        um/=(nprec)2.0;
        um=std::abs(um);
        if (um>1e-15)
        {
            timestep=std::min(timestep,s.lineslen(i)/um*par.maxCouH);
        };
        maxCo=std::max(maxCo,timestep*um/s.lineslen(i));
        minlen=std::min(minlen,s.lineslen(i));
    };
    if (timepar+timestep>par.interplot){timestep=par.interplot-timepar;};
    if (currenttime+timestep>par.tfinal){timestep=par.tfinal-currenttime;};

    par.DELTP=timestep;

    beta=std::min(beta,(minlen*minlen)/((nprec)2.0*timestep)*0.25);
    std::cout << "Max Courant Number for H transport: " << maxCo << std::endl;
    std::cout << "Avalanche term: " << beta << std::endl;
    };

    for (unsigned int i=0; i<s.nodesh.GetSize()-1; i++)
    {
        SFD1=-(nprec)1.0/(s.lineslen(i));
        SFD2=-SFD1;

        CMM(i,0)+=s.lineslen(i)/((nprec)3.0);
        CMM(i,2)+=s.lineslen(i)/((nprec)6.0);
        CMM(i+1,0)+=s.lineslen(i)/((nprec)3.0);
        CMM(i+1,1)+=s.lineslen(i)/((nprec)6.0);
        LMM(i)+=s.lineslen(i)/((nprec)2.0);
        LMM(i+1)+=s.lineslen(i)/((nprec)2.0);

        dhdx=SFD1*s.h(i)+SFD2*s.h(i+1);
        dudx=SFD1*U(i)+SFD2*U(i+1);
        Um=(U(i)+U(i+1))/((nprec)2.0);
        hm=(s.h(i)+s.h(i+1))/((nprec)2.0);

        if ((std::abs(dhdx)-sc)>(nprec)0.0)
        {
            K=beta;
        }
        else
        {
            K=(nprec)0.0;
        };

        //First order terms
        rhs(i)+=(s.lineslen(i)/((nprec)2.0))*SFD1*(U(i)*s.h(i)+U(i+1)*s.h(i+1));
        rhs(i+1)+=(s.lineslen(i)/((nprec)2.0))*SFD2*(U(i)*s.h(i)+U(i+1)*s.h(i+1));

        rhs(i)-=s.lineslen(i)*K*SFD1*dhdx;
        rhs(i+1)-=s.lineslen(i)*K*SFD2*dhdx;

        //Upwind terms
        if (Um>(nprec)0.0)
        {
            flux=Um*s.h(i);
        }
        else
        {
            flux=Um*s.h(i+1);
        };
        rhsup(i)-=flux;
        rhsup(i+1)+=flux;

        rhsup(i)-=s.lineslen(i)*K*SFD1*dhdx;
        rhsup(i+1)-=s.lineslen(i)*K*SFD2*dhdx;

        //Boundary Integrals
        if (i==0)
        {
            rhs(i)+=U(i)*s.h(i);
            rhs(i)+=(-1.0)*K*dhdx;
            rhsup(i)+=U(i)*s.h(i);
            rhsup(i)+=(-1.0)*K*dhdx;

            //Second order terms
            rhs(i)+=par.DELTP/((nprec)2.0)*(-1.0)*s.h(i)*dudx*U(i);
            rhs(i)+=par.DELTP/((nprec)2.0)*(-1.0)*U(i)*U(i)*dhdx;
        };
        if (i==(s.nodesh.GetSize()-2))
        {
            rhs(i+1)-=U(i+1)*s.h(i+1);
            rhs(i+1)+=(1.0)*K*dhdx;
            rhsup(i+1)-=U(i+1)*s.h(i+1);
            rhsup(i+1)+=(1.0)*K*dhdx;

            //Second order terms
            rhs(i+1)+=par.DELTP/((nprec)2.0)*(1.0)*s.h(i+1)*dudx*U(i+1);
            rhs(i+1)+=par.DELTP/((nprec)2.0)*(1.0)*U(i+1)*U(i+1)*dhdx;
        };

        //Second order terms
        rhs(i)+=par.DELTP/((nprec)2.0)*(s.lineslen(i)/((nprec)6.0))*dudx*dhdx*((nprec)2.0*U(i)+U(i+1));
        rhs(i+1)+=par.DELTP/((nprec)2.0)*(s.lineslen(i)/((nprec)6.0))*dudx*dhdx*((nprec)2.0*U(i+1)+U(i));

        rhs(i)-=par.DELTP/((nprec)2.0)*s.lineslen(i)*SFD1*Um*hm*dudx;
        rhs(i+1)-=par.DELTP/((nprec)2.0)*s.lineslen(i)*SFD2*Um*hm*dudx;

        rhs(i)-=par.DELTP/((nprec)2.0)*s.lineslen(i)*Um*Um*dhdx*SFD1;
        rhs(i+1)-=par.DELTP/((nprec)2.0)*s.lineslen(i)*Um*Um*dhdx*SFD2;
    };

    for (unsigned int i=0; i<s.h.GetSize(); i++)
    {
        rhsup(i)*=par.DELTP/LMM(i);
        LOsol(i)=s.h(i)+rhsup(i);
        rhs(i)*=par.DELTP;

        if (i!=0)
        {
            rhs(i)=(rhs(i)-rhs(i-1)*CMM(i,1))/(CMM(i,0)-CMM(i-1,2)*CMM(i,1));
            if (i<s.h.GetSize()-1)
            {
                CMM(i,2)/=(CMM(i,0)-CMM(i-1,2)*CMM(i,1));
            };
        }
        else
        {
            rhs(i)/=CMM(i,0);
            CMM(i,2)/=CMM(i,0);
        };
    };

    for (int i=(int)s.h.GetSize()-1; i>=0; i--)
    {
        if (i!=((int)s.h.GetSize()-1))
        {
            rhs((unsigned int)i)=rhs((unsigned int)i)-CMM((unsigned int)i,2)*rhs((unsigned int)i+1);
        };

        HOsol((unsigned int)i)=s.h((unsigned int)i)+rhs((unsigned int)i);
    };
    }

    //FCT
    {
    Vector<nprec> BMAX(s.h.GetSize());
    Vector<nprec> BMIN(s.h.GetSize());
    BMAX.Initialize(std::numeric_limits<nprec>::min());
    BMIN.Initialize(std::numeric_limits<nprec>::max());
    for (unsigned int i=0; i<s.h.GetSize()-1; i++)
    {
        BMAX(i)=std::max(std::max(std::max(BMAX(i),s.h.Read(i)),std::max(s.h.Read(i+1),LOsol.Read(i))),LOsol.Read(i+1));
        BMIN(i)=std::min(std::min(std::min(BMIN(i),s.h.Read(i)),std::min(s.h.Read(i+1),LOsol.Read(i))),LOsol.Read(i+1));
        BMAX(i+1)=std::max(std::max(std::max(BMAX(i+1),s.h.Read(i)),std::max(s.h.Read(i+1),LOsol.Read(i))),LOsol.Read(i+1));
        BMIN(i+1)=std::min(std::min(std::min(BMIN(i+1),s.h.Read(i)),std::min(s.h.Read(i+1),LOsol.Read(i))),LOsol.Read(i+1));
    };

    Vector<nprec> ASUM(s.h.GetSize());
    ASUM.Initialize((nprec)0.0);
    Vector<nprec> AREST(s.h.GetSize());
    AREST.Initialize((nprec)0.0);
    matrix<nprec> AC(s.h.GetSize()-1,4);
    Vector<nprec> ACf(s.h.GetSize());
    AC.Initialize((nprec)0.0);
    ACf.Initialize((nprec)0.0);
    Vector<nprec> ULE(s.h.GetSize());
    ULE.Initialize((nprec)0.0);
    Vector<nprec> DBCU(2);
    DBCU.Initialize((nprec)0.0);
    Vector<nprec> AEL(2);
    Vector<nprec> DBE(2);
    DBE.Initialize((nprec)0.0);

    for (unsigned int i=0; i<s.h.GetSize()-1; i++)
    {
        AEL.Initialize((nprec)0.0);

        DBE(0)=HOsol(i)-LOsol(i);
        DBE(1)=HOsol(i+1)-LOsol(i+1);

        DBCU(0)=s.lineslen(i)/((nprec)2.0)*DBE(0);
        DBCU(1)=s.lineslen(i)/((nprec)2.0)*DBE(1);

        AEL(0)+=DBCU(0);
        AEL(0)/=LMM.Read(i);
        AEL(1)+=DBCU(1);
        AEL(1)/=LMM.Read(i+1);

        if (AEL(0)>(nprec)0.0)
        {
            ASUM(i)+=AEL(0);
        }else
        {
            AREST(i)-=AEL(0);
        };
        if (AEL(1)>(nprec)0.0)
        {
            ASUM(i+1)+=AEL(1);
        }else
        {
            AREST(i+1)-=AEL(1);
        };
    };

    Vector<nprec> CP(s.h.GetSize());
    CP.Initialize((nprec)0.0);
    Vector<nprec> CM(s.h.GetSize());
    CM.Initialize((nprec)0.0);
    nprec betap, betan;

    for (unsigned int i=0; i<s.h.GetSize(); i++)
    {
        betap=(BMAX(i)-LOsol(i))/(ASUM(i)+1.0e-7);
        betan=(LOsol(i)-BMIN(i))/(AREST(i)+1.0e-7);

        CP(i)=std::min((nprec)1.0,betap);
        CM(i)=std::min((nprec)1.0,betan);
    };

    nprec cpe,cme;
    for (unsigned int i=0; i<s.h.GetSize()-1; i++)
    {
        AEL.Initialize((nprec)0.0);

        cpe=std::min(CP(i),CP(i+1));
        cme=std::min(CM(i),CM(i+1));

        DBE(0)=HOsol(i)-LOsol(i);
        DBE(1)=HOsol(i+1)-LOsol(i+1);

        DBCU(0)=s.lineslen(i)/((nprec)2.0)*DBE(0);
        DBCU(1)=s.lineslen(i)/((nprec)2.0)*DBE(1);

        AEL(0)+=DBCU(0);
        AEL(0)/=LMM.Read(i);
        AEL(1)+=DBCU(1);
        AEL(1)/=LMM.Read(i+1);

        if (AEL(0)>(nprec)0.0)
        {
            AC(i,0)=cpe*AEL(0);
        }else
        {
            AC(i,1)=cme*AEL(0);
        };
        if (AEL(1)>(nprec)0.0)
        {
            AC(i,2)=cpe*AEL(1);
        }else
        {
            AC(i,3)=cme*AEL(1);
        };
    };

    for (unsigned int i=0; i<s.h.GetSize()-1; i++)
    {
        ACf(i)+=AC(i,0)+AC(i,1);
        ACf(i+1)+=AC(i,2)+AC(i,3);
    };

    HOsol=LOsol+ACf;
    }

    for (unsigned int i=0; i<s.h.GetSize(); i++)
    {
        if (s.nodeshX(i)<=5.0||s.nodeshX(i)>=195.0)
        {
            HOsol(i)=s.h(i);
        };
    };

    s.h=HOsol;

    //Hmax update
    {
    s.Hmax=std::numeric_limits<nprec>::min();
    //nprec hmin=std::numeric_limits<nprec>::max();
    for (unsigned int i=0; i<s.h.GetSize(); i++)
    {
        s.Hmax=std::max(s.h(i),s.Hmax);
        //hmin=std::min(s.h(i),hmin);
        /*
        if (s.h(i)<0.0)
        {
            std::cout << "Negativo: Node " << i << ", position = " << s.nodeshX.Read(i) << ", hvalue = " << s.h(i) << std::endl;
        };
        */
    };
    //std::cout << "H min = " << hmin << std::endl;
    //std::cout << "H max = " << s.Hmax << std::endl;
    }
};


void NewMesh(const geometry& geo0, geometry& geo1, parameter& par, AuxVarGeo& avG, const matrix<nprec>& SFD0, matrix<nprec>& SFD1, UpwindVar& upvel0, UpwindVar& upvel1, matrix<nprec>& RHSUW1, matrix<nprec>& unknowns0, matrix<nprec>& unknowns1)
{
    unsigned int nnodes=0;
    unsigned int nelem=0;
    unsigned int nbound=0;
    nprec tol=(nprec)0.10;

    for (unsigned int i=0; i<geo0.nd.GetNnodes(); i++)
    {
        if (avG.ActivatedNodes.Read(i)!=0)
        {
            unknowns0(i,0)=unknowns1.Read((unsigned int)avG.OldNodesIndex.Read(i),0);
            unknowns0(i,1)=unknowns1.Read((unsigned int)avG.OldNodesIndex.Read(i),1);
            unknowns0(i,2)=unknowns1.Read((unsigned int)avG.OldNodesIndex.Read(i),2);
            unknowns0(i,3)=unknowns1.Read((unsigned int)avG.OldNodesIndex.Read(i),3);
        };
    };

    //Nodes Classification
    Vector<unsigned int> activatednodesprov(geo0.nd.GetNnodes());
    Vector<nprec> HMeshNodes(geo0.nd.GetNnodes());
    matrix<nprec> NodeIntersec(avG.Edges.GetRows(),2);
    avG.NodeNumberEdge.Initialize(-1);
    {
    HMeshNodes.Initialize((nprec)0.0);
    NodeIntersec.Initialize((nprec)0.0);
    activatednodesprov.Initialize(0);
    nprec X,Y,X0,X1,H0,H1,H;
    unsigned int inter;
    for (unsigned int i=0; i<geo0.nd.GetNnodes(); i++)
    {
        X=geo0.nd.ReadC(i,0);
        Y=geo0.nd.ReadC(i,1);
        inter=avG.nodeshX.GetSize();
        for (unsigned int j=0; j<avG.nodeshX.GetSize()-1; j++)
        {
            if ((avG.nodeshX(j)-eps<=X)&&(avG.nodeshX(j+1)+eps>=X))
            {
                inter=j;
                break;
            };
        };
        if (inter==avG.nodeshX.GetSize())
        {
            std::cout << "Error when searching node " << i << std::endl;
            exit(1);
        };
        X0=avG.nodeshX(inter);
        X1=avG.nodeshX(inter+1);
        H0=avG.h(inter)+geo0.nd.ReadC(avG.nodesh.Read(inter),1);
        H1=avG.h(inter+1)+geo0.nd.ReadC(avG.nodesh.Read(inter+1),1);
        H=H0+(H1-H0)/(X1-X0)*(X-X0);
        HMeshNodes(i)=H;
        if ((H-eps<Y)&&(H+eps>Y))
        {
            activatednodesprov(i)=2; //Coincident
            HMeshNodes(i)=Y;
        }
        else if (H<Y)
        {
            activatednodesprov(i)=1;
        }
        else
        {
            activatednodesprov(i)=0;
        };
    };
    }


    //Edges Classification
    {
    unsigned int n1,n2;
    nprec Y3,Y4,X3,X4,Y1,Y2,X1,X2,aux,PX,PY,LEN,L1,L2;
    bool check=true;
    while (check)
    {
    check=false;
    for (unsigned int i=0; i<avG.Edges.GetRows(); i++)
    {
        n1=avG.Edges.Read(i,0);
        n2=avG.Edges.Read(i,1);

        if (activatednodesprov.Read(n1)==1&&activatednodesprov.Read(n2)==1)
        {
            avG.EdgeCase(i)=1;
        }
        else if (activatednodesprov.Read(n1)==0&&activatednodesprov.Read(n2)==0)
        {
            avG.EdgeCase(i)=0;
        }
        else if (activatednodesprov.Read(n1)==2&&activatednodesprov.Read(n2)==2)
        {
            avG.EdgeCase(i)=3;
        }
        else if ((activatednodesprov.Read(n1)==0&&activatednodesprov.Read(n2)==1)||(activatednodesprov.Read(n1)==1&&activatednodesprov.Read(n2)==0))
        {
            avG.EdgeCase(i)=2;
            X1=geo0.nd.ReadC(n1,0);
            X2=geo0.nd.ReadC(n2,0);
            Y1=geo0.nd.ReadC(n1,1);
            Y2=geo0.nd.ReadC(n2,1);
            X3=X1;
            X4=X2;
            Y3=HMeshNodes.Read(n1);
            Y4=HMeshNodes.Read(n2);

            if ((X1-eps<X2)&&(X1+eps>X2))
            {
                PX=(X1+X2)/((nprec)2.0);
                PY=(Y3+Y4)/((nprec)2.0);

                LEN=std::sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2));
                L1=std::sqrt((PX-X1)*(PX-X1)+(PY-Y1)*(PY-Y1));
                L2=std::sqrt((PX-X2)*(PX-X2)+(PY-Y2)*(PY-Y2));

                if (L1<tol*LEN)
                {
                    check=true;
                    HMeshNodes(n1)=Y1;
                    if ((X1-eps<X2)&&(X1+eps>X2))
                    {
                        HMeshNodes(n2)=Y1;
                    };
                    activatednodesprov(n1)=2;
                    if (activatednodesprov(n2)==0)
                    {
                        avG.EdgeCase(i)=0;
                    }
                    else
                    {
                        avG.EdgeCase(i)=1;
                    };
                }
                else if (L2<tol*LEN)
                {
                    check=true;
                    HMeshNodes(n2)=Y2;
                    if ((X1-eps<X2)&&(X1+eps>X2))
                    {
                        HMeshNodes(n1)=Y2;
                    };
                    activatednodesprov(n2)=2;
                    if (activatednodesprov(n1)==0)
                    {
                        avG.EdgeCase(i)=0;
                    }
                    else
                    {
                        avG.EdgeCase(i)=1;
                    };
                }
                else
                {
                    NodeIntersec(i,0)=PX;
                    NodeIntersec(i,1)=PY;
                };
            }
            else
            {
                aux=(X1-X2)*(Y3-Y4)-(Y1-Y2)*(X3-X4);
                if (std::abs(aux)>eps)
                {
                    PX=(X1*Y2-Y1*X2)*(X3-X4)-(X1-X2)*(X3*Y4-Y3*X4);
                    PX/=aux;

                    PY=(X1*Y2-Y1*X2)*(Y3-Y4)-(Y1-Y2)*(X3*Y4-Y3*X4);
                    PY/=aux;

                    LEN=std::sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2));
                    L1=std::sqrt((PX-X1)*(PX-X1)+(PY-Y1)*(PY-Y1));
                    L2=std::sqrt((PX-X2)*(PX-X2)+(PY-Y2)*(PY-Y2));

                    if (L1<tol*LEN)
                    {
                        check=true;
                        HMeshNodes(n1)=Y1;
                        if ((X1-eps<X2)&&(X1+eps>X2))
                        {
                            HMeshNodes(n2)=Y1;
                        };
                        activatednodesprov(n1)=2;
                        if (activatednodesprov(n2)==0)
                        {
                            avG.EdgeCase(i)=0;
                        }
                        else
                        {
                            avG.EdgeCase(i)=1;
                        };
                    }
                    else if (L2<tol*LEN)
                    {
                        check=true;
                        HMeshNodes(n2)=Y2;
                        if ((X1-eps<X2)&&(X1+eps>X2))
                        {
                            HMeshNodes(n1)=Y2;
                        };
                        activatednodesprov(n2)=2;
                        if (activatednodesprov(n1)==0)
                        {
                            avG.EdgeCase(i)=0;
                        }
                        else
                        {
                            avG.EdgeCase(i)=1;
                        };
                    }
                    else
                    {
                        NodeIntersec(i,0)=PX;
                        NodeIntersec(i,1)=PY;
                    };
                }
                else
                {
                    check=true;
                    HMeshNodes(n1)=Y1;
                    activatednodesprov(n1)=2;

                    HMeshNodes(n2)=Y2;
                    activatednodesprov(n2)=2;
                };
            };
        }
        else
        {
            if (activatednodesprov.Read(n1)==0||activatednodesprov.Read(n2)==0)
            {
                avG.EdgeCase(i)=0;
            }
            else if ((activatednodesprov.Read(n1)==1||activatednodesprov.Read(n2)==1))
            {
                avG.EdgeCase(i)=1;
            }
            else
            {
                std::cout << "Error in Edges Classification" << std::endl;
                exit(1);
            };
        };
    };
    };
    }

    //Nodes counting
    {
    nnodes=0;
    for (unsigned int i=0; i<geo0.nd.GetNnodes(); i++)
    {
        if (activatednodesprov(i)!=0)
        {
            nnodes++;
        };
    };
    for (unsigned int i=0; i<avG.Edges.GetRows(); i++)
    {
        if (avG.EdgeCase.Read(i)==2)
        {
            nnodes++;
        };
    };
    }

    //Element Classification and Bound Counting (Part 1)
    {
    unsigned int e1,e2,e3;
    unsigned int n1,n2,n3;
    unsigned int aux,aux2;
    for (unsigned int i=0; i<geo0.elem.GetNelem(); i++)
    {
        e1=avG.ElemEdges.Read(i,0);
        e2=avG.ElemEdges.Read(i,1);
        e3=avG.ElemEdges.Read(i,2);
        n1=geo0.elem.ReadNode(i,0);
        n2=geo0.elem.ReadNode(i,1);
        n3=geo0.elem.ReadNode(i,2);

        if (avG.EdgeCase.Read(e1)==0&&avG.EdgeCase.Read(e2)==0&&avG.EdgeCase.Read(e3)==0)
        {
            avG.ActivatedElem(i)=0;
        }
        else if (avG.EdgeCase.Read(e1)==1&&avG.EdgeCase.Read(e2)==1&&avG.EdgeCase.Read(e3)==1)
        {
            avG.ActivatedElem(i)=1;

            nelem++;
        }
        else if (avG.EdgeCase.Read(e1)==3||avG.EdgeCase.Read(e2)==3||avG.EdgeCase.Read(e3)==3)
        {
            if (activatednodesprov.Read(n1)==1||activatednodesprov.Read(n2)==1||activatednodesprov.Read(n3)==1)
            {
                avG.ActivatedElem(i)=1;
                nelem++;
                nbound++;           //count those bound with bc=20
            }
            else
            {
                avG.ActivatedElem(i)=0;
            };
        }
        else if (avG.EdgeCase.Read(e1)==2||avG.EdgeCase.Read(e2)==2||avG.EdgeCase.Read(e3)==2)
        {
            nbound++;           //count those bounds with bc=20

            avG.ActivatedElem(i)=2;
            aux=0;
            if (avG.EdgeCase.Read(e1)==2)
            {
                aux++;
            };
            if (avG.EdgeCase.Read(e2)==2)
            {
                aux++;
            };
            if (avG.EdgeCase.Read(e3)==2)
            {
                aux++;
            };

            switch(aux)
            {
            case 1:
                nelem++;
                break;
            case 2:
                aux2=0;
                if (activatednodesprov.Read(n1)==1)
                {
                    aux2++;
                };
                if (activatednodesprov.Read(n2)==1)
                {
                    aux2++;
                };
                if (activatednodesprov.Read(n3)==1)
                {
                    aux2++;
                };

                if (aux2==1)
                {
                    nelem++;
                }
                else if (aux2==2)
                {
                    nelem+=2;
                };
                break;
            case 3:
                std::cout << "Error, element with 3 edges intersected" << std::endl;
                exit(1);
                break;
            };
        }
        else
        {
            std::cout << "Case Not Considered" << std::endl;
            exit(1);
        };
    };
    }

    //Boundary counting (Part 2)
    {
    unsigned int ed;
    for (unsigned int i=0; i<geo0.bd.GetNbounds(); i++)
    {
        ed=avG.BoundEdges.Read(i);
        if (avG.EdgeCase(ed)==1||avG.EdgeCase(ed)==2)       //We only account activated boundaries and intersected boundaries
        {
            nbound++;
        };
    };
    }

    avG.OldNodesIndex.Initialize(-1);
    avG.OldElemIndex.Initialize(-1);
    geo1.nd.SetNnodes(nnodes);
    unknowns1.SetDimensions(nnodes,4);
    unknowns1.Initialize((nprec)0.0);
    geo1.elem.SetNelem(nelem);
    geo1.bd.SetNbounds(nbound);
    SFD1.SetDimensions(nelem,10);
    avG.NeigCorresp.SetDimensions(nelem);

    //Nodes Definition
    {
    nnodes=0;
    for (unsigned int i=0; i<geo0.nd.GetNnodes(); i++)
    {
        if (activatednodesprov.Read(i)!=0)
        {
            geo1.nd.SetXY(nnodes,geo0.nd.ReadC(i,0),geo0.nd.ReadC(i,1));
            avG.OldNodesIndex(i)=(int)nnodes;
            unknowns1(nnodes,0)=unknowns0.Read(i,0);
            unknowns1(nnodes,1)=unknowns0.Read(i,1);
            unknowns1(nnodes,2)=unknowns0.Read(i,2);
            unknowns1(nnodes,3)=unknowns0.Read(i,3);
            nnodes++;
        }
        else
        {
            unknowns0(i,1)=(nprec)0.0;
            unknowns0(i,2)=(nprec)0.0;
            unknowns0(i,3)=(nprec)0.0;
        };
    };
    for (unsigned int i=0; i<avG.Edges.GetRows(); i++)
    {
        if (avG.EdgeCase.Read(i)==2)
        {

            geo1.nd.SetXY(nnodes,NodeIntersec(i,0),NodeIntersec(i,1));
            avG.NodeNumberEdge(i)=(int)nnodes;
            if (avG.PrevEdgecutted(i))
            {
                unknowns1(nnodes,0)=avG.PrevEdgeValue(i);
            }
            else
            {
                nprec Z;
                if (avG.ActivatedNodes.Read(avG.Edges.Read(i,0))!=0)
                {
                    Z=geo0.nd.ReadC(avG.Edges.Read(i,0),1)-NodeIntersec(i,1);
                    unknowns1(nnodes,0)=unknowns0(avG.Edges.Read(i,0),0)-(par.rho0/par.rhoref+unknowns0(avG.Edges.Read(i,0),3)*(par.rho1-par.rho0)/par.rhoref)*(par.g*par.Lref/(par.Uref*par.Uref))*(Z);
                    unknowns1(nnodes,0)+=(unknowns0(avG.Edges.Read(i,0),1)*unknowns0(avG.Edges.Read(i,0),1)+unknowns0(avG.Edges.Read(i,0),2)*unknowns0(avG.Edges.Read(i,0),2))/((nprec)2.0)*(par.rho0/par.rhoref+unknowns0(avG.Edges.Read(i,0),3)*(par.rho1-par.rho0)/par.rhoref);
                }
                else
                {
                    Z=geo0.nd.ReadC(avG.Edges.Read(i,1),1)-NodeIntersec(i,1);
                    unknowns1(nnodes,0)=unknowns0(avG.Edges.Read(i,1),0)-(par.rho0/par.rhoref+unknowns0(avG.Edges.Read(i,1),3)*(par.rho1-par.rho0)/par.rhoref)*(par.g*par.Lref/(par.Uref*par.Uref))*(Z);
                    unknowns1(nnodes,0)+=(unknowns0(avG.Edges.Read(i,1),1)*unknowns0(avG.Edges.Read(i,1),1)+unknowns0(avG.Edges.Read(i,1),2)*unknowns0(avG.Edges.Read(i,1),2))/((nprec)2.0)*(par.rho0/par.rhoref+unknowns0(avG.Edges.Read(i,1),3)*(par.rho1-par.rho0)/par.rhoref);
                };
            };

            unknowns1(nnodes,1)=(nprec)0.0;
            unknowns1(nnodes,2)=(nprec)0.0;
            unknowns1(nnodes,3)=(nprec)0.0;
            nnodes++;
        };
    };
    }

    avG.ActivatedNodes=activatednodesprov;

    //Elem and bound definition
    {
    nelem=0;
    nbound=0;
    unsigned int n1,n2,n3,aux,aux2,e1,e2,e3,nn1,nn2,nn3,nn4;
    nprec dia1,dia2;
    Vector<nprec> X(4);
    Vector<nprec> Y(4);
    nprec twoarea;
    for (unsigned int i=0; i<geo0.elem.GetNelem(); i++)
    {
        n1=geo0.elem.ReadNode(i,0);
        n2=geo0.elem.ReadNode(i,1);
        n3=geo0.elem.ReadNode(i,2);
        e1=avG.ElemEdges.Read(i,0);
        e2=avG.ElemEdges.Read(i,1);
        e3=avG.ElemEdges.Read(i,2);

        if (avG.ActivatedElem(i)==1)
        {
            geo1.elem.SetNode(nelem,0,(unsigned int)avG.OldNodesIndex.Read(n1));
            geo1.elem.SetNode(nelem,1,(unsigned int)avG.OldNodesIndex.Read(n2));
            geo1.elem.SetNode(nelem,2,(unsigned int)avG.OldNodesIndex.Read(n3));
            avG.OldElemIndex(i)=nelem;
            avG.NeigCorresp(nelem)=i;
            for (unsigned int j=0; j<10; j++)
            {
                SFD1(nelem,j)=SFD0.Read(i,j);
            };

            if (avG.EdgeCase.Read(e1)==3||avG.EdgeCase.Read(e2)==3||avG.EdgeCase.Read(e3)==3)
            {
                geo1.bd.SetBC(nbound,20);
                geo1.bd.SetEO(nbound,nelem);
                if (avG.EdgeCase.Read(e1)==3)
                {
                    geo1.bd.SetNode(nbound,0,(unsigned int)avG.OldNodesIndex.Read(n1));
                    geo1.bd.SetNode(nbound,1,(unsigned int)avG.OldNodesIndex.Read(n2));
                }
                else if (avG.EdgeCase.Read(e2)==3)
                {
                    geo1.bd.SetNode(nbound,0,(unsigned int)avG.OldNodesIndex.Read(n2));
                    geo1.bd.SetNode(nbound,1,(unsigned int)avG.OldNodesIndex.Read(n3));
                }
                else
                {
                    geo1.bd.SetNode(nbound,0,(unsigned int)avG.OldNodesIndex.Read(n3));
                    geo1.bd.SetNode(nbound,1,(unsigned int)avG.OldNodesIndex.Read(n1));
                };

                nbound++;
            };

            nelem++;
        }
        else if (avG.ActivatedElem(i)==2)
        {
            aux=0;
            if (avG.EdgeCase.Read(e1)==2)
            {
                aux++;
            };
            if (avG.EdgeCase.Read(e2)==2)
            {
                aux++;
            };
            if (avG.EdgeCase.Read(e3)==2)
            {
                aux++;
            };

            switch(aux)
            {
            case 1:
                if (avG.EdgeCase.Read(e1)==2)
                {
                    if (activatednodesprov.Read(n1)==1)
                    {
                        geo1.elem.SetNode(nelem,0,(unsigned int)avG.OldNodesIndex.Read(n1));
                        geo1.elem.SetNode(nelem,1,(unsigned int)avG.NodeNumberEdge(e1));
                        geo1.elem.SetNode(nelem,2,(unsigned int)avG.OldNodesIndex.Read(n3));
                    }
                    else if (activatednodesprov.Read(n2)==1)
                    {
                        geo1.elem.SetNode(nelem,0,(unsigned int)avG.OldNodesIndex.Read(n2));
                        geo1.elem.SetNode(nelem,1,(unsigned int)avG.OldNodesIndex.Read(n3));
                        geo1.elem.SetNode(nelem,2,(unsigned int)avG.NodeNumberEdge(e1));
                    }
                    else
                    {
                        std::cout << "Strange case" << std::endl;
                        exit(1);
                    };
                }
                else if (avG.EdgeCase.Read(e2)==2)
                {
                    if (activatednodesprov.Read(n2)==1)
                    {
                        geo1.elem.SetNode(nelem,0,(unsigned int)avG.OldNodesIndex.Read(n2));
                        geo1.elem.SetNode(nelem,1,(unsigned int)avG.NodeNumberEdge(e2));
                        geo1.elem.SetNode(nelem,2,(unsigned int)avG.OldNodesIndex.Read(n1));
                    }
                    else if (activatednodesprov.Read(n3)==1)
                    {
                        geo1.elem.SetNode(nelem,0,(unsigned int)avG.OldNodesIndex.Read(n3));
                        geo1.elem.SetNode(nelem,1,(unsigned int)avG.OldNodesIndex.Read(n1));
                        geo1.elem.SetNode(nelem,2,(unsigned int)avG.NodeNumberEdge(e2));
                    }
                    else
                    {
                        std::cout << "Strange case" << std::endl;
                        exit(1);
                    };
                }
                else
                {
                    if (activatednodesprov.Read(n1)==1)
                    {
                        geo1.elem.SetNode(nelem,0,(unsigned int)avG.OldNodesIndex.Read(n1));
                        geo1.elem.SetNode(nelem,1,(unsigned int)avG.OldNodesIndex.Read(n2));
                        geo1.elem.SetNode(nelem,2,(unsigned int)avG.NodeNumberEdge(e3));
                    }
                    else if (activatednodesprov.Read(n3)==1)
                    {
                        geo1.elem.SetNode(nelem,0,(unsigned int)avG.OldNodesIndex.Read(n3));
                        geo1.elem.SetNode(nelem,1,(unsigned int)avG.NodeNumberEdge(e3));
                        geo1.elem.SetNode(nelem,2,(unsigned int)avG.OldNodesIndex.Read(n2));
                    }
                    else
                    {
                        std::cout << "Strange case" << std::endl;
                        exit(1);
                    };
                };

                X(0)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,0),0);
                X(1)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,1),0);
                X(2)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,2),0);
                Y(0)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,0),1);
                Y(1)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,1),1);
                Y(2)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,2),1);
                twoarea=(X(1)-X(0))*(Y(2)-Y(0))-(X(2)-X(0))*(Y(1)-Y(0));
                SFD1(nelem,0)=(Y(1)-Y(2))/twoarea;
                SFD1(nelem,1)=(Y(2)-Y(0))/twoarea;
                SFD1(nelem,2)=(Y(0)-Y(1))/twoarea;
                SFD1(nelem,3)=(X(2)-X(1))/twoarea;
                SFD1(nelem,4)=(X(0)-X(2))/twoarea;
                SFD1(nelem,5)=(X(1)-X(0))/twoarea;
                SFD1(nelem,6)=twoarea;
                SFD1(nelem,7)=(X(1)*Y(2)-X(2)*Y(1))/twoarea;
                SFD1(nelem,8)=(X(2)*Y(0)-X(0)*Y(2))/twoarea;
                SFD1(nelem,9)=(X(0)*Y(1)-X(1)*Y(0))/twoarea;

                geo1.bd.SetBC(nbound,20);
                geo1.bd.SetEO(nbound,nelem);
                geo1.bd.SetNode(nbound,0,geo1.elem.ReadNode(nelem,1));
                geo1.bd.SetNode(nbound,1,geo1.elem.ReadNode(nelem,2));
                avG.NeigCorresp(nelem)=i;
                nbound++;
                nelem++;
                break;
            case 2:
                aux2=0;
                if (activatednodesprov.Read(n1)==1)
                {
                    aux2++;
                };
                if (activatednodesprov.Read(n2)==1)
                {
                    aux2++;
                };
                if (activatednodesprov.Read(n3)==1)
                {
                    aux2++;
                };

                if (aux2==1)
                {
                    if (activatednodesprov.Read(n1)==1)
                    {
                        geo1.elem.SetNode(nelem,0,(unsigned int)avG.OldNodesIndex.Read(n1));
                        geo1.elem.SetNode(nelem,1,(unsigned int)avG.NodeNumberEdge(e1));
                        geo1.elem.SetNode(nelem,2,(unsigned int)avG.NodeNumberEdge(e3));
                    }
                    else if (activatednodesprov.Read(n2)==1)
                    {
                        geo1.elem.SetNode(nelem,0,(unsigned int)avG.OldNodesIndex.Read(n2));
                        geo1.elem.SetNode(nelem,1,(unsigned int)avG.NodeNumberEdge(e2));
                        geo1.elem.SetNode(nelem,2,(unsigned int)avG.NodeNumberEdge(e1));
                    }
                    else
                    {
                        geo1.elem.SetNode(nelem,0,(unsigned int)avG.OldNodesIndex.Read(n3));
                        geo1.elem.SetNode(nelem,1,(unsigned int)avG.NodeNumberEdge(e3));
                        geo1.elem.SetNode(nelem,2,(unsigned int)avG.NodeNumberEdge(e2));
                    };

                    X(0)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,0),0);
                    X(1)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,1),0);
                    X(2)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,2),0);
                    Y(0)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,0),1);
                    Y(1)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,1),1);
                    Y(2)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem,2),1);
                    twoarea=(X(1)-X(0))*(Y(2)-Y(0))-(X(2)-X(0))*(Y(1)-Y(0));
                    SFD1(nelem,0)=(Y(1)-Y(2))/twoarea;
                    SFD1(nelem,1)=(Y(2)-Y(0))/twoarea;
                    SFD1(nelem,2)=(Y(0)-Y(1))/twoarea;
                    SFD1(nelem,3)=(X(2)-X(1))/twoarea;
                    SFD1(nelem,4)=(X(0)-X(2))/twoarea;
                    SFD1(nelem,5)=(X(1)-X(0))/twoarea;
                    SFD1(nelem,6)=twoarea;
                    SFD1(nelem,7)=(X(1)*Y(2)-X(2)*Y(1))/twoarea;
                    SFD1(nelem,8)=(X(2)*Y(0)-X(0)*Y(2))/twoarea;
                    SFD1(nelem,9)=(X(0)*Y(1)-X(1)*Y(0))/twoarea;

                    geo1.bd.SetBC(nbound,20);
                    geo1.bd.SetEO(nbound,nelem);
                    geo1.bd.SetNode(nbound,0,geo1.elem.ReadNode(nelem,1));
                    geo1.bd.SetNode(nbound,1,geo1.elem.ReadNode(nelem,2));
                    avG.NeigCorresp(nelem)=i;
                    nbound++;
                    nelem++;
                }
                else if (aux2==2)
                {
                    if (activatednodesprov.Read(n1)==0)
                    {
                        nn1=(unsigned int)avG.NodeNumberEdge(e1);
                        nn2=(unsigned int)avG.OldNodesIndex.Read(n2);
                        nn3=(unsigned int)avG.OldNodesIndex.Read(n3);
                        nn4=(unsigned int)avG.NodeNumberEdge(e3);
                    }
                    else if (activatednodesprov.Read(n2)==0)
                    {
                        nn1=(unsigned int)avG.NodeNumberEdge(e2);
                        nn2=(unsigned int)avG.OldNodesIndex.Read(n3);
                        nn3=(unsigned int)avG.OldNodesIndex.Read(n1);
                        nn4=(unsigned int)avG.NodeNumberEdge(e1);
                    }
                    else
                    {
                        nn1=(unsigned int)avG.NodeNumberEdge(e3);
                        nn2=(unsigned int)avG.OldNodesIndex.Read(n1);
                        nn3=(unsigned int)avG.OldNodesIndex.Read(n2);
                        nn4=(unsigned int)avG.NodeNumberEdge(e2);
                    };

                    X(0)=geo1.nd.ReadC(nn1,0);
                    X(1)=geo1.nd.ReadC(nn2,0);
                    X(2)=geo1.nd.ReadC(nn3,0);
                    X(3)=geo1.nd.ReadC(nn4,0);
                    Y(0)=geo1.nd.ReadC(nn1,1);
                    Y(1)=geo1.nd.ReadC(nn2,1);
                    Y(2)=geo1.nd.ReadC(nn3,1);
                    Y(3)=geo1.nd.ReadC(nn4,1);
                    dia1=std::sqrt((X(0)-X(2))*(X(0)-X(2))+(Y(0)-Y(2))*(Y(0)-Y(2)));
                    dia2=std::sqrt((X(1)-X(3))*(X(1)-X(3))+(Y(1)-Y(3))*(Y(1)-Y(3)));

                    if (dia1<=dia2)
                    {
                        geo1.elem.SetNode(nelem,0,nn1);
                        geo1.elem.SetNode(nelem,1,nn2);
                        geo1.elem.SetNode(nelem,2,nn3);

                        geo1.elem.SetNode(nelem+1,0,nn1);
                        geo1.elem.SetNode(nelem+1,1,nn3);
                        geo1.elem.SetNode(nelem+1,2,nn4);

                        geo1.bd.SetBC(nbound,20);
                        geo1.bd.SetEO(nbound,nelem+1);
                        geo1.bd.SetNode(nbound,0,nn4);
                        geo1.bd.SetNode(nbound,1,nn1);
                        nbound++;
                    }
                    else
                    {
                        geo1.elem.SetNode(nelem,0,nn1);
                        geo1.elem.SetNode(nelem,1,nn2);
                        geo1.elem.SetNode(nelem,2,nn4);

                        geo1.elem.SetNode(nelem+1,0,nn2);
                        geo1.elem.SetNode(nelem+1,1,nn3);
                        geo1.elem.SetNode(nelem+1,2,nn4);

                        geo1.bd.SetBC(nbound,20);
                        geo1.bd.SetEO(nbound,nelem);
                        geo1.bd.SetNode(nbound,0,nn4);
                        geo1.bd.SetNode(nbound,1,nn1);
                        nbound++;
                    };
                    for (unsigned int j=0; j<2; j++)
                    {
                        X(0)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem+j,0),0);
                        X(1)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem+j,1),0);
                        X(2)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem+j,2),0);
                        Y(0)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem+j,0),1);
                        Y(1)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem+j,1),1);
                        Y(2)=geo1.nd.ReadC(geo1.elem.ReadNode(nelem+j,2),1);
                        twoarea=(X(1)-X(0))*(Y(2)-Y(0))-(X(2)-X(0))*(Y(1)-Y(0));
                        SFD1(nelem+j,0)=(Y(1)-Y(2))/twoarea;
                        SFD1(nelem+j,1)=(Y(2)-Y(0))/twoarea;
                        SFD1(nelem+j,2)=(Y(0)-Y(1))/twoarea;
                        SFD1(nelem+j,3)=(X(2)-X(1))/twoarea;
                        SFD1(nelem+j,4)=(X(0)-X(2))/twoarea;
                        SFD1(nelem+j,5)=(X(1)-X(0))/twoarea;
                        SFD1(nelem+j,6)=twoarea;
                        SFD1(nelem+j,7)=(X(1)*Y(2)-X(2)*Y(1))/twoarea;
                        SFD1(nelem+j,8)=(X(2)*Y(0)-X(0)*Y(2))/twoarea;
                        SFD1(nelem+j,9)=(X(0)*Y(1)-X(1)*Y(0))/twoarea;
                        avG.NeigCorresp(nelem+j)=i;
                    };

                    nelem+=2;
                };
                break;
            };
        };
    };
    }

    //Bound Definition (Part 2)
    {
    unsigned int ed,elold,n0,n1,auxn;
    bool found;
    for (unsigned int i=0; i<geo0.bd.GetNbounds(); i++)
    {
        ed=avG.BoundEdges.Read(i);
        elold=geo0.bd.ReadEO(i);
        if (avG.EdgeCase(ed)==1)       //We only account activated boundaries and intersected boundaries
        {
            if (avG.ActivatedElem(elold)==1)
            {
                geo1.bd.SetBC(nbound,geo0.bd.ReadBC(i));
                geo1.bd.SetEO(nbound,avG.OldElemIndex.Read(elold));
                geo1.bd.SetNode(nbound,0,avG.OldNodesIndex.Read(geo0.bd.ReadNode(i,0)));
                geo1.bd.SetNode(nbound,1,avG.OldNodesIndex.Read(geo0.bd.ReadNode(i,1)));
                nbound++;
            }
            else if (avG.ActivatedElem(elold)==2)
            {
                geo1.bd.SetBC(nbound,geo0.bd.ReadBC(i));
                geo1.bd.SetNode(nbound,0,avG.OldNodesIndex.Read(geo0.bd.ReadNode(i,0)));
                geo1.bd.SetNode(nbound,1,avG.OldNodesIndex.Read(geo0.bd.ReadNode(i,1)));
                found=false;
                n0=avG.OldNodesIndex.Read(geo0.bd.ReadNode(i,0));
                n1=avG.OldNodesIndex.Read(geo0.bd.ReadNode(i,1));
                for (unsigned int j=0; j<geo1.elem.GetNelem(); j++)
                {
                    if ((n0==geo1.elem.ReadNode(j,0)||n0==geo1.elem.ReadNode(j,1)||n0==geo1.elem.ReadNode(j,2))\
                        &&(n1==geo1.elem.ReadNode(j,0)||n1==geo1.elem.ReadNode(j,1)||n1==geo1.elem.ReadNode(j,2)))
                    {
                        geo1.bd.SetEO(nbound,j);
                        if (!((n0==geo1.elem.ReadNode(j,0)&&n1==geo1.elem.ReadNode(j,1))||\
                              (n0==geo1.elem.ReadNode(j,1)&&n1==geo1.elem.ReadNode(j,2))||\
                              (n0==geo1.elem.ReadNode(j,2)&&n1==geo1.elem.ReadNode(j,0))))
                        {
                            auxn=geo1.bd.ReadNode(nbound,0);
                            geo1.bd.SetNode(nbound,0,geo1.bd.ReadNode(nbound,1));
                            geo1.bd.SetNode(nbound,1,auxn);
                        };
                        found=true;
                        break;
                    };
                };
                if (!found)
                {
                    //std::cout << n0 << " " << n1 << std::endl;
                    std::cout << "Error: bound not found" << std::endl;
                    exit(1);
                };
                nbound++;
            };
        }
        else if (avG.EdgeCase(ed)==2)
        {
            if (activatednodesprov.Read(geo0.bd.ReadNode(i,0))==1)
            {
                n0=avG.OldNodesIndex.Read(geo0.bd.ReadNode(i,0));
                n1=(unsigned int)avG.NodeNumberEdge.Read(ed);
            }
            else if (activatednodesprov.Read(geo0.bd.ReadNode(i,1))==1)
            {
                n0=(unsigned int)avG.NodeNumberEdge.Read(ed);
                n1=avG.OldNodesIndex.Read(geo0.bd.ReadNode(i,1));
            }
            else
            {
                std::cout << "Error Found" << std::endl;
                exit(1);
            };
            geo1.bd.SetBC(nbound,geo0.bd.ReadBC(i));
            geo1.bd.SetNode(nbound,0,n0);
            geo1.bd.SetNode(nbound,1,n1);
            found=false;
            for (unsigned int j=0; j<geo1.elem.GetNelem(); j++)
            {
                if ((n0==geo1.elem.ReadNode(j,0)||n0==geo1.elem.ReadNode(j,1)||n0==geo1.elem.ReadNode(j,2))\
                    &&(n1==geo1.elem.ReadNode(j,0)||n1==geo1.elem.ReadNode(j,1)||n1==geo1.elem.ReadNode(j,2)))
                {
                    geo1.bd.SetEO(nbound,j);
                    if (!((n0==geo1.elem.ReadNode(j,0)&&n1==geo1.elem.ReadNode(j,1))||\
                          (n0==geo1.elem.ReadNode(j,1)&&n1==geo1.elem.ReadNode(j,2))||\
                          (n0==geo1.elem.ReadNode(j,2)&&n1==geo1.elem.ReadNode(j,0))))
                    {
                        auxn=geo1.bd.ReadNode(nbound,0);
                        geo1.bd.SetNode(nbound,0,geo1.bd.ReadNode(nbound,1));
                        geo1.bd.SetNode(nbound,1,auxn);
                    };
                    found=true;
                    break;
                };
            };
            if (!found)
            {
                //std::cout << n0 << " " << n1 << std::endl;
                std::cout << "Error: bound not found" << std::endl;
                exit(1);
            };
            nbound++;
        };
    };
    }

    {
    nprec dx,dy;
    for (unsigned int i=0; i<geo1.bd.GetNbounds(); i++)
    {
        dx = geo1.nd.ReadC(geo1.bd.ReadNode(i,1),0)-geo1.nd.ReadC(geo1.bd.ReadNode(i,0),0);
        dy = geo1.nd.ReadC(geo1.bd.ReadNode(i,1),1)-geo1.nd.ReadC(geo1.bd.ReadNode(i,0),1);
        geo1.bd.SetNormal(i,0,dy/(std::sqrt(dx*dx+dy*dy)));
        geo1.bd.SetNormal(i,1,-dx/(std::sqrt(dx*dx+dy*dy)));
        geo1.bd.SetLen(i,std::sqrt(dx*dx+dy*dy));
    };
    }

    if (par.indexfctu!=0)
    {
        Upwindgeo(geo1,upvel1);
        RHSUW1.SetDimensions(geo1.nd.GetNnodes(),3);//La columna 0 se queda libre por si se quiere meter ahi algo relativo a las presiones en el futuro
    };


    //Auxiliary vector for nodes with wall condition
    {
    unsigned int nwall=0;
    Vector<int> vhelp(geo1.nd.GetNnodes());
    vhelp.Initialize(-1);
    geo1.wnodes.SetDimensions(geo1.bd.GetNbounds()*2,3);
    geo1.wnodes.Initialize(-1);
    int nn,j;
    for (unsigned int i=0; i<geo1.bd.GetNbounds(); i++)
    {
        if (geo1.bd.ReadBC(i)==2||geo1.bd.ReadBC(i)==6||geo1.bd.ReadBC(i)==20)
        {
            nn=geo1.bd.ReadNode(i,0);
            j=vhelp(nn);
            if (j==-1)
            {
                nwall++;
                geo1.wnodes(nwall-1,0) = nn;
                geo1.wnodes(nwall-1,1) = i;
                vhelp(nn) = nwall-1;
            }
            else
            {
                geo1.wnodes(j,2) = i;
            };

            nn=geo1.bd.ReadNode(i,1);
            j=vhelp(nn);
            if (j==-1)
            {
                nwall++;
                geo1.wnodes(nwall-1,0) = nn;
                geo1.wnodes(nwall-1,1) = i;
                vhelp(nn) = nwall-1;
            }
            else
            {
                geo1.wnodes(j,2) = i;
            };
        };
    };

    int ib1,ib2;
    nprec ach;
    for (unsigned int i=0; i<nwall; i++)
    {
        ib1=geo1.wnodes.Read(i,1);
        ib2=geo1.wnodes.Read(i,2);
        if ((ib1!=-1)&&(ib2!=-1))
        {
            ach=geo1.bd.ReadNormal(ib1,0)*geo1.bd.ReadNormal(ib2,0)+geo1.bd.ReadNormal(ib1,1)*geo1.bd.ReadNormal(ib2,1);
            if (ach<(nprec)-0.2)
            {
                //std::cout << "Node " << geo1.wnodes.Read(i,0)+1 << " is excessive acute edge" << std::endl;
            };
        };
    };
    //Resize wnodes
    if (nwall==0)
    {
        geo1.wnodes.SetDimensions(0,0);
    }
    else
    {
        matrix<int> tmp(nwall,3);
        for (unsigned int i=0; i<nwall; i++)
        {
            tmp(i,0)=geo1.wnodes.Read(i,0);
            tmp(i,1)=geo1.wnodes.Read(i,1);
            tmp(i,2)=geo1.wnodes.Read(i,2);
        };
        geo1.wnodes.SetDimensions(nwall,3);
        geo1.wnodes=tmp;
    };
    }

    //Auxiliary vector for nodes with slip condition
    {
    unsigned int nslip=0;
    Vector<int> vhelp(geo1.nd.GetNnodes());
    vhelp.Initialize(-1);
    geo1.slipnodes.SetDimensions(geo1.bd.GetNbounds()*2,3);
    geo1.slipnodes.Initialize(-1);
    int nn,j;
    for (unsigned int i=0; i<geo1.bd.GetNbounds(); i++)
    {
        if (geo1.bd.ReadBC(i)==3)
        {
            nn=geo1.bd.ReadNode(i,0);
            j=vhelp(nn);
            if (j==-1)
            {
                nslip++;
                geo1.slipnodes(nslip-1,0) = nn;
                geo1.slipnodes(nslip-1,1) = i;
                vhelp(nn) = nslip-1;
            }
            else
            {
                geo1.slipnodes(j,2) = i;
            };

            nn=geo1.bd.ReadNode(i,1);
            j=vhelp(nn);
            if (j==-1)
            {
                nslip++;
                geo1.slipnodes(nslip-1,0) = nn;
                geo1.slipnodes(nslip-1,1) = i;
                vhelp(nn) = nslip-1;
            }
            else
            {
                geo1.slipnodes(j,2) = i;
            };
        };
    };
    //Resize slipnodes
    if (nslip==0)
    {
        geo1.slipnodes.SetDimensions(0,0);
    }
    else
    {
        matrix<int> tmp(nslip,3);
        for (unsigned int i=0; i<nslip; i++)
        {
            tmp(i,0)=geo1.slipnodes.Read(i,0);
            tmp(i,1)=geo1.slipnodes.Read(i,1);
            tmp(i,2)=geo1.slipnodes.Read(i,2);
        };
        geo1.slipnodes.SetDimensions(nslip,3);
        geo1.slipnodes=tmp;
    };

    //Calculate normals for each node with slip condition
    geo1.slipnormal.SetDimensions(nslip,2);
    int bd1,bd2;
    nprec aux;
    for (unsigned int i=0; i<nslip; i++)
    {
        bd1=geo1.slipnodes.Read(i,1);
        bd2=geo1.slipnodes.Read(i,2);
        if (bd2!=-1)
        {
            if ((geo1.bd.ReadNormal(bd1,0)*geo1.bd.ReadNormal(bd2,0)+geo1.bd.ReadNormal(bd1,1)*geo1.bd.ReadNormal(bd2,1))>0.05)
            {
                geo1.slipnormal(i,0)=geo1.bd.ReadNormal(bd1,0)+geo1.bd.ReadNormal(bd2,0);
                geo1.slipnormal(i,1)=geo1.bd.ReadNormal(bd1,1)+geo1.bd.ReadNormal(bd2,1);
                aux=std::sqrt(geo1.slipnormal.Read(i,0)*geo1.slipnormal.Read(i,0)+geo1.slipnormal.Read(i,1)*geo1.slipnormal.Read(i,1));
                geo1.slipnormal(i,0)/=aux;
                geo1.slipnormal(i,1)/=aux;
            } else
            {
                geo1.slipnormal(i,0)=par.slipnx;
                geo1.slipnormal(i,1)=par.slipny;
            };
        } else
        {
            geo1.slipnormal(i,0)=(nprec)0.0;//geo1.bd.ReadNormal(bd1,0);
            geo1.slipnormal(i,1)=(nprec)0.0;//geo1.bd.ReadNormal(bd1,1);
        };
    };
    }

     //Auxiliary vector for nodes with open condition
    {
    unsigned int nopen=0;
    Vector<int> vhelp(geo1.nd.GetNnodes());
    vhelp.Initialize(-1);
    geo1.opennodes.SetDimensions(geo1.bd.GetNbounds()*2,3);
    geo1.opennodes.Initialize(-1);
    int nn,j;
    for (unsigned int i=0; i<geo1.bd.GetNbounds(); i++)
    {
        if (geo1.bd.ReadBC(i)==1)
        {
            nn=geo1.bd.ReadNode(i,0);
            j=vhelp(nn);
            if (j==-1)
            {
                nopen++;
                geo1.opennodes(nopen-1,0) = nn;
                geo1.opennodes(nopen-1,1) = i;
                vhelp(nn) = nopen-1;
            }
            else
            {
                geo1.opennodes(j,2) = i;
            };

            nn=geo1.bd.ReadNode(i,1);
            j=vhelp(nn);
            if (j==-1)
            {
                nopen++;
                geo1.opennodes(nopen-1,0) = nn;
                geo1.opennodes(nopen-1,1) = i;
                vhelp(nn) = nopen-1;
            }
            else
            {
                geo1.opennodes(j,2) = i;
            };
        };
    };
    //Resize opennodes
    if (nopen==0)
    {
        geo1.opennodes.SetDimensions(0,0);
        par.open=false;
    }
    else
    {
        par.open=true;
        matrix<int> tmp(nopen,3);
        for (unsigned int i=0; i<nopen; i++)
        {
            tmp(i,0)=geo1.opennodes.Read(i,0);
            tmp(i,1)=geo1.opennodes.Read(i,1);
            tmp(i,2)=geo1.opennodes.Read(i,2);
        };
        geo1.opennodes.SetDimensions(nopen,3);
        geo1.opennodes=tmp;
    };

    //Calculate normals for each node with open condition
    geo1.opennormal.SetDimensions(nopen,2);
    int bd1,bd2;
    nprec aux;
    for (unsigned int i=0; i<nopen; i++)
    {
        bd1=geo1.opennodes.Read(i,1);
        bd2=geo1.opennodes.Read(i,2);
        if (bd2!=-1)
        {
            if ((geo1.bd.ReadNormal(bd1,0)*geo1.bd.ReadNormal(bd2,0)+geo1.bd.ReadNormal(bd1,1)*geo1.bd.ReadNormal(bd2,1))>0.05)
            {
                geo1.opennormal(i,0)=geo1.bd.ReadNormal(bd1,0)+geo1.bd.ReadNormal(bd2,0);
                geo1.opennormal(i,1)=geo1.bd.ReadNormal(bd1,1)+geo1.bd.ReadNormal(bd2,1);
                aux=std::sqrt(geo1.opennormal(i,0)*geo1.opennormal(i,0)+geo1.opennormal(i,1)*geo1.opennormal(i,1));
                geo1.opennormal(i,0)/=aux;
                geo1.opennormal(i,1)/=aux;
            } else
            {
                geo1.opennormal(i,0)=par.slipnx;
                geo1.opennormal(i,1)=par.slipny;
            };
        } else
        {
            geo1.opennormal(i,0)=geo1.bd.ReadNormal(bd1,0);
            geo1.opennormal(i,1)=geo1.bd.ReadNormal(bd1,1);
        };

    };
    }
    if (par.nsp>0)
    {
        bool encontrado;
        nprec eps2=1.0e-5;
        for (unsigned int i=0; i<par.nsp; i++)
        {
            encontrado=false;
            for (unsigned int j=0; j<geo1.slipnodes.GetRows(); j++)
            {
                if ((par.SPc.Read(i,0)<geo1.nd.ReadC(geo1.slipnodes.Read(j,0),0)+eps2)&&(par.SPc.Read(i,0)>geo1.nd.ReadC(geo1.slipnodes.Read(j,0),0)-eps2)&&\
                    (par.SPc.Read(i,1)<geo1.nd.ReadC(geo1.slipnodes.Read(j,0),1)+eps2)&&(par.SPc.Read(i,1)>geo1.nd.ReadC(geo1.slipnodes.Read(j,0),1)-eps2))
                {
                    encontrado=true;
                    geo1.slipnormal(j,0)=par.SPn.Read(i,0);// /std::sqrt(par.SPn.Read(i,0)*par.SPn.Read(i,0)+par.SPn.Read(i,1)*par.SPn.Read(i,1));
                    geo1.slipnormal(j,1)=par.SPn.Read(i,1);// /std::sqrt(par.SPn.Read(i,0)*par.SPn.Read(i,0)+par.SPn.Read(i,1)*par.SPn.Read(i,1));;
                    break;
                };
            };
            if (!encontrado)
            {
                std::cout << "Special Point " << i << " not found" << std::endl;
            };
        };
    };

    //Corners with outlet condition, we assign them 0 velocity condition
    {
    unsigned int corner=0;
    unsigned int nout=0;
    Vector<int> vhelp(geo1.nd.GetNnodes());
    vhelp.Initialize(-1);
    matrix<int> mhelp(geo1.bd.GetNbounds()*2,3);
    mhelp.Initialize(-1);
    int nn,j;
    for (unsigned int i=0; i<geo1.bd.GetNbounds(); i++)
    {
        if (geo1.bd.ReadBC(i)==5)
        {
            nn=geo1.bd.ReadNode(i,0);
            j=vhelp(nn);
            if (j==-1)
            {
                nout++;
                mhelp(nout-1,0) = nn;
                mhelp(nout-1,1) = i;
                vhelp(nn) = nout-1;
            }
            else
            {
                mhelp(j,2) = i;
            };

            nn=geo1.bd.ReadNode(i,1);
            j=vhelp(nn);
            if (j==-1)
            {
                nout++;
                mhelp(nout-1,0) = nn;
                mhelp(nout-1,1) = i;
                vhelp(nn) = nout-1;
            }
            else
            {
                mhelp(j,2) = i;
            };
        };
    };

    //Resize
    if (nout==0)
    {
        geo1.nodesu0.SetDimensions(0);
    }
    else
    {
        Vector<unsigned int> tmp(nout);
        //Calculate normals for each node with open condition
        int bd1,bd2;
        for (unsigned int i=0; i<nout; i++)
        {
            bd1=mhelp.Read(i,1);
            bd2=mhelp.Read(i,2);
            if (bd2!=-1)
            {
                if ((geo1.bd.ReadNormal(bd1,0)*geo1.bd.ReadNormal(bd2,0)+geo1.bd.ReadNormal(bd1,1)*geo1.bd.ReadNormal(bd2,1))<0.05)
                {
                    corner++;
                    tmp(corner-1)=mhelp.Read(i,0);
                };
            };
        };
        geo1.nodesu0.SetDimensions(corner);
        for(unsigned int i=0; i<corner; i++)
        {
            geo1.nodesu0(i)=tmp(i);
        };
    };
    }
}
