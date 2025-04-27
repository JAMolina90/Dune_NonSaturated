#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include <stdlib.h>
#include <limits>
#include "Globals.h"
#include "matrix.h"
#include "myvector.h"
#include "geometry.h"
#include "parameter.h"

void Getgeometry(geometry& geo, parameter& par, matrix<nprec>& unknowns)
{
    if(par.GInput==0)
    {
        std::ifstream readfile("nsdat10_geome");
        assert(readfile.is_open());
        unsigned int nelem,npoint,nbound;

        //nelem, npoint and nbound reading
        {
            for (int i=0 ; i<3 ; i++)
            {
                readfile.ignore(1000,'\n');
            };

            readfile >> nelem >> npoint >> nbound;
            geo.elem.SetNelem(nelem);
            geo.nd.SetNnodes(npoint);
            geo.bd.SetNbounds(nbound);
            unknowns.SetDimensions(npoint,4);
        }

        //Elements reading
        {
            for (int i=0 ; i<2 ; i++)
            {
                readfile.ignore(1000,'\n');
            };

            unsigned int tmp1, tmp2, tmp3, tmp4;
            for (unsigned int i=0 ; i<nelem ; i++)
            {
                readfile >> tmp1 >> tmp2 >> tmp3 >> tmp4;
                readfile.ignore(1000,'\n');
                geo.elem.SetNode(tmp1-1,0,tmp2-1);
                geo.elem.SetNode(tmp1-1,1,tmp3-1);
                geo.elem.SetNode(tmp1-1,2,tmp4-1);
            };
        }

        //Nodes reading
        {
            readfile.ignore(1000,'\n');
            unsigned int tmp1;
            nprec tmp2, tmp3;
            for (unsigned int i=0; i<npoint; i++)
            {
                readfile >> tmp1 >> tmp2 >> tmp3;
                readfile.ignore(1000,'\n');
                geo.nd.SetX(tmp1-1,tmp2);
                geo.nd.SetY(tmp1-1,tmp3);
            };
        }

        //Unknowns reading
        {
            readfile.ignore(1000,'\n');
            int tmp1;
            readfile >> tmp1;
            if (tmp1==0)
            {
                unknowns.Initialize((nprec)0.0);
            };
        }

        //Bounds reading
        {
            for (int i=0; i<2; i++) {readfile.ignore(1000,'\n');};
            unsigned int tmp1, tmp2, tmp3, tmp4;
            for (unsigned int i=0; i<nbound; i++)
            {
                readfile >> tmp1 >> tmp2 >> tmp3 >> tmp4;
                readfile.ignore(1000,'\n');
                geo.bd.SetNode(i,0,tmp1-1);
                geo.bd.SetNode(i,1,tmp2-1);
                geo.bd.SetEO(i,tmp3-1);
                geo.bd.SetBC(i,tmp4);
            };
        }
        readfile.close();
    }
    else
    {
        std::ifstream readfile("MESH.unv");
        if (!readfile.is_open()){std::cout << "Mesh file not found" << std::endl;exit(1);};
        Vector<unsigned int> bdname;

        int check;
        for (int i=1; i<=18; i++)
        {
            readfile.ignore(1000,'\n');
        };
        readfile >> check;
        if (check!=2411)
        {
            std::cout << "UNV file not expected (1)" << std::endl;exit(1);
        };

        //Nodes count
        int a; //First number of the row
        unsigned int npoint=0;
        readfile >> a;
        while (a!=-1)
        {
            readfile.ignore(1000,'\n');
            readfile.ignore(1000,'\n');
            readfile >> a;
            npoint+=1;
        };

        readfile.ignore(1000,'\n');
        readfile.ignore(1000,'\n');
        readfile >> check;
        if (check!=2412)
        {
            std::cout << "UNV file not expected (2)" << std::endl;exit(1);
        };

        //Elem and bound count
        int n[8];
        unsigned int nbound=0;
        unsigned int nelem=0;
        readfile >> n[0];
        while (n[0]!=-1)
        {
            for (int i=1; i<6; i++){readfile >> n[i];};

            switch(n[5])
            {
                case 2:
                    nbound+=1;
                    readfile.ignore(1000,'\n');
                    readfile.ignore(1000,'\n');
                    readfile.ignore(1000,'\n');
                    readfile >> n[0];
                    break;
                case 3:
                    nelem+=1;
                    readfile.ignore(1000,'\n');
                    readfile.ignore(1000,'\n');
                    readfile >> n[0];
                    break;
                default:
                    std::cout << "UNV file not expected (3) " << n[0] <<std::endl;exit(1);
            };
        };

        geo.elem.SetNelem(nelem);
        geo.nd.SetNnodes(npoint);
        geo.bd.SetNbounds(nbound);
        bdname.SetDimensions(nbound);
        unknowns.SetDimensions(npoint,4);
        unknowns.Initialize((nprec)0.0);

        readfile.seekg(std::ios::beg);

        for (int i=1; i<=19; i++)
        {
            readfile.ignore(1000,'\n');
        };

        //Nodes Storage
        for (unsigned int i=0; i<npoint; i++)
        {
            readfile.ignore(1000,'\n');
            readfile >> geo.nd.SetX(i) >> geo.nd.SetY(i);
            readfile.ignore(1000,'\n');
        };

        for (int i=0; i<3; i++) {readfile.ignore(1000,'\n');};

        //Bound Storage
        unsigned int aux1,aux2,aux3;
        for (unsigned int i=0; i<nbound; i++)
        {
            readfile >> aux1;
            bdname(i)=aux1-1;
            readfile.ignore(1000,'\n');
            readfile.ignore(1000,'\n');
            readfile >> aux1 >> aux2;
            readfile.ignore(1000,'\n');
            geo.bd.SetNode(i,0,aux1-1);
            geo.bd.SetNode(i,1,aux2-1);
        };

        //Elements Storage
        for (unsigned int i=0; i<nelem; i++)
        {
            readfile.ignore(1000,'\n');
            readfile >> aux1 >> aux2 >> aux3;
            readfile.ignore(1000,'\n');
            geo.elem.SetNode(i,0,aux1-1);
            geo.elem.SetNode(i,1,aux2-1);
            geo.elem.SetNode(i,2,aux3-1);
        };
        for (int i=0; i<2; i++) {readfile.ignore(1000,'\n');};
        readfile >> check;
        if (check!=2467)
        {
            std::cout << "UNV file not expected (4)" << std::endl;exit(1);
        };

        //Check Element Orientation
        {
        Vector<nprec> X(3);
        Vector<nprec> Y(3);
        nprec twoarea;
        unsigned int na,nb;
        for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
        {
            X(0)=geo.nd.ReadC(geo.elem.ReadNode(i,0),0);
            X(1)=geo.nd.ReadC(geo.elem.ReadNode(i,1),0);
            X(2)=geo.nd.ReadC(geo.elem.ReadNode(i,2),0);
            Y(0)=geo.nd.ReadC(geo.elem.ReadNode(i,0),1);
            Y(1)=geo.nd.ReadC(geo.elem.ReadNode(i,1),1);
            Y(2)=geo.nd.ReadC(geo.elem.ReadNode(i,2),1);
            twoarea=(X(1)-X(0))*(Y(2)-Y(0))-(X(2)-X(0))*(Y(1)-Y(0));
            if (twoarea<(nprec)0.0)
            {
                na=geo.elem.ReadNode(i,1);
                nb=geo.elem.ReadNode(i,2);
                geo.elem.SetNode(i,1,nb);
                geo.elem.SetNode(i,2,na);
            };
        };

        }

        //Boundary Condition Storage
        int aux[4];
        int name;
        int ncomprob;
        unsigned int nBC=0;
        readfile >> check;
        while (check!=-1)
        {
            for (int j=0; j<7; j++) {readfile >> n[j];};
            nBC+=n[6];
            readfile >> name;
            ncomprob=0;
            for (int j=0; j<n[6]; j++)
            {
                readfile >> aux[0] >> aux[1] >> aux[2] >> aux[3];
                for (unsigned int k=0; k<nbound; k++)
                {
                    if (bdname(k)==(unsigned int)aux[1]-1)
                    {
                        geo.bd.SetBC(k,name);
                        ncomprob++;
                        break;
                    };
                };
            };
            if (ncomprob!=n[6]){std::cout << "Boundaries not found" << std::endl; exit(1);};
            readfile >> check;
        };
        if (nBC!=nbound) {std::cout << "Boundary Conditions not match with bound number" << std::endl;exit(1);};

        //Element owned of each bound and bound orientation correction
        unsigned int n0,n1,auxn;
        ncomprob=0;
        for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
        {
            n0=geo.bd.ReadNode(i,0);
            n1=geo.bd.ReadNode(i,1);
            for (unsigned int j=0; j<geo.elem.GetNelem(); j++)
            {
                if ((n0==geo.elem.ReadNode(j,0)||n0==geo.elem.ReadNode(j,1)||n0==geo.elem.ReadNode(j,2))\
                    &&(n1==geo.elem.ReadNode(j,0)||n1==geo.elem.ReadNode(j,1)||n1==geo.elem.ReadNode(j,2)))
                {
                    geo.bd.SetEO(i,j);
                    ncomprob++;
                    if (!((n0==geo.elem.ReadNode(j,0)&&n1==geo.elem.ReadNode(j,1))||\
                          (n0==geo.elem.ReadNode(j,1)&&n1==geo.elem.ReadNode(j,2))||\
                          (n0==geo.elem.ReadNode(j,2)&&n1==geo.elem.ReadNode(j,0))))
                    {
                        auxn=geo.bd.ReadNode(i,0);
                        geo.bd.SetNode(i,0,geo.bd.ReadNode(i,1));
                        geo.bd.SetNode(i,1,auxn);
                    };
                    break;
                };
            };
        };
        if ((unsigned int)ncomprob!=geo.bd.GetNbounds()){std::cout << "There are boundaries without their elements" << std::endl; exit(1);};
        readfile.close();
    };

    //Alpha: initial value reading
    if (!par.restart)
    {
        nprec df;
        int aux;
        nprec c1x,c1y,c2x,c2y;
        nprec eps=(nprec)1e-9;
        std::ifstream readfile("alpha_t0.dat");
        if (!readfile.is_open()){std::cout << "Initial alpha values file not found" << std::endl;exit(1);};
        readfile >> df;
        readfile.ignore(1000,'\n');
        for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
        {
            unknowns(i,3)=df;
        };
        readfile >> aux;
        readfile.ignore(1000,'\n');
        for (int i=0; i<aux; i++)
        {
            readfile >> c1x >> c1y;
            readfile.ignore(1000,'\n');
            readfile >> c2x >> c2y;
            readfile.ignore(1000,'\n');
            readfile >> df;
            readfile.ignore(1000,'\n');
            for (unsigned int j=0; j<geo.nd.GetNnodes(); j++)
            {
                if ((geo.nd.ReadC(j,0)>=(c1x-eps)) && (geo.nd.ReadC(j,0)<=(c2x+eps)) && (geo.nd.ReadC(j,1)>=(c1y-eps)) && (geo.nd.ReadC(j,1)<=(c2y+eps)))
                {
                    unknowns(j,3)=df;
                };
            };
            par.f1height=c2y;
        };
        readfile.close();
    };

    if (par.restart)
    {
        std::ifstream mf2 ("restart.bin", std::ios::in | std::ios::binary);
        if (!mf2.is_open()){std::cout << "Restart file not found" << std::endl;exit(1);};
        mf2.read( (char*)&unknowns[0], sizeof(nprec)*unknowns.GetRows()*unknowns.GetColumns());
        mf2.close();
    };

    if((unsigned int)par.itpcgmax<=geo.nd.GetNnodes()){par.itpcgmax=(int)geo.nd.GetNnodes();};

    //Normal calculation
    nprec dx,dy;
    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        dx = geo.nd.ReadC(geo.bd.ReadNode(i,1),0)-geo.nd.ReadC(geo.bd.ReadNode(i,0),0);
        dy = geo.nd.ReadC(geo.bd.ReadNode(i,1),1)-geo.nd.ReadC(geo.bd.ReadNode(i,0),1);
        geo.bd.SetNormal(i,0,dy/(std::sqrt(dx*dx+dy*dy)));
        geo.bd.SetNormal(i,1,-dx/(std::sqrt(dx*dx+dy*dy)));
        geo.bd.SetLen(i,std::sqrt(dx*dx+dy*dy));
    };

    //Auxiliary vector for nodes with wall condition
    {
    unsigned int nwall=0;
    Vector<int> vhelp(geo.nd.GetNnodes());
    vhelp.Initialize(-1);
    geo.wnodes.SetDimensions(geo.bd.GetNbounds()*2,3);
    geo.wnodes.Initialize(-1);
    int nn,j;
    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        if (geo.bd.ReadBC(i)==2||geo.bd.ReadBC(i)==6||geo.bd.ReadBC(i)==20)
        {
            nn=geo.bd.ReadNode(i,0);
            j=vhelp(nn);
            if (j==-1)
            {
                nwall++;
                geo.wnodes(nwall-1,0) = nn;
                geo.wnodes(nwall-1,1) = i;
                vhelp(nn) = nwall-1;
            }
            else
            {
                geo.wnodes(j,2) = i;
            };

            nn=geo.bd.ReadNode(i,1);
            j=vhelp(nn);
            if (j==-1)
            {
                nwall++;
                geo.wnodes(nwall-1,0) = nn;
                geo.wnodes(nwall-1,1) = i;
                vhelp(nn) = nwall-1;
            }
            else
            {
                geo.wnodes(j,2) = i;
            };
        };
    };

    int ib1,ib2;
    nprec ach;
    for (unsigned int i=0; i<nwall; i++)
    {
        ib1=geo.wnodes.Read(i,1);
        ib2=geo.wnodes.Read(i,2);
        if ((ib1!=-1)&&(ib2!=-1))
        {
            ach=geo.bd.ReadNormal(ib1,0)*geo.bd.ReadNormal(ib2,0)+geo.bd.ReadNormal(ib1,1)*geo.bd.ReadNormal(ib2,1);
            if (ach<(nprec)-0.2)
            {
                std::cout << "Node " << geo.wnodes.Read(i,0)+1 << " is excessive acute edge" << std::endl;
            };
        };
    };
    //Resize wnodes
    if (nwall==0)
    {
        geo.wnodes.SetDimensions(0,0);
    }
    else
    {
        matrix<int> tmp(nwall,3);
        for (unsigned int i=0; i<nwall; i++)
        {
            tmp(i,0)=geo.wnodes.Read(i,0);
            tmp(i,1)=geo.wnodes.Read(i,1);
            tmp(i,2)=geo.wnodes.Read(i,2);
        };
        geo.wnodes.SetDimensions(nwall,3);
        geo.wnodes=tmp;
    };
    }

    //Auxiliary vector for nodes with slip condition
    {
    unsigned int nslip=0;
    Vector<int> vhelp(geo.nd.GetNnodes());
    vhelp.Initialize(-1);
    geo.slipnodes.SetDimensions(geo.bd.GetNbounds()*2,3);
    geo.slipnodes.Initialize(-1);
    int nn,j;
    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        if (geo.bd.ReadBC(i)==3)
        {
            nn=geo.bd.ReadNode(i,0);
            j=vhelp(nn);
            if (j==-1)
            {
                nslip++;
                geo.slipnodes(nslip-1,0) = nn;
                geo.slipnodes(nslip-1,1) = i;
                vhelp(nn) = nslip-1;
            }
            else
            {
                geo.slipnodes(j,2) = i;
            };

            nn=geo.bd.ReadNode(i,1);
            j=vhelp(nn);
            if (j==-1)
            {
                nslip++;
                geo.slipnodes(nslip-1,0) = nn;
                geo.slipnodes(nslip-1,1) = i;
                vhelp(nn) = nslip-1;
            }
            else
            {
                geo.slipnodes(j,2) = i;
            };
        };
    };
    //Resize slipnodes
    if (nslip==0)
    {
        geo.slipnodes.SetDimensions(0,0);
    }
    else
    {
        matrix<int> tmp(nslip,3);
        for (unsigned int i=0; i<nslip; i++)
        {
            tmp(i,0)=geo.slipnodes.Read(i,0);
            tmp(i,1)=geo.slipnodes.Read(i,1);
            tmp(i,2)=geo.slipnodes.Read(i,2);
        };
        geo.slipnodes.SetDimensions(nslip,3);
        geo.slipnodes=tmp;
    };

    //Calculate normals for each node with slip condition
    geo.slipnormal.SetDimensions(nslip,2);
    int bd1,bd2;
    nprec aux;
    for (unsigned int i=0; i<nslip; i++)
    {
        bd1=geo.slipnodes.Read(i,1);
        bd2=geo.slipnodes.Read(i,2);
        if (bd2!=-1)
        {
            if ((geo.bd.ReadNormal(bd1,0)*geo.bd.ReadNormal(bd2,0)+geo.bd.ReadNormal(bd1,1)*geo.bd.ReadNormal(bd2,1))>0.05)
            {
                geo.slipnormal(i,0)=geo.bd.ReadNormal(bd1,0)+geo.bd.ReadNormal(bd2,0);
                geo.slipnormal(i,1)=geo.bd.ReadNormal(bd1,1)+geo.bd.ReadNormal(bd2,1);
                aux=std::sqrt(geo.slipnormal.Read(i,0)*geo.slipnormal.Read(i,0)+geo.slipnormal.Read(i,1)*geo.slipnormal.Read(i,1));
                geo.slipnormal(i,0)/=aux;
                geo.slipnormal(i,1)/=aux;
            } else
            {
                geo.slipnormal(i,0)=par.slipnx;
                geo.slipnormal(i,1)=par.slipny;
            };
        } else
        {
            geo.slipnormal(i,0)=(nprec)0.0;//geo.bd.ReadNormal(bd1,0);
            geo.slipnormal(i,1)=(nprec)0.0;//geo.bd.ReadNormal(bd1,1);
        };
    };
    }

     //Auxiliary vector for nodes with open condition
    {
    unsigned int nopen=0;
    Vector<int> vhelp(geo.nd.GetNnodes());
    vhelp.Initialize(-1);
    geo.opennodes.SetDimensions(geo.bd.GetNbounds()*2,3);
    geo.opennodes.Initialize(-1);
    int nn,j;
    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        if (geo.bd.ReadBC(i)==1)
        {
            nn=geo.bd.ReadNode(i,0);
            j=vhelp(nn);
            if (j==-1)
            {
                nopen++;
                geo.opennodes(nopen-1,0) = nn;
                geo.opennodes(nopen-1,1) = i;
                vhelp(nn) = nopen-1;
            }
            else
            {
                geo.opennodes(j,2) = i;
            };

            nn=geo.bd.ReadNode(i,1);
            j=vhelp(nn);
            if (j==-1)
            {
                nopen++;
                geo.opennodes(nopen-1,0) = nn;
                geo.opennodes(nopen-1,1) = i;
                vhelp(nn) = nopen-1;
            }
            else
            {
                geo.opennodes(j,2) = i;
            };
        };
    };
    //Resize opennodes
    if (nopen==0)
    {
        geo.opennodes.SetDimensions(0,0);
        par.open=false;
    }
    else
    {
        par.open=true;
        matrix<int> tmp(nopen,3);
        for (unsigned int i=0; i<nopen; i++)
        {
            tmp(i,0)=geo.opennodes.Read(i,0);
            tmp(i,1)=geo.opennodes.Read(i,1);
            tmp(i,2)=geo.opennodes.Read(i,2);
        };
        geo.opennodes.SetDimensions(nopen,3);
        geo.opennodes=tmp;
    };

    //Calculate normals for each node with open condition
    geo.opennormal.SetDimensions(nopen,2);
    int bd1,bd2;
    nprec aux;
    for (unsigned int i=0; i<nopen; i++)
    {
        bd1=geo.opennodes.Read(i,1);
        bd2=geo.opennodes.Read(i,2);
        if (bd2!=-1)
        {
            if ((geo.bd.ReadNormal(bd1,0)*geo.bd.ReadNormal(bd2,0)+geo.bd.ReadNormal(bd1,1)*geo.bd.ReadNormal(bd2,1))>0.05)
            {
                geo.opennormal(i,0)=geo.bd.ReadNormal(bd1,0)+geo.bd.ReadNormal(bd2,0);
                geo.opennormal(i,1)=geo.bd.ReadNormal(bd1,1)+geo.bd.ReadNormal(bd2,1);
                aux=std::sqrt(geo.opennormal(i,0)*geo.opennormal(i,0)+geo.opennormal(i,1)*geo.opennormal(i,1));
                geo.opennormal(i,0)/=aux;
                geo.opennormal(i,1)/=aux;
            } else
            {
                geo.opennormal(i,0)=par.slipnx;
                geo.opennormal(i,1)=par.slipny;
            };
        } else
        {
            geo.opennormal(i,0)=geo.bd.ReadNormal(bd1,0);
            geo.opennormal(i,1)=geo.bd.ReadNormal(bd1,1);
        };

    };
    }
    if (par.nsp>0)
    {
        bool encontrado;
        nprec eps=1.0e-5;
        for (unsigned int i=0; i<par.nsp; i++)
        {
            encontrado=false;
            for (unsigned int j=0; j<geo.slipnodes.GetRows(); j++)
            {
                if ((par.SPc.Read(i,0)<geo.nd.ReadC(geo.slipnodes.Read(j,0),0)+eps)&&(par.SPc.Read(i,0)>geo.nd.ReadC(geo.slipnodes.Read(j,0),0)-eps)&&\
                    (par.SPc.Read(i,1)<geo.nd.ReadC(geo.slipnodes.Read(j,0),1)+eps)&&(par.SPc.Read(i,1)>geo.nd.ReadC(geo.slipnodes.Read(j,0),1)-eps))
                {
                    encontrado=true;
                    geo.slipnormal(j,0)=par.SPn.Read(i,0);// /std::sqrt(par.SPn.Read(i,0)*par.SPn.Read(i,0)+par.SPn.Read(i,1)*par.SPn.Read(i,1));
                    geo.slipnormal(j,1)=par.SPn.Read(i,1);// /std::sqrt(par.SPn.Read(i,0)*par.SPn.Read(i,0)+par.SPn.Read(i,1)*par.SPn.Read(i,1));;
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
    Vector<int> vhelp(geo.nd.GetNnodes());
    vhelp.Initialize(-1);
    matrix<int> mhelp(geo.bd.GetNbounds()*2,3);
    mhelp.Initialize(-1);
    int nn,j;
    for (unsigned int i=0; i<geo.bd.GetNbounds(); i++)
    {
        if (geo.bd.ReadBC(i)==5)
        {
            nn=geo.bd.ReadNode(i,0);
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

            nn=geo.bd.ReadNode(i,1);
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
        geo.nodesu0.SetDimensions(0);
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
                if ((geo.bd.ReadNormal(bd1,0)*geo.bd.ReadNormal(bd2,0)+geo.bd.ReadNormal(bd1,1)*geo.bd.ReadNormal(bd2,1))<0.05)
                {
                    corner++;
                    tmp(corner-1)=mhelp.Read(i,0);
                };
            };
        };
        geo.nodesu0.SetDimensions(corner);
        for(unsigned int i=0; i<corner; i++)
        {
            geo.nodesu0(i)=tmp(i);
        };
    };


    }
    /*
    //Set yref
    par.yref=std::numeric_limits<nprec>::max();
    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        par.yref=std::min(par.yref,geo.nd.ReadC(i,1));
    };
    */

}
