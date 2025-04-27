#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <limits>
#include "Globals.h"
#include "geometry.h"
#include "myvector.h"
#include "matrix.h"
#include "parameter.h"
#include "Neighbour.h"

void FCTU(const parameter& par, const geometry& geo, const matrix<nprec>& SFD, const Vector<nprec>& DMAT, const Vector<nprec>& varPre1, const Vector<nprec>& varPre2, const Vector<nprec>& varL1, const Vector<nprec>& varL2, const Vector<nprec>& varH1, const Vector<nprec>& varH2, matrix<nprec>& FS, const Vector<nprec>& Phi)
{
    bool sync=true;
    bool prelimiter=false;

    matrix<nprec> CP(geo.nd.GetNnodes(),2);
    CP.Initialize((nprec)0.0);
    matrix<nprec> CM(geo.nd.GetNnodes(),2);
    CM.Initialize((nprec)0.0);
    matrix<nprec> ULE(geo.nd.GetNnodes(),2);
    ULE.Initialize((nprec)0.0);

    Vector<int> NULE(geo.nd.GetNnodes());
    NULE.Initialize(0);
    matrix<nprec> ASUM(geo.nd.GetNnodes(),2);
    ASUM.Initialize((nprec)0.0);
    matrix<nprec> AREST(geo.nd.GetNnodes(),2);
    AREST.Initialize((nprec)0.0);
    matrix<nprec> AC(geo.nd.GetNnodes(),2);
    AC.Initialize((nprec)0.0);
    matrix<nprec> DBCU(3,2);
    DBCU.Initialize((nprec)0.0);
    matrix<nprec> AEL(3,2);
    matrix<nprec> DBE(3,2);

    nprec betap, betan;
    matrix<nprec> BMAX(geo.nd.GetNnodes(),2);
    matrix<nprec> BMIN(geo.nd.GetNnodes(),2);
    unsigned int i1,i2,i3;
    //nprec half=(nprec)0.5;
    nprec dif=(nprec)1.0/((nprec)30.0*par.deps);
    //HAY QUE RETOQUETEAR ESTE PARÁMETRO
    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        BMAX(i,0)=std::max(varPre1.Read(i),varL1.Read(i));
        BMIN(i,0)=std::min(varPre1.Read(i),varL1.Read(i));
        BMAX(i,1)=std::max(varPre2.Read(i),varL2.Read(i));
        BMIN(i,1)=std::min(varPre2.Read(i),varL2.Read(i));
    };

    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);
        //if ((Phi.Read(i1)<=half&&Phi.Read(i2)<=half)||(Phi.Read(i1)>half&&Phi.Read(i2)>half))
        if(std::abs(Phi.Read(i1)-Phi.Read(i2))<=dif)
        {
            BMAX(i1,0)=std::max(BMAX(i1,0),std::max(varPre1.Read(i2),varL1.Read(i2)));
            BMIN(i1,0)=std::min(BMIN(i1,0),std::min(varPre1.Read(i2),varL1.Read(i2)));
            BMAX(i1,1)=std::max(BMAX(i1,1),std::max(varPre2.Read(i2),varL2.Read(i2)));
            BMIN(i1,1)=std::min(BMIN(i1,1),std::min(varPre2.Read(i2),varL2.Read(i2)));

            BMAX(i2,0)=std::max(BMAX(i2,0),std::max(varPre1.Read(i1),varL1.Read(i1)));
            BMIN(i2,0)=std::min(BMIN(i2,0),std::min(varPre1.Read(i1),varL1.Read(i1)));
            BMAX(i2,1)=std::max(BMAX(i2,1),std::max(varPre2.Read(i1),varL2.Read(i1)));
            BMIN(i2,1)=std::min(BMIN(i2,1),std::min(varPre2.Read(i1),varL2.Read(i1)));
        };

        //if ((Phi.Read(i2)<=half&&Phi.Read(i3)<=half)||(Phi.Read(i2)>half&&Phi.Read(i3)>half))
        if(std::abs(Phi.Read(i2)-Phi.Read(i3))<=dif)
        {
            BMAX(i2,0)=std::max(BMAX(i2,0),std::max(varPre1.Read(i3),varL1.Read(i3)));
            BMIN(i2,0)=std::min(BMIN(i2,0),std::min(varPre1.Read(i3),varL1.Read(i3)));
            BMAX(i2,1)=std::max(BMAX(i2,1),std::max(varPre2.Read(i3),varL2.Read(i3)));
            BMIN(i2,1)=std::min(BMIN(i2,1),std::min(varPre2.Read(i3),varL2.Read(i3)));

            BMAX(i3,0)=std::max(BMAX(i3,0),std::max(varPre1.Read(i2),varL1.Read(i2)));
            BMIN(i3,0)=std::min(BMIN(i3,0),std::min(varPre1.Read(i2),varL1.Read(i2)));
            BMAX(i3,1)=std::max(BMAX(i3,1),std::max(varPre2.Read(i2),varL2.Read(i2)));
            BMIN(i3,1)=std::min(BMIN(i3,1),std::min(varPre2.Read(i2),varL2.Read(i2)));
        };

        //if ((Phi.Read(i1)<=half&&Phi.Read(i3)<=half)||(Phi.Read(i1)>half&&Phi.Read(i3)>half))
        if(std::abs(Phi.Read(i1)-Phi.Read(i3))<=dif)
        {
            BMAX(i1,0)=std::max(BMAX(i1,0),std::max(varPre1.Read(i3),varL1.Read(i3)));
            BMIN(i1,0)=std::min(BMIN(i1,0),std::min(varPre1.Read(i3),varL1.Read(i3)));
            BMAX(i1,1)=std::max(BMAX(i1,1),std::max(varPre2.Read(i3),varL2.Read(i3)));
            BMIN(i1,1)=std::min(BMIN(i1,1),std::min(varPre2.Read(i3),varL2.Read(i3)));

            BMAX(i3,0)=std::max(BMAX(i3,0),std::max(varPre1.Read(i1),varL1.Read(i1)));
            BMIN(i3,0)=std::min(BMIN(i3,0),std::min(varPre1.Read(i1),varL1.Read(i1)));
            BMAX(i3,1)=std::max(BMAX(i3,1),std::max(varPre2.Read(i1),varL2.Read(i1)));
            BMIN(i3,1)=std::min(BMIN(i3,1),std::min(varPre2.Read(i1),varL2.Read(i1)));
        };
    };

    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);
        DBE(0,0)=varH1.Read(i1)-varL1.Read(i1);
        DBE(1,0)=varH1.Read(i2)-varL1.Read(i2);
        DBE(2,0)=varH1.Read(i3)-varL1.Read(i3);

        DBE(0,1)=varH2.Read(i1)-varL2.Read(i1);
        DBE(1,1)=varH2.Read(i2)-varL2.Read(i2);
        DBE(2,1)=varH2.Read(i3)-varL2.Read(i3);

        if (prelimiter)
        {
            ULE(i1,0)+=(varL1.Read(i1)+varL1.Read(i2)+varL1.Read(i3))/((nprec)3.0);
            ULE(i2,0)+=(varL1.Read(i1)+varL1.Read(i2)+varL1.Read(i3))/((nprec)3.0);
            ULE(i3,0)+=(varL1.Read(i1)+varL1.Read(i2)+varL1.Read(i3))/((nprec)3.0);
            ULE(i1,1)+=(varL2.Read(i1)+varL2.Read(i2)+varL2.Read(i3))/((nprec)3.0);
            ULE(i2,1)+=(varL2.Read(i1)+varL2.Read(i2)+varL2.Read(i3))/((nprec)3.0);
            ULE(i3,1)+=(varL2.Read(i1)+varL2.Read(i2)+varL2.Read(i3))/((nprec)3.0);
            NULE(i1)++;
            NULE(i2)++;
            NULE(i3)++;
        };

        DBCU(0,0)=SFD.Read(i,6)/((nprec)6.0)*DBE(0,0);
        DBCU(1,0)=SFD.Read(i,6)/((nprec)6.0)*DBE(1,0);
        DBCU(2,0)=SFD.Read(i,6)/((nprec)6.0)*DBE(2,0);

        DBCU(0,1)=SFD.Read(i,6)/((nprec)6.0)*DBE(0,1);
        DBCU(1,1)=SFD.Read(i,6)/((nprec)6.0)*DBE(1,1);
        DBCU(2,1)=SFD.Read(i,6)/((nprec)6.0)*DBE(2,1);

        AEL(0,0)=DBCU(0,0)*DMAT.Read(i1);
        AEL(1,0)=DBCU(1,0)*DMAT.Read(i2);
        AEL(2,0)=DBCU(2,0)*DMAT.Read(i3);

        AEL(0,1)=DBCU(0,1)*DMAT.Read(i1);
        AEL(1,1)=DBCU(1,1)*DMAT.Read(i2);
        AEL(2,1)=DBCU(2,1)*DMAT.Read(i3);

        if (AEL(0,0)>(nprec)0.0)
        {
            ASUM(i1,0)+=AEL(0,0);
        }else
        {
            AREST(i1,0)-=AEL(0,0);
        };
        if (AEL(1,0)>(nprec)0.0)
        {
            ASUM(i2,0)+=AEL(1,0);
        }else
        {
            AREST(i2,0)-=AEL(1,0);
        };
        if (AEL(2,0)>(nprec)0.0)
        {
            ASUM(i3,0)+=AEL(2,0);
        }else
        {
            AREST(i3,0)-=AEL(2,0);
        };

        if (AEL(0,1)>(nprec)0.0)
        {
            ASUM(i1,1)+=AEL(0,1);
        }else
        {
            AREST(i1,1)-=AEL(0,1);
        };
        if (AEL(1,1)>(nprec)0.0)
        {
            ASUM(i2,1)+=AEL(1,1);
        }else
        {
            AREST(i2,1)-=AEL(1,1);
        };
        if (AEL(2,1)>(nprec)0.0)
        {
            ASUM(i3,1)+=AEL(2,1);
        }else
        {
            AREST(i3,1)-=AEL(2,1);
        };
    };

    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        if (prelimiter)
        {
            ULE(i,0)/=((nprec)NULE(i)+1.0e-9);
            if (varL1.Read(i)>ULE(i,0)&&(ASUM(i,0)-AREST(i,0))<((nprec)0.0))
            {
                AREST(i,0)=(nprec)0.0;
                ASUM(i,0)=(nprec)0.0;
            };
            if (varL1.Read(i)<ULE(i,0)&&(ASUM(i,0)-AREST(i,0))>((nprec)0.0))
            {
                AREST(i,0)=(nprec)0.0;
                ASUM(i,0)=(nprec)0.0;
            };

            ULE(i,1)/=((nprec)NULE(i)+1.0e-9);
            if (varL2.Read(i)>ULE(i,1)&&(ASUM(i,1)-AREST(i,1))<((nprec)0.0))
            {
                AREST(i,1)=(nprec)0.0;
                ASUM(i,1)=(nprec)0.0;
            };
            if (varL2.Read(i)<ULE(i,1)&&(ASUM(i,1)-AREST(i,1))>((nprec)0.0))
            {
                AREST(i,1)=(nprec)0.0;
                ASUM(i,1)=(nprec)0.0;
            };
        };

        betap=(BMAX(i,0)-varL1.Read(i))/(ASUM(i,0)+1.0e-9);
        betan=(varL1.Read(i)-BMIN(i,0))/(AREST(i,0)+1.0e-9);
        CP(i,0)=std::min((nprec)1.0,betap);
        CM(i,0)=std::min((nprec)1.0,betan);

        betap=(BMAX(i,1)-varL2.Read(i))/(ASUM(i,1)+1.0e-9);
        betan=(varL2.Read(i)-BMIN(i,1))/(AREST(i,1)+1.0e-9);
        CP(i,1)=std::min((nprec)1.0,betap);
        CM(i,1)=std::min((nprec)1.0,betan);

        //NULE(i)=0; no sirve para nada
    };

    nprec cpe1,cme1,cpe2,cme2;
    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        AEL.Initialize((nprec)0.0);
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);

        cpe1=std::min(std::min(CP(i1,0),CP(i2,0)),CP(i3,0));
        cme1=std::min(std::min(CM(i1,0),CM(i2,0)),CM(i3,0));
        cpe2=std::min(std::min(CP(i1,1),CP(i2,1)),CP(i3,1));
        cme2=std::min(std::min(CM(i1,1),CM(i2,1)),CM(i3,1));

        //cpe1=std::min(cpe1,cme1);
        //cme1=cpe1;
        //cpe2=std::min(cpe2,cme2);
        //cme2=cpe2;

        if (sync)
        {
            cpe1=std::min(cpe1,cpe2);
            cpe2=cpe1;
            cme1=std::min(cme1,cme2);
            cme2=cme1;
	    //Si estas lineas de abajo se quitan, el FCT introduce menos difusion
	    cpe1=std::min(cpe1,cme1);
	    cpe2=cpe1;
	    cme1=cpe1;
	    cme2=cpe1;
        };

        DBE(0,0)=varH1.Read(i1)-varL1.Read(i1);
        DBE(1,0)=varH1.Read(i2)-varL1.Read(i2);
        DBE(2,0)=varH1.Read(i3)-varL1.Read(i3);

        DBE(0,1)=varH2.Read(i1)-varL2.Read(i1);
        DBE(1,1)=varH2.Read(i2)-varL2.Read(i2);
        DBE(2,1)=varH2.Read(i3)-varL2.Read(i3);

        DBCU(0,0)=SFD.Read(i,6)/((nprec)6.0)*DBE(0,0);
        DBCU(1,0)=SFD.Read(i,6)/((nprec)6.0)*DBE(1,0);
        DBCU(2,0)=SFD.Read(i,6)/((nprec)6.0)*DBE(2,0);

        DBCU(0,1)=SFD.Read(i,6)/((nprec)6.0)*DBE(0,1);
        DBCU(1,1)=SFD.Read(i,6)/((nprec)6.0)*DBE(1,1);
        DBCU(2,1)=SFD.Read(i,6)/((nprec)6.0)*DBE(2,1);

        AEL(0,0)=DBCU(0,0)*DMAT.Read(i1);
        AEL(1,0)=DBCU(1,0)*DMAT.Read(i2);
        AEL(2,0)=DBCU(2,0)*DMAT.Read(i3);

        AEL(0,1)=DBCU(0,1)*DMAT.Read(i1);
        AEL(1,1)=DBCU(1,1)*DMAT.Read(i2);
        AEL(2,1)=DBCU(2,1)*DMAT.Read(i3);

        if (AEL(0,0)>(nprec)0.0)
        {
            AC(i1,0)+=cpe1*AEL(0,0);
        }else
        {
            AC(i1,0)+=cme1*AEL(0,0);
        };
        if (AEL(1,0)>(nprec)0.0)
        {
            AC(i2,0)+=cpe1*AEL(1,0);
        }else
        {
            AC(i2,0)+=cme1*AEL(1,0);
        };
        if (AEL(2,0)>(nprec)0.0)
        {
            AC(i3,0)+=cpe1*AEL(2,0);
        }else
        {
            AC(i3,0)+=cme1*AEL(2,0);
        };

        if (AEL(0,1)>(nprec)0.0)
        {
            AC(i1,1)+=cpe2*AEL(0,1);
        }else
        {
            AC(i1,1)+=cme2*AEL(0,1);
        };
        if (AEL(1,1)>(nprec)0.0)
        {
            AC(i2,1)+=cpe2*AEL(1,1);
        }else
        {
            AC(i2,1)+=cme2*AEL(1,1);
        };
        if (AEL(2,1)>(nprec)0.0)
        {
            AC(i3,1)+=cpe2*AEL(2,1);
        }else
        {
            AC(i3,1)+=cme2*AEL(2,1);
        };
    };

    if (prelimiter)
    {
        for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
        {
            if (varL1.Read(i)>ULE(i,0)&&AC(i,0)<(nprec)0.0)
            {
                AC(i,0)=(nprec)0.0;
            };
            if (varL1.Read(i)<ULE(i,0)&&AC(i,0)>(nprec)0.0)
            {
                AC(i,0)=(nprec)0.0;
            };

            if (varL2.Read(i)>ULE(i,1)&&AC(i,1)<(nprec)0.0)
            {
                AC(i,1)=(nprec)0.0;
            };
            if (varL2.Read(i)<ULE(i,1)&&AC(i,1)>(nprec)0.0)
            {
                AC(i,1)=(nprec)0.0;
            };
        };
    };

    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        FS(i,0)=varL1.Read(i)+AC(i,0);
        FS(i,1)=varL2.Read(i)+AC(i,1);
    };
}

void FCTU2(const geometry& geo, const matrix<nprec>& SFD, const Vector<nprec>& DMAT, const Vector<nprec>& varPre1, const Vector<nprec>& varPre2, const Vector<nprec>& varL1, const Vector<nprec>& varL2, const Vector<nprec>& varH1, const Vector<nprec>& varH2, matrix<nprec>& FS, const Vector<nprec>& dUho, const Vector<nprec>& dVho)
{
    bool sync=true;

    nprec zero=(nprec)0.0;
    matrix<nprec> CP(geo.nd.GetNnodes(),2);
    CP.Initialize(zero);
    matrix<nprec> CM(geo.nd.GetNnodes(),2);
    CM.Initialize(zero);
    matrix<nprec> ElContU (geo.elem.GetNelem(),3);
    ElContU.Initialize(zero);
    matrix<nprec> ElContV (geo.elem.GetNelem(),3);
    ElContV.Initialize(zero);
    matrix<nprec> FluxU(geo.nd.GetNnodes(),2);
    FluxU.Initialize(zero);
    matrix<nprec> FluxV(geo.nd.GetNnodes(),2);
    FluxV.Initialize(zero);
    matrix<nprec> AC(geo.nd.GetNnodes(),2);
    AC.Initialize(zero);

    nprec betap, betan, rj, dmeanu, dmeanv;
    matrix<nprec> BMAX(geo.nd.GetNnodes(),2);
    matrix<nprec> BMIN(geo.nd.GetNnodes(),2);
    unsigned int i1,i2,i3;

    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        BMAX(i,0)=std::max(varPre1.Read(i),varL1.Read(i));
        BMIN(i,0)=std::min(varPre1.Read(i),varL1.Read(i));
        BMAX(i,1)=std::max(varPre2.Read(i),varL2.Read(i));
        BMIN(i,1)=std::min(varPre2.Read(i),varL2.Read(i));
    };

    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);
        BMAX(i1,0)=std::max(BMAX(i1,0),std::max(varPre1.Read(i2),varL1.Read(i2)));
        BMIN(i1,0)=std::min(BMIN(i1,0),std::min(varPre1.Read(i2),varL1.Read(i2)));
        BMAX(i1,1)=std::max(BMAX(i1,1),std::max(varPre2.Read(i2),varL2.Read(i2)));
        BMIN(i1,1)=std::min(BMIN(i1,1),std::min(varPre2.Read(i2),varL2.Read(i2)));

        BMAX(i2,0)=std::max(BMAX(i2,0),std::max(varPre1.Read(i1),varL1.Read(i1)));
        BMIN(i2,0)=std::min(BMIN(i2,0),std::min(varPre1.Read(i1),varL1.Read(i1)));
        BMAX(i2,1)=std::max(BMAX(i2,1),std::max(varPre2.Read(i1),varL2.Read(i1)));
        BMIN(i2,1)=std::min(BMIN(i2,1),std::min(varPre2.Read(i1),varL2.Read(i1)));

        BMAX(i2,0)=std::max(BMAX(i2,0),std::max(varPre1.Read(i3),varL1.Read(i3)));
        BMIN(i2,0)=std::min(BMIN(i2,0),std::min(varPre1.Read(i3),varL1.Read(i3)));
        BMAX(i2,1)=std::max(BMAX(i2,1),std::max(varPre2.Read(i3),varL2.Read(i3)));
        BMIN(i2,1)=std::min(BMIN(i2,1),std::min(varPre2.Read(i3),varL2.Read(i3)));

        BMAX(i3,0)=std::max(BMAX(i3,0),std::max(varPre1.Read(i2),varL1.Read(i2)));
        BMIN(i3,0)=std::min(BMIN(i3,0),std::min(varPre1.Read(i2),varL1.Read(i2)));
        BMAX(i3,1)=std::max(BMAX(i3,1),std::max(varPre2.Read(i2),varL2.Read(i2)));
        BMIN(i3,1)=std::min(BMIN(i3,1),std::min(varPre2.Read(i2),varL2.Read(i2)));

        BMAX(i1,0)=std::max(BMAX(i1,0),std::max(varPre1.Read(i3),varL1.Read(i3)));
        BMIN(i1,0)=std::min(BMIN(i1,0),std::min(varPre1.Read(i3),varL1.Read(i3)));
        BMAX(i1,1)=std::max(BMAX(i1,1),std::max(varPre2.Read(i3),varL2.Read(i3)));
        BMIN(i1,1)=std::min(BMIN(i1,1),std::min(varPre2.Read(i3),varL2.Read(i3)));

        BMAX(i3,0)=std::max(BMAX(i3,0),std::max(varPre1.Read(i1),varL1.Read(i1)));
        BMIN(i3,0)=std::min(BMIN(i3,0),std::min(varPre1.Read(i1),varL1.Read(i1)));
        BMAX(i3,1)=std::max(BMAX(i3,1),std::max(varPre2.Read(i1),varL2.Read(i1)));
        BMIN(i3,1)=std::min(BMIN(i3,1),std::min(varPre2.Read(i1),varL2.Read(i1)));

    };

    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);
        rj=SFD.Read(i,6)/((nprec)24.0);
        dmeanu=dUho.Read(i1)+dUho.Read(i2)+dUho.Read(i3);
        dmeanv=dVho.Read(i1)+dVho.Read(i2)+dVho.Read(i3);
        ElContU(i,0)=DMAT.Read(i1)*(((nprec)3.0)*dUho.Read(i1)-dmeanu)*rj;
        ElContU(i,1)=DMAT.Read(i2)*(((nprec)3.0)*dUho.Read(i2)-dmeanu)*rj;
        ElContU(i,2)=DMAT.Read(i3)*(((nprec)3.0)*dUho.Read(i3)-dmeanu)*rj;
        ElContV(i,0)=DMAT.Read(i1)*(((nprec)3.0)*dVho.Read(i1)-dmeanv)*rj;
        ElContV(i,1)=DMAT.Read(i2)*(((nprec)3.0)*dVho.Read(i2)-dmeanv)*rj;
        ElContV(i,2)=DMAT.Read(i3)*(((nprec)3.0)*dVho.Read(i3)-dmeanv)*rj;

        if (ElContU(i,0)>=zero){FluxU(i1,0)+=ElContU(i,0);}else{FluxU(i1,1)-=ElContU(i,0);};
        if (ElContU(i,1)>=zero){FluxU(i2,0)+=ElContU(i,1);}else{FluxU(i2,1)-=ElContU(i,1);};
        if (ElContU(i,2)>=zero){FluxU(i3,0)+=ElContU(i,2);}else{FluxU(i3,1)-=ElContU(i,2);};
        if (ElContV(i,0)>=zero){FluxV(i1,0)+=ElContV(i,0);}else{FluxV(i1,1)-=ElContV(i,0);};
        if (ElContV(i,1)>=zero){FluxV(i2,0)+=ElContV(i,1);}else{FluxV(i2,1)-=ElContV(i,1);};
        if (ElContV(i,2)>=zero){FluxV(i3,0)+=ElContV(i,2);}else{FluxV(i3,1)-=ElContV(i,2);};
    };

    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        betap=(BMAX(i,0)-varL1.Read(i))/(FluxU(i,0)+1.0e-9);
        betan=(varL1.Read(i)-BMIN(i,0))/(FluxU(i,1)+1.0e-9);
        CP(i,0)=std::min((nprec)1.0,betap);
        CM(i,0)=std::min((nprec)1.0,betan);

        betap=(BMAX(i,1)-varL2.Read(i))/(FluxV(i,0)+1.0e-9);
        betan=(varL2.Read(i)-BMIN(i,1))/(FluxV(i,1)+1.0e-9);
        CP(i,1)=std::min((nprec)1.0,betap);
        CM(i,1)=std::min((nprec)1.0,betan);
    };

    nprec cpe1,cme1,cpe2,cme2;
    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);

        cpe1=std::min(std::min(CP(i1,0),CP(i2,0)),CP(i3,0));
        cme1=std::min(std::min(CM(i1,0),CM(i2,0)),CM(i3,0));
        cpe2=std::min(std::min(CP(i1,1),CP(i2,1)),CP(i3,1));
        cme2=std::min(std::min(CM(i1,1),CM(i2,1)),CM(i3,1));

        cpe1=std::min(cpe1,cme1);
        cme1=cpe1;
        cpe2=std::min(cpe2,cme2);
        cme2=cpe2;

        if (sync)
        {
            cpe1=std::min(cpe1,cpe2);
            cpe2=cpe1;
            cme1=std::min(cme1,cme2);
            cme2=cme1;
        };

        if (ElContU(i,0)>zero)
        {
            AC(i1,0)+=cpe1*ElContU(i,0);
        }else
        {
            AC(i1,0)+=cme1*ElContU(i,0);
        };
        if (ElContU(i,1)>zero)
        {
            AC(i2,0)+=cpe1*ElContU(i,1);
        }else
        {
            AC(i2,0)+=cme1*ElContU(i,1);
        };
        if (ElContU(i,2)>zero)
        {
            AC(i3,0)+=cpe1*ElContU(i,2);
        }else
        {
            AC(i3,0)+=cme1*ElContU(i,2);
        };

        if (ElContV(i,0)>zero)
        {
            AC(i1,1)+=cpe2*ElContV(i,0);
        }else
        {
            AC(i1,1)+=cme2*ElContV(i,0);
        };
        if (ElContV(i,1)>zero)
        {
            AC(i2,1)+=cpe2*ElContV(i,1);
        }else
        {
            AC(i2,1)+=cme2*ElContV(i,1);
        };
        if (ElContV(i,2)>zero)
        {
            AC(i3,1)+=cpe2*ElContV(i,2);
        }else
        {
            AC(i3,1)+=cme2*ElContV(i,2);
        };
    };

    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
        FS(i,0)=varL1.Read(i)+AC(i,0);
        FS(i,1)=varL2.Read(i)+AC(i,1);
    };
}
