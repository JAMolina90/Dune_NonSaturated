#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <limits>
#include "Globals.h"
#include "geometry.h"
#include "myvector.h"
#include "UpwindVar.h"

void FCT(const geometry& geo, const Vector<nprec>& varPre, UpwindVar& up)
{
    bool prelimiter=false;
    Vector<nprec> BMAX(up.LOS.GetSize());
    Vector<nprec> BMIN(up.LOS.GetSize());
    BMAX.Initialize(std::numeric_limits<nprec>::min());
    BMIN.Initialize(std::numeric_limits<nprec>::max());
    unsigned int i1,i2,i3;
    /*
    nprec totalmass=(nprec)0.0;
    for (unsigned int i=0; i<up.LOS.GetSize(); i++)
    {
        totalmass+=varPre.Read(i)*up.VOL.Read(i);
    };
    */
    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);
        BMAX(i1)=std::max(std::max(std::max(BMAX(i1),varPre.Read(i1)),std::max(varPre.Read(i2),varPre.Read(i3))),std::max(std::max(up.LOS.Read(i1),up.LOS.Read(i2)),up.LOS.Read(i3)));
        BMIN(i1)=std::min(std::min(std::min(BMIN(i1),varPre.Read(i1)),std::min(varPre.Read(i2),varPre.Read(i3))),std::min(std::min(up.LOS.Read(i1),up.LOS.Read(i2)),up.LOS.Read(i3)));
        BMAX(i2)=std::max(std::max(std::max(BMAX(i2),varPre.Read(i1)),std::max(varPre.Read(i2),varPre.Read(i3))),std::max(std::max(up.LOS.Read(i1),up.LOS.Read(i2)),up.LOS.Read(i3)));
        BMIN(i2)=std::min(std::min(std::min(BMIN(i2),varPre.Read(i1)),std::min(varPre.Read(i2),varPre.Read(i3))),std::min(std::min(up.LOS.Read(i1),up.LOS.Read(i2)),up.LOS.Read(i3)));
        BMAX(i3)=std::max(std::max(std::max(BMAX(i3),varPre.Read(i1)),std::max(varPre.Read(i2),varPre.Read(i3))),std::max(std::max(up.LOS.Read(i1),up.LOS.Read(i2)),up.LOS.Read(i3)));
        BMIN(i3)=std::min(std::min(std::min(BMIN(i3),varPre.Read(i1)),std::min(varPre.Read(i2),varPre.Read(i3))),std::min(std::min(up.LOS.Read(i1),up.LOS.Read(i2)),up.LOS.Read(i3)));
    };

    Vector<int> NULE(up.LOS.GetSize());
    NULE.Initialize(0);
    Vector<nprec> ASUM(up.LOS.GetSize());
    ASUM.Initialize((nprec)0.0);
    Vector<nprec> AREST(up.LOS.GetSize());
    AREST.Initialize((nprec)0.0);
    matrix<nprec> AC(geo.elem.GetNelem(),6);
    Vector<nprec> ACf(up.LOS.GetSize());
    AC.Initialize((nprec)0.0);
    ACf.Initialize((nprec)0.0);
    Vector<nprec> ULE(up.LOS.GetSize());
    ULE.Initialize((nprec)0.0);
    Vector<nprec> DBCU(3);
    DBCU.Initialize((nprec)0.0);
    Vector<nprec> AEL(3);
    Vector<nprec> DBE(3);
    DBE.Initialize((nprec)0.0);
    nprec massp, massn;
    massp=(nprec)0.0;
    massn=(nprec)0.0;
    Vector<nprec> massNp(up.LOS.GetSize());
    Vector<nprec> massNn(up.LOS.GetSize());
    massNp.Initialize((nprec)0.0);
    massNn.Initialize((nprec)0.0);

    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        AEL.Initialize((nprec)0.0);
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);
        DBE(0)=up.HOS(i1)-up.LOS(i1);
        DBE(1)=up.HOS(i2)-up.LOS(i2);
        DBE(2)=up.HOS(i3)-up.LOS(i3);

        if (prelimiter)
        {
            ULE(i1)+=(up.LOS(i1)+up.LOS(i2)+up.LOS(i3))/((nprec)3.0);
            ULE(i2)+=(up.LOS(i1)+up.LOS(i2)+up.LOS(i3))/((nprec)3.0);
            ULE(i3)+=(up.LOS(i1)+up.LOS(i2)+up.LOS(i3))/((nprec)3.0);
            NULE(i1)++;
            NULE(i2)++;
            NULE(i3)++;
        };

        DBCU(0)=up.GEOEL.Read(i,9)*DBE(0);
        DBCU(1)=up.GEOEL.Read(i,10)*DBE(1);
        DBCU(2)=up.GEOEL.Read(i,11)*DBE(2);

        AEL(0)+=DBCU(0);
        AEL(0)/=up.VOL.Read(i1);
        AEL(1)+=DBCU(1);
        AEL(1)/=up.VOL.Read(i2);
        AEL(2)+=DBCU(2);
        AEL(2)/=up.VOL.Read(i3);

        if (AEL(0)>(nprec)0.0)
        {
            ASUM(i1)+=AEL(0);
        }else
        {
            AREST(i1)-=AEL(0);
        };
        if (AEL(1)>(nprec)0.0)
        {
            ASUM(i2)+=AEL(1);
        }else
        {
            AREST(i2)-=AEL(1);
        };
        if (AEL(2)>(nprec)0.0)
        {
            ASUM(i3)+=AEL(2);
        }else
        {
            AREST(i3)-=AEL(2);
        };
    };

    Vector<nprec> CP(up.LOS.GetSize());
    CP.Initialize((nprec)0.0);
    Vector<nprec> CM(up.LOS.GetSize());
    CM.Initialize((nprec)0.0);
    nprec betap, betan;

    for (unsigned int i=0; i<up.LOS.GetSize(); i++)
    {
        if (prelimiter)
        {
            ULE(i)/=((nprec)NULE(i)+1.0e-7);
            if (up.LOS(i)>ULE(i)&&(ASUM(i)-AREST(i))<((nprec)0.0))
            {
                AREST(i)=(nprec)0.0;
                ASUM(i)=(nprec)0.0;
            };
            if (up.LOS(i)<ULE(i)&&(ASUM(i)-AREST(i))>((nprec)0.0))
            {
                AREST(i)=(nprec)0.0;
                ASUM(i)=(nprec)0.0;
            };
        };

        betap=(BMAX(i)-up.LOS(i))/(ASUM(i)+1.0e-7);
        betan=(up.LOS(i)-BMIN(i))/(AREST(i)+1.0e-7);

        CP(i)=std::min((nprec)1.0,betap);
        CM(i)=std::min((nprec)1.0,betan);

        NULE(i)=0;
    };

    nprec cpe,cme;
    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        AEL.Initialize((nprec)0.0);
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);

        cpe=std::min(std::min(CP(i1),CP(i2)),CP(i3));
        cme=std::min(std::min(CM(i1),CM(i2)),CM(i3));

        DBE(0)=up.HOS(i1)-up.LOS(i1);
        DBE(1)=up.HOS(i2)-up.LOS(i2);
        DBE(2)=up.HOS(i3)-up.LOS(i3);

        DBCU(0)=up.GEOEL.Read(i,9)*DBE(0);
        DBCU(1)=up.GEOEL.Read(i,10)*DBE(1);
        DBCU(2)=up.GEOEL.Read(i,11)*DBE(2);

        AEL(0)+=DBCU(0);
        AEL(0)/=up.VOL.Read(i1);
        AEL(1)+=DBCU(1);
        AEL(1)/=up.VOL.Read(i2);
        AEL(2)+=DBCU(2);
        AEL(2)/=up.VOL.Read(i3);

        if (AEL(0)>(nprec)0.0)
        {
            AC(i,0)=cpe*AEL(0);
            massNp(i1)+=cpe*AEL(0);
        }else
        {
            AC(i,1)=cme*AEL(0);
            massNn(i1)-=cme*AEL(0);
        };
        if (AEL(1)>(nprec)0.0)
        {
            AC(i,2)=cpe*AEL(1);
            massNp(i2)+=cpe*AEL(1);
        }else
        {
            AC(i,3)=cme*AEL(1);
            massNn(i2)-=cme*AEL(1);
        };
        if (AEL(2)>(nprec)0.0)
        {
            AC(i,4)=cpe*AEL(2);
            massNp(i3)+=cpe*AEL(2);
        }else
        {
            AC(i,5)=cme*AEL(2);
            massNn(i3)-=cme*AEL(2);
        };
    };

    for (unsigned int i=0; i<up.LOS.GetSize(); i++)
    {
        massn+=massNn(i)*up.VOL.Read(i);
        massp+=massNp(i)*up.VOL.Read(i);
    };

    nprec cpe2,cme2;
    if (std::abs(massp-massn)>(nprec)0.0)
    {
        if (massp>massn)
        {
            cme2=(nprec)1.0;
            cpe2=massn/massp;
        }else
        {
            cme2=massp/massn;
            cpe2=(nprec)1.0;
        };

    }else
    {
        cpe2=(nprec)1.0;
        cme2=(nprec)1.0;
    };

    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        i1=geo.elem.ReadNode(i,0); i2=geo.elem.ReadNode(i,1); i3=geo.elem.ReadNode(i,2);
        ACf(i1)+=AC(i,0)*cpe2+AC(i,1)*cme2;
        ACf(i2)+=AC(i,2)*cpe2+AC(i,3)*cme2;
        ACf(i3)+=AC(i,4)*cpe2+AC(i,5)*cme2;
    };


    if (prelimiter)
    {
        for (unsigned int i=0; i<up.LOS.GetSize(); i++)
        {
            if (up.LOS(i)>ULE(i)&&ACf(i)<(nprec)0.0)
            {
                ACf(i)=(nprec)0.0;
            };
            if (up.LOS(i)<ULE(i)&&ACf(i)>(nprec)0.0)
            {
                ACf(i)=(nprec)0.0;
            };
        };
    };

    up.HOS=up.LOS+ACf;
    /*{
    nprec masscheck=0.0;
    nprec mass1=0.0;
    nprec mass2=0.0;
    for (unsigned int i=0; i<geo.nd.GetNnodes(); i++)
    {
    	mass1+=up.LOS.Read(i)*up.VOL.Read(i);
    	mass2+=up.HOS.Read(i)*up.VOL.Read(i);
    };
    masscheck=(mass2-mass1)/mass1;
    std::cout << "Mass Check FCT: " << masscheck << std::endl;
    }*/
}
