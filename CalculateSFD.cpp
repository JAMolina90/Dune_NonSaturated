#include <cmath>
#include "Globals.h"
#include "matrix.h"
#include "myvector.h"
#include "geometry.h"

void CalculateSFD(matrix<nprec>& SFD, const geometry& geo)
{
    SFD.SetDimensions(geo.elem.GetNelem(),10);
    Vector<nprec> X(3);
    Vector<nprec> Y(3);
    nprec twoarea;
    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        X(0)=geo.nd.ReadC(geo.elem.ReadNode(i,0),0);
        X(1)=geo.nd.ReadC(geo.elem.ReadNode(i,1),0);
        X(2)=geo.nd.ReadC(geo.elem.ReadNode(i,2),0);
        Y(0)=geo.nd.ReadC(geo.elem.ReadNode(i,0),1);
        Y(1)=geo.nd.ReadC(geo.elem.ReadNode(i,1),1);
        Y(2)=geo.nd.ReadC(geo.elem.ReadNode(i,2),1);
        twoarea=(X(1)-X(0))*(Y(2)-Y(0))-(X(2)-X(0))*(Y(1)-Y(0));
        SFD(i,0)=(Y(1)-Y(2))/twoarea;
        SFD(i,1)=(Y(2)-Y(0))/twoarea;
        SFD(i,2)=(Y(0)-Y(1))/twoarea;
        SFD(i,3)=(X(2)-X(1))/twoarea;
        SFD(i,4)=(X(0)-X(2))/twoarea;
        SFD(i,5)=(X(1)-X(0))/twoarea;
        SFD(i,6)=twoarea;
        SFD(i,7)=(X(1)*Y(2)-X(2)*Y(1))/twoarea;
        SFD(i,8)=(X(2)*Y(0)-X(0)*Y(2))/twoarea;
        SFD(i,9)=(X(0)*Y(1)-X(1)*Y(0))/twoarea;
    };
}
