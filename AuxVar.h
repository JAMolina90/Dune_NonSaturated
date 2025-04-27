#ifndef AuxVar_H_INCLUDED
#define AuxVar_H_INCLUDED

#include <stdlib.h>
#include <cmath>
#include <limits>
#include "Globals.h"
#include "geometry.h"
#include "myvector.h"
#include "matrix.h"
#include "parameter.h"
#include "UpwindVar.h"
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/SparseCholesky"

struct AuxVar {
    Vector<nprec> TT;
    matrix<nprec> RHS;
    matrix<nprec> ConvTerm;
    Vector<nprec> FRHS;
    Vector<nprec> u2a;
    Vector<nprec> v2a;
    Vector<nprec> TwoWayForce;
    //matrix<nprec> SSM;
    Vector<nprec> DMAT;
    Vector<nprec> rhon;
    Vector<nprec> dviscon;
    Vector<nprec> rhon1;
    Vector<nprec> dviscon1;
    Vector<nprec> epsilon;
    Vector<nprec> rhomel;
    //matrix<nprec> LHSM;
    matrix<nprec> VelPred;
    //Vector<bool> dalphadn0;
    Vector<bool> dudn0;
    Vector<bool> dpdn0;
    Vector<bool> dtaudn0;
    Vector<bool> slipelem;
    Vector<bool> reinitact;
    matrix<nprec> normalel;
    nprec dtm;
    nprec tol;
    nprec totalarea;
    bool inletoutlet;
    Vector<nprec> walllawY;
    Vector<bool> walllawind;
    Vector<unsigned int>walllawnode;    //here we storage the node which is not at boundary
    int ipfix;
    Eigen::SparseMatrix<nprec> SSM2;
    Eigen::SparseMatrix<nprec> LHSM2;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<nprec>> solverCMM,solverPress;
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<nprec>, Eigen::Lower|Eigen::Upper> solverCMM,solverPress;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<nprec>,Eigen::IncompleteLUT<nprec>> solverCMM,solverPress;

};

struct AuxVarGeo {

    matrix<int> Edges; //Edge number and nodes of each edge
    matrix<int> SharedEdges; //for each edge, we store elements which contains that edge
    matrix<int> ElemEdges; //For each element, we store each edge number
    Vector<int> BoundEdges; //Each boundary is a Edge
    Vector<unsigned int> EdgeCase;
    //matrix<unsigned int> geoantes; //Position0: interface element yes or no, Position1-3:second to fourth elements position resulting of dividing into 4 element, Position4-6: position of new nodes
    //matrix<unsigned int> geodespues;
    //matrix<unsigned int> EdgeRefinedantes; //First Position --> Yes or Not; Second postion --> New node number
    //matrix<unsigned int> EdgeRefineddespues;


    Vector<unsigned int> nodesh;
    Vector<nprec> h;
    Vector<nprec> nodeshX;
    nprec Hmax;
    matrix<unsigned int> linesh; //Pos0: Node 0; Pos1; Node 1;
    Vector<nprec> lineslen;
    Vector<int> linklist; //Position: real node name; Value: virtual node name
    matrix<unsigned int> conectivity; //Each position is one node, column 0 is number of elements connected with node, column 1,2 are the elements connected with node
    Vector<unsigned int> nodesbd;
    Vector<nprec> nodeshqs;
    Vector<nprec> nodeshMeanqs;
    Vector<nprec> nodeshveff;
    Vector<nprec> nodeshsedrho;
    Vector<nprec> nodeshuast;

    Vector<int> NodeNumberEdge;
    Vector<bool> PrevEdgecutted;
    Vector<nprec> PrevEdgeValue;

    Vector<unsigned int> ActivatedNodes; //0--> deactivated, 1 activated, 2 activated and bed passing through it
    Vector<unsigned int> ActivatedElem;  //0--> deactivated, 1 activated, 2 divided
    Vector<int> OldNodesIndex;
    Vector<int> OldElemIndex;
    Vector<int> NeigCorresp;
};

inline void TimeStep(const geometry& geo, const matrix<nprec>& SFD, const matrix<nprec>& unknowns, parameter& par, nprec timepar, nprec currenttime)
{
    unsigned int n1,n2,n3;
    nprec maxtime=par.delt*par.Uref/par.Lref;
    nprec len1,len2,len3,len,um,vm,ut;
    nprec one=(nprec)1.0;
    nprec maxCourant=std::numeric_limits<nprec>::min();
    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        n1=geo.elem.ReadNode(i,0);
        n2=geo.elem.ReadNode(i,1);
        n3=geo.elem.ReadNode(i,2);
        len1=one/(std::sqrt(SFD.Read(i,0)*SFD.Read(i,0)+SFD.Read(i,3)*SFD.Read(i,3)));
        len2=one/(std::sqrt(SFD.Read(i,1)*SFD.Read(i,1)+SFD.Read(i,4)*SFD.Read(i,4)));
        len3=one/(std::sqrt(SFD.Read(i,2)*SFD.Read(i,2)+SFD.Read(i,5)*SFD.Read(i,5)));
        len=std::min(std::min(len1,len2),len3);
        um=(unknowns.Read(n1,1)+unknowns.Read(n2,1)+unknowns.Read(n3,1))/((nprec)3.0);
        vm=(unknowns.Read(n1,2)+unknowns.Read(n2,2)+unknowns.Read(n3,2))/((nprec)3.0);
        ut=std::sqrt(um*um+vm*vm);
        if (ut>1e-9){maxtime=std::min(maxtime,len/ut*par.maxCou);};
        maxCourant=std::max(maxCourant,ut*par.DELTP/len);
    };

    if (timepar+maxtime>par.interplot){maxtime=par.interplot-timepar;};
    if (currenttime+maxtime>par.tfinal){maxtime=par.tfinal-currenttime;};

    par.DELTP=maxtime;
    std::cout << "Max Courant Number: " << maxCourant << std::endl;
}


inline void AuxVarGeoIni(AuxVarGeo& avG, const geometry& geo0)
{
    avG.ActivatedNodes.SetDimensions(geo0.nd.GetNnodes());
    avG.ActivatedNodes.Initialize(0);
    avG.ActivatedElem.SetDimensions(geo0.elem.GetNelem());
    avG.ActivatedElem.Initialize(0);
    avG.OldNodesIndex.SetDimensions(geo0.nd.GetNnodes());
    avG.OldElemIndex.SetDimensions(geo0.elem.GetNelem());
};

inline void AuxVarIni(AuxVar& s, const geometry& geo, const matrix<nprec>& SFD, const parameter& par)
{
    unsigned int nnodes = geo.nd.GetNnodes();
    unsigned int nelem = geo.elem.GetNelem();
    nprec one=(nprec)1.0;
    //s.P0.SetDimensions(nnodes); //Pressure values
    //s.P0.Initialize((nprec)0.0);
    s.TT.SetDimensions(nnodes); //Vector with 1.0 for free nodes and 0.0 for prescribed p nodes
    s.TT.Initialize((nprec)1.0);
    s.RHS.SetDimensions(nnodes,4);
    s.FRHS.SetDimensions(nnodes);
    //s.P1.SetDimensions(nnodes);
    s.u2a.SetDimensions(nnodes);
    s.v2a.SetDimensions(nnodes);
    s.TwoWayForce.SetDimensions(geo.nd.GetNnodes());
    s.TwoWayForce.Initialize((nprec)0.0);
    s.DMAT.SetDimensions(nnodes);
    s.DMAT.Initialize((nprec)0.0);
    s.rhon.SetDimensions(nnodes);
    s.dviscon.SetDimensions(nnodes);
    s.rhon1.SetDimensions(nnodes);
    s.dviscon1.SetDimensions(nnodes);
    //s.SSM.SetDimensions(nelem,6);
    //s.LHSM.SetDimensions(nelem,6);
    s.VelPred.SetDimensions(nnodes,2);
    s.VelPred.Initialize((nprec)0.0);
    s.epsilon.SetDimensions(nelem);
    s.rhomel.SetDimensions(nelem);
    //s.dalphadn0.SetDimensions(nelem);
    //s.dalphadn0.Initialize(false);
    s.dudn0.SetDimensions(nelem);
    s.dudn0.Initialize(false);
    s.dpdn0.SetDimensions(nelem);
    s.dpdn0.Initialize(false);
    s.dtaudn0.SetDimensions(nelem);
    s.dtaudn0.Initialize(false);
    s.slipelem.SetDimensions(nelem);
    s.slipelem.Initialize(false);
    s.inletoutlet=false;
    s.normalel.SetDimensions(nelem,2);
    s.normalel.Initialize((nprec)0.0);
    s.reinitact.SetDimensions(nelem);
    s.reinitact.Initialize(true);
    s.walllawind.SetDimensions(nelem);
    s.walllawind.Initialize(false);
    s.walllawY.SetDimensions(nelem);
    s.walllawY.Initialize(0.0);
    s.walllawnode.SetDimensions(nelem);
    s.SSM2.resize(nnodes,nnodes);
    s.LHSM2.resize(nnodes,nnodes);
    s.ConvTerm.SetDimensions(nnodes,2);

    nprec area;
    s.totalarea=(nprec)0.0;
    //unsigned int minelem=0;
    //nprec minarea=std::numeric_limits<nprec>::max();

    for (unsigned int i=0; i<nelem; i++)
    {
        area = SFD.Read(i,6)/((nprec)2.0);
        s.totalarea+=area;
        s.DMAT(geo.elem.ReadNode(i,0))+=area/((nprec)3.0);
        s.DMAT(geo.elem.ReadNode(i,1))+=area/((nprec)3.0);
        s.DMAT(geo.elem.ReadNode(i,2))+=area/((nprec)3.0);
    };

    for (unsigned int i=0; i<s.DMAT.GetSize(); i++)
    {
        s.DMAT(i)=((nprec)1.0)/s.DMAT(i);
    };

    //s.LHSM.Initialize((nprec)0.0);
    nprec two=(nprec)2.0;
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(geo.elem.GetNelem()*9);
    unsigned int n1,n2,n3;
    for (unsigned int i=0; i<geo.elem.GetNelem(); i++)
    {
        n1=geo.elem.ReadNode(i,0);
        n2=geo.elem.ReadNode(i,1);
        n3=geo.elem.ReadNode(i,2);
        area=SFD.Read(i,6)/((nprec)24.0);
        /*
        s.LHSM(i,0)=two*area;
        s.LHSM(i,1)=one*area;
        s.LHSM(i,2)=two*area;
        s.LHSM(i,3)=one*area;
        s.LHSM(i,4)=one*area;
        s.LHSM(i,5)=two*area;
        */
        tripletList.push_back(T(n1,n1,two*area));
        tripletList.push_back(T(n2,n2,two*area));
        tripletList.push_back(T(n3,n3,two*area));
        tripletList.push_back(T(n1,n2,one*area));
        tripletList.push_back(T(n2,n1,one*area));
        tripletList.push_back(T(n1,n3,one*area));
        tripletList.push_back(T(n3,n1,one*area));
        tripletList.push_back(T(n2,n3,one*area));
        tripletList.push_back(T(n3,n2,one*area));
    };
    s.LHSM2.setFromTriplets(tripletList.begin(), tripletList.end());
    s.solverCMM.compute(s.LHSM2);

    s.rhon.Initialize(par.rho0/par.rhoref);
    s.rhon1.Initialize(par.rho0/par.rhoref);
    s.dviscon.Initialize(par.dvisco0/par.dviscoref);
    s.dviscon1.Initialize(par.dvisco0/par.dviscoref);
    s.rhomel.Initialize(par.rho0/par.rhoref);

};

inline void PressureInterpolation(const geometry& geoRef, const geometry& geo1, matrix<nprec>& unknownsRef, const matrix<nprec>& unknowns1, AuxVarGeo& avG)
{
    /*
    Vector<unsigned int> listaN(geo1.nd.GetNnodes());
    listaN.Initialize(0);
    unsigned int contN=0;
    {
    Vector<bool> check(geo1.nd.GetNnodes());
    check.Initialize(false);
    for (unsigned int i=0; i<geo1.bd.GetNbounds(); i++)
    {
        if (geo1.bd.ReadBC(i)==20)
        {
            if (check(geo1.bd.ReadNode(i,0))==false)
            {
                check(geo1.bd.ReadNode(i,0))=true;
                listaN(contN)=geo1.bd.ReadNode(i,0);
                contN++;
            };
            if (check(geo1.bd.ReadNode(i,1))==false)
            {
                check(geo1.bd.ReadNode(i,1))=true;
                listaN(contN)=geo1.bd.ReadNode(i,1);
                contN++;
            };
        };
    }
    }

    unsigned int i,j,nearestnode;
    nprec distmin,dist,X,Y;
    #pragma omp parallel for default(none) shared(geoRef,geo1,avG,unknowns1,unknownsRef,listaN,contN) private(i,j,nearestnode,dist,distmin,X,Y)
    for (i=0; i<geoRef.nd.GetNnodes(); i++)
    {
        if (avG.ActivatedNodes.Read(i)==0)
        {
            X=geoRef.nd.ReadC(i,0);
            Y=geoRef.nd.ReadC(i,1);
            nearestnode=0;
            distmin=std::numeric_limits<nprec>::max();
            for (j=0; j<contN; j++)
            {
                dist=std::sqrt((X-geo1.nd.ReadC(listaN(j),0))*(X-geo1.nd.ReadC(listaN(j),0))+(Y-geo1.nd.ReadC(listaN(j),1))*(Y-geo1.nd.ReadC(listaN(j),1)));
                if (dist<distmin)
                {
                    distmin=dist;
                    nearestnode=j;
                };
            };
            unknownsRef(i,0)=unknowns1.Read(listaN(nearestnode),0);
        };
    };
    */
    Vector<unsigned int> Cont(geoRef.nd.GetNnodes());
    Cont.Initialize(0);
    Vector<nprec> Value(geoRef.nd.GetNnodes());
    Value.Initialize((nprec)0.0);
    for (unsigned int i=0; i<avG.Edges.GetRows(); i++)
    {
        if (avG.EdgeCase.Read(i)==2)
        {
            if (avG.ActivatedNodes.Read((unsigned int)avG.Edges.Read(i,0))==0)
            {
                Value((unsigned int)avG.Edges.Read(i,0))+=unknowns1.Read((unsigned int)avG.NodeNumberEdge.Read(i),0);
                Cont((unsigned int)avG.Edges.Read(i,0))+=1;
            }
            else if (avG.ActivatedNodes.Read((unsigned int)avG.Edges.Read(i,1))==0)
            {
                Value((unsigned int)avG.Edges.Read(i,1))+=unknowns1.Read((unsigned int)avG.NodeNumberEdge.Read(i),0);
                Cont((unsigned int)avG.Edges.Read(i,1))+=1;
            }
            else
            {
                std::cout << "Something is wrong" << std::endl;
                exit(1);
            };
        };
    };

    for (unsigned int i=0; i<geoRef.nd.GetNnodes(); i++)
    {
        if (Cont(i)>0)
        {
            unknownsRef(i,0)=Value(i)/((nprec)Cont(i));
        };
    };
};

#endif // AuxVar_H_INCLUDED
