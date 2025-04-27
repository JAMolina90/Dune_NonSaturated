#ifndef BOUND_H
#define BOUND_H

#include "Globals.h"

class bound
{
    public:
        bound();
        bound(unsigned int n);
        ~bound();
        void SetNbounds(unsigned int n);
        unsigned int GetNbounds() const;
        void SetBC(unsigned int bound,unsigned int boundaryc);
        unsigned int ReadBC(unsigned int bound) const;
        void SetNode(unsigned int bound,unsigned int node, unsigned int N);
        unsigned int ReadNode(unsigned int bound,unsigned int node) const;
        void SetEO(unsigned int bound,unsigned int e);
        unsigned int ReadEO(unsigned int bound) const;
        void SetNormal(unsigned int bound, unsigned int coor, nprec N);
        nprec ReadNormal(unsigned int bound, unsigned int coor) const;
        void SetLen(unsigned int bound, nprec l);
        nprec ReadLen(unsigned int bound) const;

    private:
        unsigned int nbound;
        unsigned int* nodes;
        unsigned int* EO;
        unsigned int* BC;
        nprec* normals;
        nprec* lengths;
};

#endif // BOUND_H
