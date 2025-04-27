#ifndef GEOMETRY_H_INCLUDED
#define GEOMETRY_H_INCLUDED

#include "Globals.h"
#include "bound.h"
#include "node.h"
#include "triangle.h"
#include "myvector.h"
#include "matrix.h"

struct geometry {
    triangle elem;
    node nd;
    bound bd;
    matrix<int> wnodes;
    matrix<int> slipnodes;
    matrix<nprec> slipnormal;
    matrix<int> opennodes;
    matrix<nprec> opennormal;
    Vector<unsigned int> nodesu0;
};


#endif // GEOMETRY_H_INCLUDED
