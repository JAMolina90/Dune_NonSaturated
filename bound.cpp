#include "bound.h"
#include <cassert>

bound::bound()
{
    nbound = 0;
    nodes = new unsigned int [nbound*2];
    EO = new unsigned int [nbound];
    BC = new unsigned int [nbound];
    normals = new nprec [nbound*2];
    lengths = new nprec [nbound];
}

bound::bound(unsigned int n)
{
    nbound = n;
    nodes = new unsigned int [nbound*2];
    EO = new unsigned int [nbound];
    BC = new unsigned int [nbound];
    normals = new nprec [nbound*2];
    lengths = new nprec [nbound];
}

bound::~bound()
{
    delete[] nodes;
    delete[] EO;
    delete[] BC;
    delete[] normals;
    delete[] lengths;
}

void bound::SetNbounds(unsigned int n)
{
    nbound = n;
    delete[] nodes;
    delete[] EO;
    delete[] BC;
    delete[] normals;
    delete[] lengths;
    nodes = new unsigned int [nbound*2];
    EO = new unsigned int [nbound];
    BC = new unsigned int [nbound];
    normals = new nprec [nbound*2];
    lengths = new nprec [nbound];
}

unsigned int bound::GetNbounds() const
{
    return nbound;
}

void bound::SetBC(unsigned int bound,unsigned int boundaryc)
{
    assert(bound<nbound);
    BC[bound]=boundaryc;
}

unsigned int bound::ReadBC(unsigned int bound) const
{
    assert(bound<nbound);
    return BC[bound];
}

void bound::SetNode(unsigned int bound,unsigned int node, unsigned int N)
{
    assert(bound<nbound);
    assert(node<2);
    nodes[bound*2+node]=N;
}

unsigned int bound::ReadNode(unsigned int bound,unsigned int node) const
{
    assert(bound<nbound);
    assert(node<2);
    return nodes[bound*2+node];
}

void bound::SetEO(unsigned int bound,unsigned int e)
{
    assert(bound<nbound);
    EO[bound]=e;
}

unsigned int bound::ReadEO(unsigned int bound) const
{
    assert(bound<nbound);
    return EO[bound];
}

void bound::SetNormal(unsigned int bound, unsigned int coor, nprec N)
{
    assert(bound<nbound);
    assert(coor<2);
    normals[bound*2+coor]=N;
}

nprec bound::ReadNormal(unsigned int bound, unsigned int coor) const
{
    assert(bound<nbound);
    assert(coor<2);
    return normals[bound*2+coor];
}

void bound::SetLen(unsigned int bound, nprec l)
{
    assert(bound<nbound);
    lengths[bound]=l;
}

nprec bound::ReadLen(unsigned int bound) const
{
    assert(bound<nbound);
    return lengths[bound];
}
