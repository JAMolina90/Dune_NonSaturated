#ifndef FIRSTSTEP_H
#define FIRSTSTEP_H

void ConvTerm(const geometry& geo, const parameter& par, AuxVar& av, const matrix<nprec>& SFD, const matrix<nprec>& unknowns, AuxVarGeo& avG);
void FirstStep(const geometry& geo, const parameter& par, AuxVar& av, const matrix<nprec>& SFD, const matrix<nprec>& unknowns, AuxVarGeo& avG, const int& ind);

#endif // FIRSTSTEP_H
