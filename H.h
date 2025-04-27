#ifndef H_H_INCLUDED
#define H_H_INCLUDED

void HvaluesIni(AuxVarGeo& s, const geometry& geo);

void HnodesCalculationIni(AuxVarGeo& s, const geometry& geo);

void UastCalculation(AuxVarGeo& s, const geometry& geo, const matrix<nprec>& unknowns1, parameter& par, const matrix<nprec>& SFD, Neig& NB);

void Sediment(AuxVar& av, AuxVarGeo& s, const geometry& geo, const matrix<nprec>& unknowns1, parameter& par, const matrix<nprec>& SFD);

void HnodesCalculation(AuxVarGeo& s, const geometry& geo, const matrix<nprec>& unknowns1, parameter& par, const matrix<nprec>& SFD, nprec timepar, nprec currenttime);

void NewMesh(const geometry& geo0, geometry& geo1, parameter& par, AuxVarGeo& avG, const matrix<nprec>& SFD0, matrix<nprec>& SFD1, UpwindVar& upvel0, UpwindVar& upvel1, matrix<nprec>& RHSUW1, matrix<nprec>& unknowns0, matrix<nprec>& unknowns1);

#endif // H_H_INCLUDED
