#ifndef FIRSTSTEPUPWIND_H
#define FIRSTSTEPUPWIND_H

void FirstStepUpwind(const geometry& geo, const parameter& par, AuxVar& av, const matrix<nprec>& unknowns, matrix<nprec>& RHSUW, UpwindVar& upvel);

#endif // FIRSTSTEPUPWIND_H
