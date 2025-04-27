#ifndef FCTU_H_INCLUDED
#define FCTU_H_INCLUDED

void FCTU(const parameter& par, const geometry& geo, const matrix<nprec>& SFD, const Vector<nprec>& DMAT, const Vector<nprec>& varPre1, const Vector<nprec>& varPre2, const Vector<nprec>& varL1, const Vector<nprec>& varL2, const Vector<nprec>& varH1, const Vector<nprec>& varH2, matrix<nprec>& FS, const Vector<nprec>& Phi);
void FCTU2(const geometry& geo, const matrix<nprec>& SFD, const Vector<nprec>& DMAT, const Vector<nprec>& varPre1, const Vector<nprec>& varPre2, const Vector<nprec>& varL1, const Vector<nprec>& varL2, const Vector<nprec>& varH1, const Vector<nprec>& varH2, matrix<nprec>& FS, const Vector<nprec>& dUho, const Vector<nprec>& dVho);

#endif // FCTU_H_INCLUDED
