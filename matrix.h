#ifndef MATRIX2_H
#define MATRIX2_H

#include <cassert>
#include <cmath>
#include <string>
#include "Globals.h"
#include "myvector.h"

template<class T>
class matrix
{
    public:
        matrix();
        matrix(unsigned int row, unsigned int column);                         //Constructor
        matrix(const matrix<T>& otherMatrix);              //Copy
        ~matrix();                              //Destructor
        T& operator[] (unsigned int i);
        T& operator() (unsigned int i, unsigned int j);                        //Reader
        matrix<T> operator () (unsigned int i1, unsigned int i2, unsigned int j1, unsigned int j2);
        matrix<T>& operator=(const matrix<T>& otherMatrix);   //operator =
        matrix<T> operator+() const;
        matrix<T> operator-() const;
        matrix<T> operator+(const matrix<T>& m1) const;
        matrix<T>& operator+=(const matrix<T>& m1);
        matrix<T> operator-(const matrix<T>& m1) const;
        matrix<T>& operator-=(const matrix<T>& m1);
        matrix<T> operator*(T& a) const;
        matrix<T> operator*(const matrix<T>& m2) const;
        T Read(unsigned int i, unsigned int j) const;
        //T Read(unsigned int i) const;
        unsigned int GetRows() const;
        unsigned int GetColumns() const;
        void Initialize(T number);
        void SetDimensions(unsigned int row, unsigned int column);

        //Template friend operators have to be declared inline
        friend Vector<T> operator*(const matrix<T>& m, const Vector<T>& v)
        {
            assert(m.C==v.vSize);
            Vector<T> res(m.C);
            res.Initialize((T)0);
            for(unsigned int i=0; i<m.R; i++)
            {
                for (unsigned int j=0; j<m.C; j++)
                {
                    res.vData[i]+=m.mData[i*m.C+j]*v.vData[j];
                };
            };
            return res;
        };

        friend Vector<T> operator*(const Vector<T>& v,const matrix<T>& m)
        {
            assert(m.R==v.vSize);
            Vector<T> res(m.C);
            res.Initialize((T)0);
            //#pragma omp parallel for shared(res,R,C,v,m) private(i,j)
            for(unsigned int i=0; i<m.C; i++)
            {
                for (unsigned int j=0; j<m.R; j++)
                {
                    res.vData[i]+=v.vData[j]*m.mData[j*m.C+i];
                };
            };
            return res;
        };

        friend Vector<T> vectorize(const matrix<T>& m, unsigned int rc, std::string rowcolum)
        {
            assert(rowcolum=="row"||rowcolum=="column");
            Vector<T> res;
            if (rowcolum=="row")
            {
                assert(rc<m.R);
                res.SetDimensions(m.C);
                //#pragma omp parallel for shared(res,C,m,rc) private (i)
                for (unsigned int i=0; i<m.C; i++)
                {
                    res(i)=m.mData[rc*m.C+i];
                };
            }
            else
            {
                assert(rc<m.GetColumns());
                res.SetDimensions(m.GetRows());
                //#pragma omp parallel for shared(res,R,m,rc) private (i)
                for (unsigned int i=0; i<m.R; i++)
                {
                    res(i)=m.mData[i*m.C+rc];
                };
            };
            return res;
        };

    private:
        T* mData;                           //Pointer to data
        unsigned int R;
        unsigned int C;
        unsigned int D;
};

//template<class T> Vector<T> operator*(const matrix<T>& m, const Vector<T>& v);
//template<class T> Vector<T> operator*(const Vector<T>& v,const matrix<T>& m);

template <class T>
matrix<T>::matrix()
{
    R = 0;
    C = 0;
    D = R*C;
    mData = new T[D];

};

template <class T>
matrix<T>::matrix(unsigned int row, unsigned int column)
{
    R = row;
    C = column;
    D = R*C;
    mData = new T[D];

};

template <class T>
matrix<T>::matrix(const matrix<T>& otherMatrix)
{
    R = otherMatrix.GetRows();
    C = otherMatrix.GetColumns();
    D = otherMatrix.D;
    mData = new T[D];
    //#pragma omp parallel for shared(otherMatrix) private(i)
    for (unsigned int i=0; i<D; i++)
    {
        mData[i] = otherMatrix.mData[i];
    };

};

template <class T>
matrix<T>::~matrix()
{
    delete[] mData;
};

template <class T>
unsigned int matrix<T>::GetRows() const
{
    return R;
};

template <class T>
unsigned int matrix<T>::GetColumns() const
{
    return C;
};

template <class T>
T& matrix<T>::operator[] (unsigned int i)
{
    assert(i < D);
    return mData[i];
};

template <class T>
T& matrix<T>::operator() (unsigned int i, unsigned int j)
{
    assert(i < R);
    assert(j < C);
    return mData[i*C+j];
};

template <class T>
matrix<T> matrix<T>::operator () (unsigned int i1, unsigned int i2, unsigned int j1, unsigned int j2)
{
    assert(i1<R);
    assert(i2<R);
    assert(j1<C);
    assert(j2<C);
    assert(i1<=i2);
    assert(j1<=j2);
    matrix<T> res(i2-i1+1,j2-j1+1);
    unsigned row=0;
    unsigned col=0;
    for (unsigned int i=i1; i<i2+1; i++)
    {
        for (unsigned int j=j1; j<j2+1; j++)
        {
            res(row,col)=mData[i*C+j];
            col++;
        };
        row++;
    };
    return res;
};

template <class T>
matrix<T>& matrix<T>::operator=(const matrix<T>& otherMatrix)
{
    assert(D == otherMatrix.D);
    unsigned int i;
    //#pragma omp parallel for shared(otherMatrix) private(i)
    for (i=0; i<D; i++)
    {
        mData[i] = otherMatrix.mData[i];
    };
    return *this;
};

template <class T>
void matrix<T>::Initialize(T number)
{
    unsigned int i;
    //#pragma omp parallel for shared(number) private(i)
    for (i=0; i<D; i++)
    {
        mData[i] = number;
    };
};

template <class T>
void matrix<T>::SetDimensions(unsigned int row, unsigned int column)
{
    R=row;
    C=column;
    D=R*C;
    delete[] mData;
    mData = new T[D];
};

template <class T>
matrix<T> matrix<T>::operator+() const
{
    matrix<T> mat(R,C);
    unsigned int i;
    //#pragma omp parallel for shared(mat) private(i)
    for (i=0; i<D; i++)
    {
        mat.mData[i] = mData[i];
    };

    return mat;
};

template <class T>
matrix<T> matrix<T>::operator-() const
{
    matrix<T> mat(R,C);
    unsigned int i;
    //#pragma omp parallel for shared(mat) private(i)
    for (i=0; i<D; i++)
    {
        mat.mData[i] = -mData[i];
    };

    return mat;
};

template <class T>
matrix<T> matrix<T>::operator+(const matrix<T>& m1) const
{
    assert(R==m1.GetRows());
    assert(C==m1.GetColumns());
    matrix<T> res(R,C);
    unsigned int i;
    //#pragma omp parallel for shared(res,m1) private(i)
    for (i=0; i<D; i++)
    {
        res.mData[i] = mData[i]+m1.mData[i];
    };
    return res;
};

template <class T>
matrix<T> matrix<T>::operator-(const matrix<T>& m1) const
{
    assert(R==m1.GetRows());
    assert(C==m1.GetColumns());
    matrix<T> res(R,C);
    unsigned int i;
    //#pragma omp parallel for shared(res,m1) private(i)
    for (i=0; i<D; i++)
    {
        res.mData[i] = mData[i]-m1.mData[i];
    };
    return res;
};

template <class T>
matrix<T> matrix<T>::operator*(T& a) const
{
    matrix<T> res(R,C);
    unsigned int i;
    //#pragma omp parallel for shared(res,a) private(i)
    for (i=0; i<D; i++)
    {
        res.mData[i] = mData[i]*a;
    };
    return res;
};

template <class T>
T matrix<T>::Read(unsigned int i, unsigned int j) const
{
    assert(i<R);
    assert(j<C);
    return mData[i*C+j];
}

/*
template <class T>
T matrix<T>::Read(unsigned int i) const
{
    assert(i<D);
    return mData[i];
}
*/

template <class T>
matrix<T> matrix<T>::operator*(const matrix<T>& m2) const
{
    assert(C==m2.R);
    matrix<T> res(R,m2.C);
    res.Initialize((T)0);
    unsigned int i,j,k;
    unsigned int m2c=m2.C;
    //#pragma omp parallel for shared(m2,res,m2c) private(i,j,k)
    for (i=0; i<R; i++)
    {
        for (j=0; j<m2c; j++)
        {
            for (k=0; k<C; k++)
            {
                res.mData[i*m2c+j]+=mData[i*C+k]*m2.mData[k*m2c+j];
            };
        };
    };
    return res;
};

template <class T>
matrix<T>& matrix<T>::operator+=(const matrix<T>& m1)
{
    assert(D==m1.D);
    for (unsigned int i=0; i<D; i++)
    {
        this->mData[i]+=m1.mData[i];
    };
    return *this;
}

template <class T>
matrix<T>& matrix<T>::operator-=(const matrix<T>& m1)
{
    assert(D==m1.D);
    for (unsigned int i=0; i<D; i++)
    {
        this->mData[i]-=m1.mData[i];
    };
    return *this;
}

#endif // MATRIX2_H
