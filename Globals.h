#ifndef GLOBALS_H_INCLUDED
#define GLOBALS_H_INCLUDED

typedef double nprec;

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_num_threads() 1
    #define omp_get_max_threads() 1
#endif

#endif // GLOBALS_H_INCLUDED
