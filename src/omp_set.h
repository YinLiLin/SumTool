#ifndef OMP_SET_H_
#define OMP_SET_H_

#if defined(_OPENMP)
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#else
#endif

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

static int omp_setup(int threads);
static inline int omp_setup(int threads=0) {
    int t = 1;
#ifdef _OPENMP
    if (threads == 0) {
        t = omp_get_num_procs() - 1;
        t = t > 0 ? t : 1;
    } else {
        t = threads > 0 ? threads : 1;
    }
    omp_set_num_threads(t);
#else
#endif
    return t;
}

#endif
