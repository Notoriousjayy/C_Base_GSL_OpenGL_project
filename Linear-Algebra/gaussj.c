//
// Created by jorda on 7/15/2024.
//
#include "gaussj.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define SWAP macro for variables
#define SWAP(a, b) { double temp = (a); (a) = (b); (b) = temp; }

void gaussj(gsl_matrix* a, gsl_matrix* b) {
    size_t n = a->size1;
    size_t m = b->size2;
    gsl_vector* indxc = gsl_vector_alloc(n);
    gsl_vector* indxr = gsl_vector_alloc(n);
    gsl_vector* ipiv = gsl_vector_alloc(n);

    for (size_t j = 0; j < n; ++j) gsl_vector_set(ipiv, j, 0);

    for (size_t i = 0; i < n; ++i) {
        double big = 0.0;
        size_t irow = 0, icol = 0;
        for (size_t j = 0; j < n; ++j) {
            if (gsl_vector_get(ipiv, j) != 1) {
                for (size_t k = 0; k < n; ++k) {
                    if (gsl_vector_get(ipiv, k) == 0) {
                        double a_jk = gsl_matrix_get(a, j, k);
                        if (fabs(a_jk) >= big) {
                            big = fabs(a_jk);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }

        gsl_vector_set(ipiv, icol, gsl_vector_get(ipiv, icol) + 1);

        if (irow != icol) {
            for (size_t l = 0; l < n; ++l) {
                double temp = gsl_matrix_get(a, irow, l);
                gsl_matrix_set(a, irow, l, gsl_matrix_get(a, icol, l));
                gsl_matrix_set(a, icol, l, temp);
            }
            for (size_t l = 0; l < m; ++l) {
                double temp = gsl_matrix_get(b, irow, l);
                gsl_matrix_set(b, irow, l, gsl_matrix_get(b, icol, l));
                gsl_matrix_set(b, icol, l, temp);
            }
        }

        gsl_vector_set(indxr, i, irow);
        gsl_vector_set(indxc, i, icol);

        if (gsl_matrix_get(a, icol, icol) == 0.0) {
            fprintf(stderr, "gaussj: Singular Matrix\n");
            exit(EXIT_FAILURE);
        }

        double pivinv = 1.0 / gsl_matrix_get(a, icol, icol);
        gsl_matrix_set(a, icol, icol, 1.0);

        for (size_t l = 0; l < n; ++l) gsl_matrix_set(a, icol, l, gsl_matrix_get(a, icol, l) * pivinv);
        for (size_t l = 0; l < m; ++l) gsl_matrix_set(b, icol, l, gsl_matrix_get(b, icol, l) * pivinv);

        for (size_t ll = 0; ll < n; ++ll) {
            if (ll != icol) {
                double dum = gsl_matrix_get(a, ll, icol);
                gsl_matrix_set(a, ll, icol, 0.0);
                for (size_t l = 0; l < n; ++l) gsl_matrix_set(a, ll, l, gsl_matrix_get(a, ll, l) - gsl_matrix_get(a, icol, l) * dum);
                for (size_t l = 0; l < m; ++l) gsl_matrix_set(b, ll, l, gsl_matrix_get(b, ll, l) - gsl_matrix_get(b, icol, l) * dum);
            }
        }
    }

    for (size_t l = n; l-- > 0;) {
        if (gsl_vector_get(indxr, l) != gsl_vector_get(indxc, l)) {
            for (size_t k = 0; k < n; ++k) {
                double temp = gsl_matrix_get(a, k, gsl_vector_get(indxr, l));
                gsl_matrix_set(a, k, gsl_vector_get(indxr, l), gsl_matrix_get(a, k, gsl_vector_get(indxc, l)));
                gsl_matrix_set(a, k, gsl_vector_get(indxc, l), temp);
            }
        }
    }

    gsl_vector_free(indxc);
    gsl_vector_free(indxr);
    gsl_vector_free(ipiv);
}

void gaussj_single(gsl_matrix* a) {
    gsl_matrix* b = gsl_matrix_alloc(a->size1, 0);
    gaussj(a, b);
    gsl_matrix_free(b);
}
