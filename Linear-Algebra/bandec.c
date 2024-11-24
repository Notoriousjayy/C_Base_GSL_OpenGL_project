//
// Created by jorda on 7/15/2024.
//
#include "../include/bandec.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

void banmul(const gsl_matrix *a, int m1, int m2, const gsl_vector *x, gsl_vector *b) {
    int n = a->size1;
    for (int i = 0; i < n; i++) {
        int k = i - m1;
        int tmploop = ((m1 + m2 + 1) < (n - k)) ? (m1 + m2 + 1) : (n - k);
        gsl_vector_set(b, i, 0.0);
        for (int j = ((0) > (-k)) ? (0) : (-k); j < tmploop; j++) {
            int x_index = j + k;
            if (x_index >= 0 && x_index < x->size) {
                gsl_vector_set(b, i, gsl_vector_get(b, i) + gsl_matrix_get(a, i, j) * gsl_vector_get(x, x_index));
            }
        }
    }
}

void solve_band_system(const gsl_matrix *a, const gsl_vector *b, gsl_vector *x, int m1, int m2) {
    int n = a->size1;
    gsl_matrix *ab = gsl_matrix_alloc(n, 2 * m1 + m2 + 1);
    gsl_vector_uint *p = gsl_vector_uint_alloc(n);
    int signum;

    // Pack the matrix 'a' into banded format 'ab'
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m1 + m2 + 1; j++) {
            int col = i + j - m1;
            if (col >= 0 && col < n) {
                gsl_matrix_set(ab, i, j, gsl_matrix_get(a, i, col));
            } else {
                gsl_matrix_set(ab, i, j, 0.0);
            }
        }
    }

    gsl_linalg_LU_band_decomp(n, m1, m2, ab, p);
    gsl_linalg_LU_band_solve(m1, m2, ab, p, b, x);

    gsl_matrix_free(ab);
    gsl_vector_uint_free(p);
}

double determinant_band_matrix(const gsl_matrix *a, int m1, int m2) {
    int n = a->size1;
    gsl_matrix *ab = gsl_matrix_alloc(n, 2 * m1 + m2 + 1);
    gsl_vector_uint *p = gsl_vector_uint_alloc(n);
    int signum;

    // Pack the matrix 'a' into banded format 'ab'
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m1 + m2 + 1; j++) {
            int col = i + j - m1;
            if (col >= 0 && col < n) {
                gsl_matrix_set(ab, i, j, gsl_matrix_get(a, i, col));
            } else {
                gsl_matrix_set(ab, i, j, 0.0);
            }
        }
    }

    gsl_linalg_LU_band_decomp(n, m1, m2, ab, p);
    double det = gsl_linalg_LU_det(ab, signum);

    gsl_matrix_free(ab);
    gsl_vector_uint_free(p);

    return det;
}
