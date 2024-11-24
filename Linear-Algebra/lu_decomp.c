//
// Created by jorda on 7/15/2024.
//
#include "../include/lu_decomp.h"
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

LUdcmp* lu_decomp_create(const gsl_matrix* a) {
    LUdcmp* ludcmp = (LUdcmp*)malloc(sizeof(LUdcmp));
    ludcmp->n = a->size1;
    ludcmp->lu = gsl_matrix_alloc(ludcmp->n, ludcmp->n);
    gsl_matrix_memcpy(ludcmp->lu, a);
    ludcmp->perm = gsl_permutation_alloc(ludcmp->n);
    gsl_linalg_LU_decomp(ludcmp->lu, ludcmp->perm, &(ludcmp->d));
    return ludcmp;
}

void lu_decomp_free(LUdcmp* ludcmp) {
    gsl_matrix_free(ludcmp->lu);
    gsl_permutation_free(ludcmp->perm);
    free(ludcmp);
}

void lu_decomp_solve(const LUdcmp* ludcmp, const gsl_vector* b, gsl_vector* x) {
    gsl_linalg_LU_solve(ludcmp->lu, ludcmp->perm, b, x);
}

void lu_decomp_solve_matrix(const LUdcmp* ludcmp, const gsl_matrix* b, gsl_matrix* x) {
    gsl_matrix* xx = gsl_matrix_alloc(ludcmp->n, 1);
    gsl_vector* col_b = gsl_vector_alloc(ludcmp->n);
    gsl_vector* col_x = gsl_vector_alloc(ludcmp->n);
    for (size_t j = 0; j < b->size2; j++) {
        gsl_matrix_get_col(col_b, b, j);
        lu_decomp_solve(ludcmp, col_b, col_x);
        gsl_matrix_set_col(x, j, col_x);
    }
    gsl_vector_free(col_b);
    gsl_vector_free(col_x);
    gsl_matrix_free(xx);
}

void lu_decomp_inverse(const LUdcmp* ludcmp, gsl_matrix* ainv) {
    gsl_matrix* identity = gsl_matrix_alloc(ludcmp->n, ludcmp->n);
    gsl_matrix_set_identity(identity);
    lu_decomp_solve_matrix(ludcmp, identity, ainv);
    gsl_matrix_free(identity);
}

double lu_decomp_det(const LUdcmp* ludcmp) {
    return gsl_linalg_LU_det(ludcmp->lu, ludcmp->d);
}

void lu_decomp_mprove(const LUdcmp* ludcmp, const gsl_vector* b, gsl_vector* x) {
    gsl_vector* r = gsl_vector_alloc(ludcmp->n);
    gsl_blas_dgemv(CblasNoTrans, -1.0, ludcmp->lu, x, 1.0, r);
    gsl_vector_sub(r, b);
    lu_decomp_solve(ludcmp, r, r);
    gsl_vector_sub(x, r);
    gsl_vector_free(r);
}
