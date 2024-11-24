//
// Created by jorda on 7/15/2024.
//

#ifndef LU_DECOMP_H
#define LU_DECOMP_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>

typedef struct {
    int n;
    gsl_matrix* lu;
    gsl_permutation* perm;
    double d;
} LUdcmp;

LUdcmp* lu_decomp_create(const gsl_matrix* a);
void lu_decomp_free(LUdcmp* ludcmp);
void lu_decomp_solve(const LUdcmp* ludcmp, const gsl_vector* b, gsl_vector* x);
void lu_decomp_solve_matrix(const LUdcmp* ludcmp, const gsl_matrix* b, gsl_matrix* x);
void lu_decomp_inverse(const LUdcmp* ludcmp, gsl_matrix* ainv);
double lu_decomp_det(const LUdcmp* ludcmp);
void lu_decomp_mprove(const LUdcmp* ludcmp, const gsl_vector* b, gsl_vector* x);

#endif // LU_DECOMP_H
