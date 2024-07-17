//
// Created by jorda on 7/15/2024.
//

#ifndef SVD_SPARSE_H
#define SVD_SPARSE_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>

typedef struct {
    int m, n;
    gsl_matrix *u;
    gsl_matrix *v;
    gsl_vector *w;
    double eps, tsh;
} SVD_Sparse;

SVD_Sparse* svd_sparse_create(const gsl_matrix *a);
void svd_sparse_free(SVD_Sparse *svd);
void svd_sparse_solve(const SVD_Sparse *svd, const gsl_vector *b, gsl_vector *x, double thresh);
void svd_sparse_solve_matrix(const SVD_Sparse *svd, const gsl_matrix *b, gsl_matrix *x, double thresh);
int svd_sparse_rank(const SVD_Sparse *svd, double thresh);
int svd_sparse_nullity(const SVD_Sparse *svd, double thresh);
gsl_matrix* svd_sparse_range(const SVD_Sparse *svd, double thresh);
gsl_matrix* svd_sparse_nullspace(const SVD_Sparse *svd, double thresh);
double svd_sparse_inv_condition(const SVD_Sparse *svd);

#endif // SVD_SPARSE_H
