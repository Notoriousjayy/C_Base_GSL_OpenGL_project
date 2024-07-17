//
// Created by jorda on 7/15/2024.
//

#ifndef SVD_H
#define SVD_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct {
    int m, n;
    gsl_matrix *u, *v;
    gsl_vector *w;
    double eps, tsh;
} SVD;

SVD* svd_create(const gsl_matrix *a);
void svd_free(SVD *svd);
void svd_solve(const SVD *svd, const gsl_vector *b, gsl_vector *x, double thresh);
void svd_solve_matrix(const SVD *svd, const gsl_matrix *b, gsl_matrix *x, double thresh);
int svd_rank(const SVD *svd, double thresh);
int svd_nullity(const SVD *svd, double thresh);
gsl_matrix* svd_range(const SVD *svd, double thresh);
gsl_matrix* svd_nullspace(const SVD *svd, double thresh);
double svd_inv_condition(const SVD *svd);

#endif // SVD_H
