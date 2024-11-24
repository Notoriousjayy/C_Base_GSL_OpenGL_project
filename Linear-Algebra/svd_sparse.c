//
// Created by jorda on 7/15/2024.
//
#include "../include/svd_sparse.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <gsl/gsl_linalg.h>

static double pythag(double a, double b) {
    double absa = fabs(a), absb = fabs(b);
    if (absa > absb) return absa * sqrt(1.0 + (absb / absa) * (absb / absa));
    else return (absb == 0.0) ? 0.0 : absb * sqrt(1.0 + (absa / absb) * (absa / absb));
}

static void svd_sparse_decompose(SVD_Sparse *svd_sparse) {
    gsl_matrix *A_copy = gsl_matrix_alloc(svd_sparse->m, svd_sparse->n);
    gsl_matrix_memcpy(A_copy, svd_sparse->u);
    gsl_vector *work = gsl_vector_alloc(svd_sparse->n);
    gsl_linalg_SV_decomp(A_copy, svd_sparse->v, svd_sparse->w, work);
    gsl_matrix_memcpy(svd_sparse->u, A_copy);
    gsl_matrix_free(A_copy);
    gsl_vector_free(work);
}

static void svd_sparse_reorder(SVD_Sparse *svd_sparse) {
    int i, j, k, inc = 1;
    double sw;
    gsl_vector *su = gsl_vector_alloc(svd_sparse->m);
    gsl_vector *sv = gsl_vector_alloc(svd_sparse->n);

    do { inc *= 3; inc++; } while (inc <= svd_sparse->n);
    do {
        inc /= 3;
        for (i = inc; i < svd_sparse->n; i++) {
            sw = gsl_vector_get(svd_sparse->w, i);
            for (k = 0; k < svd_sparse->m; k++) gsl_vector_set(su, k, gsl_matrix_get(svd_sparse->u, k, i));
            for (k = 0; k < svd_sparse->n; k++) gsl_vector_set(sv, k, gsl_matrix_get(svd_sparse->v, k, i));
            j = i;
            while (j >= inc && gsl_vector_get(svd_sparse->w, j - inc) < sw) {
                gsl_vector_set(svd_sparse->w, j, gsl_vector_get(svd_sparse->w, j - inc));
                for (k = 0; k < svd_sparse->m; k++) gsl_matrix_set(svd_sparse->u, k, j, gsl_matrix_get(svd_sparse->u, k, j - inc));
                for (k = 0; k < svd_sparse->n; k++) gsl_matrix_set(svd_sparse->v, k, j, gsl_matrix_get(svd_sparse->v, k, j - inc));
                j -= inc;
            }
            gsl_vector_set(svd_sparse->w, j, sw);
            for (k = 0; k < svd_sparse->m; k++) gsl_matrix_set(svd_sparse->u, k, j, gsl_vector_get(su, k));
            for (k = 0; k < svd_sparse->n; k++) gsl_matrix_set(svd_sparse->v, k, j, gsl_vector_get(sv, k));
        }
    } while (inc > 1);

    for (k = 0; k < svd_sparse->n; k++) {
        int s = 0;
        for (i = 0; i < svd_sparse->m; i++) if (gsl_matrix_get(svd_sparse->u, i, k) < 0.0) s++;
        for (j = 0; j < svd_sparse->n; j++) if (gsl_matrix_get(svd_sparse->v, j, k) < 0.0) s++;
        if (s > (svd_sparse->m + svd_sparse->n) / 2) {
            for (i = 0; i < svd_sparse->m; i++) gsl_matrix_set(svd_sparse->u, i, k, -gsl_matrix_get(svd_sparse->u, i, k));
            for (j = 0; j < svd_sparse->n; j++) gsl_matrix_set(svd_sparse->v, j, k, -gsl_matrix_get(svd_sparse->v, j, k));
        }
    }

    gsl_vector_free(su);
    gsl_vector_free(sv);
}

SVD_Sparse* svd_sparse_create(const gsl_matrix *a) {
    SVD_Sparse *svd_sparse = (SVD_Sparse*)malloc(sizeof(SVD_Sparse));
    svd_sparse->m = a->size1;
    svd_sparse->n = a->size2;
    svd_sparse->u = gsl_matrix_alloc(svd_sparse->m, svd_sparse->n);
    svd_sparse->v = gsl_matrix_alloc(svd_sparse->n, svd_sparse->n);
    svd_sparse->w = gsl_vector_alloc(svd_sparse->n);
    svd_sparse->eps = DBL_EPSILON;

    gsl_matrix_memcpy(svd_sparse->u, a);
    svd_sparse_decompose(svd_sparse);
    svd_sparse_reorder(svd_sparse);

    svd_sparse->tsh = 0.5 * sqrt((double)(svd_sparse->m + svd_sparse->n + 1)) * gsl_vector_get(svd_sparse->w, 0) * svd_sparse->eps;

    return svd_sparse;
}

void svd_sparse_free(SVD_Sparse *svd_sparse) {
    gsl_matrix_free(svd_sparse->u);
    gsl_matrix_free(svd_sparse->v);
    gsl_vector_free(svd_sparse->w);
    free(svd_sparse);
}

void svd_sparse_solve(const SVD_Sparse *svd_sparse, const gsl_vector *b, gsl_vector *x, double thresh) {
    int i, j;
    double s;
    gsl_vector *tmp = gsl_vector_alloc(svd_sparse->n);

    if (thresh < 0.0) thresh = svd_sparse->tsh;
    for (j = 0; j < svd_sparse->n; j++) {
        s = 0.0;
        if (gsl_vector_get(svd_sparse->w, j) > thresh) {
            for (i = 0; i < svd_sparse->m; i++) s += gsl_matrix_get(svd_sparse->u, i, j) * gsl_vector_get(b, i);
            s /= gsl_vector_get(svd_sparse->w, j);
        }
        gsl_vector_set(tmp, j, s);
    }

    for (j = 0; j < svd_sparse->n; j++) {
        s = 0.0;
        for (i = 0; i < svd_sparse->n; i++) s += gsl_matrix_get(svd_sparse->v, j, i) * gsl_vector_get(tmp, i);
        gsl_vector_set(x, j, s);
    }

    gsl_vector_free(tmp);
}

void svd_sparse_solve_matrix(const SVD_Sparse *svd_sparse, const gsl_matrix *b, gsl_matrix *x, double thresh) {
    int i, j, p = b->size2;
    gsl_vector *xx = gsl_vector_alloc(svd_sparse->n);
    gsl_vector *bcol = gsl_vector_alloc(svd_sparse->m);

    if (b->size1 != svd_sparse->m || x->size1 != svd_sparse->n || x->size2 != p) {
        fprintf(stderr, "svd_sparse_solve_matrix: bad sizes\n");
        exit(EXIT_FAILURE);
    }

    for (j = 0; j < p; j++) {
        for (i = 0; i < svd_sparse->m; i++) gsl_vector_set(bcol, i, gsl_matrix_get(b, i, j));
        svd_sparse_solve(svd_sparse, bcol, xx, thresh);
        for (i = 0; i < svd_sparse->n; i++) gsl_matrix_set(x, i, j, gsl_vector_get(xx, i));
    }

    gsl_vector_free(xx);
    gsl_vector_free(bcol);
}

int svd_sparse_rank(const SVD_Sparse *svd_sparse, double thresh) {
    int j, nr = 0;
    if (thresh < 0.0) thresh = svd_sparse->tsh;
    for (j = 0; j < svd_sparse->n; j++) if (gsl_vector_get(svd_sparse->w, j) > thresh) nr++;
    return nr;
}

int svd_sparse_nullity(const SVD_Sparse *svd_sparse, double thresh) {
    int j, nn = 0;
    if (thresh < 0.0) thresh = svd_sparse->tsh;
    for (j = 0; j < svd_sparse->n; j++) if (gsl_vector_get(svd_sparse->w, j) <= thresh) nn++;
    return nn;
}

gsl_matrix* svd_sparse_range(const SVD_Sparse *svd_sparse, double thresh) {
    int i, j, nr = 0;
    gsl_matrix *range;
    if (thresh < 0.0) thresh = svd_sparse->tsh;

    range = gsl_matrix_alloc(svd_sparse->m, svd_sparse_rank(svd_sparse, thresh));
    for (j = 0; j < svd_sparse->n; j++) {
        if (gsl_vector_get(svd_sparse->w, j) > thresh) {
            for (i = 0; i < svd_sparse->m; i++) gsl_matrix_set(range, i, nr, gsl_matrix_get(svd_sparse->u, i, j));
            nr++;
        }
    }
    return range;
}

gsl_matrix* svd_sparse_nullspace(const SVD_Sparse *svd_sparse, double thresh) {
    int j, jj, nn = 0;
    gsl_matrix *nullsp;
    if (thresh < 0.0) thresh = svd_sparse->tsh;

    nullsp = gsl_matrix_alloc(svd_sparse->n, svd_sparse_nullity(svd_sparse, thresh));
    for (j = 0; j < svd_sparse->n; j++) {
        if (gsl_vector_get(svd_sparse->w, j) <= thresh) {
            for (jj = 0; jj < svd_sparse->n; jj++) gsl_matrix_set(nullsp, jj, nn, gsl_matrix_get(svd_sparse->v, jj, j));
            nn++;
        }
    }
    return nullsp;
}

double svd_sparse_inv_condition(const SVD_Sparse *svd_sparse) {
    if (gsl_vector_get(svd_sparse->w, 0) <= 0.0 || gsl_vector_get(svd_sparse->w, svd_sparse->n - 1) <= 0.0) return 0.0;
    return gsl_vector_get(svd_sparse->w, svd_sparse->n - 1) / gsl_vector_get(svd_sparse->w, 0);
}
