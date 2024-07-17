//
// Created by jorda on 7/15/2024.
//

#include "svd.h"
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

static void svd_decompose(SVD *svd) {
    gsl_vector *work = gsl_vector_alloc(svd->n);
    gsl_linalg_SV_decomp(svd->u, svd->v, svd->w, work);
    gsl_vector_free(work);
}

static void svd_reorder(SVD *svd) {
    int i, j, k, inc = 1;
    double sw;
    gsl_vector *su = gsl_vector_alloc(svd->m);
    gsl_vector *sv = gsl_vector_alloc(svd->n);

    do { inc *= 3; inc++; } while (inc <= svd->n);
    do {
        inc /= 3;
        for (i = inc; i < svd->n; i++) {
            sw = gsl_vector_get(svd->w, i);
            for (k = 0; k < svd->m; k++) gsl_vector_set(su, k, gsl_matrix_get(svd->u, k, i));
            for (k = 0; k < svd->n; k++) gsl_vector_set(sv, k, gsl_matrix_get(svd->v, k, i));
            j = i;
            while (j >= inc && gsl_vector_get(svd->w, j - inc) < sw) {
                gsl_vector_set(svd->w, j, gsl_vector_get(svd->w, j - inc));
                for (k = 0; k < svd->m; k++) gsl_matrix_set(svd->u, k, j, gsl_matrix_get(svd->u, k, j - inc));
                for (k = 0; k < svd->n; k++) gsl_matrix_set(svd->v, k, j, gsl_matrix_get(svd->v, k, j - inc));
                j -= inc;
            }
            gsl_vector_set(svd->w, j, sw);
            for (k = 0; k < svd->m; k++) gsl_matrix_set(svd->u, k, j, gsl_vector_get(su, k));
            for (k = 0; k < svd->n; k++) gsl_matrix_set(svd->v, k, j, gsl_vector_get(sv, k));
        }
    } while (inc > 1);

    for (k = 0; k < svd->n; k++) {
        int s = 0;
        for (i = 0; i < svd->m; i++) if (gsl_matrix_get(svd->u, i, k) < 0.0) s++;
        for (j = 0; j < svd->n; j++) if (gsl_matrix_get(svd->v, j, k) < 0.0) s++;
        if (s > (svd->m + svd->n) / 2) {
            for (i = 0; i < svd->m; i++) gsl_matrix_set(svd->u, i, k, -gsl_matrix_get(svd->u, i, k));
            for (j = 0; j < svd->n; j++) gsl_matrix_set(svd->v, j, k, -gsl_matrix_get(svd->v, j, k));
        }
    }

    gsl_vector_free(su);
    gsl_vector_free(sv);
}

SVD* svd_create(const gsl_matrix *a) {
    SVD *svd = (SVD*)malloc(sizeof(SVD));
    svd->m = a->size1;
    svd->n = a->size2;
    svd->u = gsl_matrix_alloc(svd->m, svd->n);
    svd->v = gsl_matrix_alloc(svd->n, svd->n);
    svd->w = gsl_vector_alloc(svd->n);
    svd->eps = DBL_EPSILON;

    gsl_matrix_memcpy(svd->u, a);
    svd_decompose(svd);
    svd_reorder(svd);

    svd->tsh = 0.5 * sqrt((double)(svd->m + svd->n + 1)) * gsl_vector_get(svd->w, 0) * svd->eps;

    return svd;
}

void svd_free(SVD *svd) {
    gsl_matrix_free(svd->u);
    gsl_matrix_free(svd->v);
    gsl_vector_free(svd->w);
    free(svd);
}

void svd_solve(const SVD *svd, const gsl_vector *b, gsl_vector *x, double thresh) {
    int i, j;
    double s;
    gsl_vector *tmp = gsl_vector_alloc(svd->n);

    if (thresh < 0.0) thresh = svd->tsh;
    for (j = 0; j < svd->n; j++) {
        s = 0.0;
        if (gsl_vector_get(svd->w, j) > thresh) {
            for (i = 0; i < svd->m; i++) s += gsl_matrix_get(svd->u, i, j) * gsl_vector_get(b, i);
            s /= gsl_vector_get(svd->w, j);
        }
        gsl_vector_set(tmp, j, s);
    }

    for (j = 0; j < svd->n; j++) {
        s = 0.0;
        for (i = 0; i < svd->n; i++) s += gsl_matrix_get(svd->v, j, i) * gsl_vector_get(tmp, i);
        gsl_vector_set(x, j, s);
    }

    gsl_vector_free(tmp);
}

void svd_solve_matrix(const SVD *svd, const gsl_matrix *b, gsl_matrix *x, double thresh) {
    int i, j, p = b->size2;
    gsl_vector *xx = gsl_vector_alloc(svd->n);
    gsl_vector *bcol = gsl_vector_alloc(svd->m);

    if (b->size1 != svd->m || x->size1 != svd->n || x->size2 != p) {
        fprintf(stderr, "svd_solve_matrix: bad sizes\n");
        exit(EXIT_FAILURE);
    }

    for (j = 0; j < p; j++) {
        for (i = 0; i < svd->m; i++) gsl_vector_set(bcol, i, gsl_matrix_get(b, i, j));
        svd_solve(svd, bcol, xx, thresh);
        for (i = 0; i < svd->n; i++) gsl_matrix_set(x, i, j, gsl_vector_get(xx, i));
    }

    gsl_vector_free(xx);
    gsl_vector_free(bcol);
}

int svd_rank(const SVD *svd, double thresh) {
    int j, nr = 0;
    if (thresh < 0.0) thresh = svd->tsh;
    for (j = 0; j < svd->n; j++) if (gsl_vector_get(svd->w, j) > thresh) nr++;
    return nr;
}

int svd_nullity(const SVD *svd, double thresh) {
    int j, nn = 0;
    if (thresh < 0.0) thresh = svd->tsh;
    for (j = 0; j < svd->n; j++) if (gsl_vector_get(svd->w, j) <= thresh) nn++;
    return nn;
}

gsl_matrix* svd_range(const SVD *svd, double thresh) {
    int i, j, nr = 0;
    gsl_matrix *range;
    if (thresh < 0.0) thresh = svd->tsh;

    range = gsl_matrix_alloc(svd->m, svd_rank(svd, thresh));
    for (j = 0; j < svd->n; j++) {
        if (gsl_vector_get(svd->w, j) > thresh) {
            for (i = 0; i < svd->m; i++) gsl_matrix_set(range, i, nr, gsl_matrix_get(svd->u, i, j));
            nr++;
        }
    }
    return range;
}

gsl_matrix* svd_nullspace(const SVD *svd, double thresh) {
    int j, jj, nn = 0;
    gsl_matrix *nullsp;
    if (thresh < 0.0) thresh = svd->tsh;

    nullsp = gsl_matrix_alloc(svd->n, svd_nullity(svd, thresh));
    for (j = 0; j < svd->n; j++) {
        if (gsl_vector_get(svd->w, j) <= thresh) {
            for (jj = 0; jj < svd->n; jj++) gsl_matrix_set(nullsp, jj, nn, gsl_matrix_get(svd->v, jj, j));
            nn++;
        }
    }
    return nullsp;
}

double svd_inv_condition(const SVD *svd) {
    if (gsl_vector_get(svd->w, 0) <= 0.0 || gsl_vector_get(svd->w, svd->n - 1) <= 0.0) return 0.0;
    return gsl_vector_get(svd->w, svd->n - 1) / gsl_vector_get(svd->w, 0);
}
