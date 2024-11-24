//
// Created by jorda on 7/16/2024.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "../include/cholesky.h"

void cholesky_decomp(Cholesky *ch, const gsl_matrix *a) {
    int i, j, k;
    double sum;
    ch->n = a->size1;
    ch->el = gsl_matrix_alloc(ch->n, ch->n);
    gsl_matrix_memcpy(ch->el, a);

    for (i = 0; i < ch->n; i++) {
        for (j = i; j < ch->n; j++) {
            sum = gsl_matrix_get(ch->el, i, j);
            for (k = i - 1; k >= 0; k--) {
                sum -= gsl_matrix_get(ch->el, i, k) * gsl_matrix_get(ch->el, j, k);
            }
            if (i == j) {
                if (sum <= 0.0) {
                    fprintf(stderr, "Cholesky decomposition failed\n");
                    exit(EXIT_FAILURE);
                }
                gsl_matrix_set(ch->el, i, i, sqrt(sum));
            } else {
                gsl_matrix_set(ch->el, j, i, sum / gsl_matrix_get(ch->el, i, i));
            }
        }
    }
    for (i = 0; i < ch->n; i++) {
        for (j = 0; j < i; j++) {
            gsl_matrix_set(ch->el, j, i, 0.0);
        }
    }
}

void cholesky_solve(const Cholesky *ch, const gsl_vector *b, gsl_vector *x) {
    int i, k;
    double sum;
    if (b->size != ch->n || x->size != ch->n) {
        fprintf(stderr, "bad lengths in Cholesky solve\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < ch->n; i++) {
        sum = gsl_vector_get(b, i);
        for (k = i - 1; k >= 0; k--) {
            sum -= gsl_matrix_get(ch->el, i, k) * gsl_vector_get(x, k);
        }
        gsl_vector_set(x, i, sum / gsl_matrix_get(ch->el, i, i));
    }
    for (i = ch->n - 1; i >= 0; i--) {
        sum = gsl_vector_get(x, i);
        for (k = i + 1; k < ch->n; k++) {
            sum -= gsl_matrix_get(ch->el, k, i) * gsl_vector_get(x, k);
        }
        gsl_vector_set(x, i, sum / gsl_matrix_get(ch->el, i, i));
    }
}

void cholesky_elmult(const Cholesky *ch, const gsl_vector *y, gsl_vector *b) {
    int i, j;
    if (b->size != ch->n || y->size != ch->n) {
        fprintf(stderr, "bad lengths in Cholesky elmult\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < ch->n; i++) {
        gsl_vector_set(b, i, 0.0);
        for (j = 0; j <= i; j++) {
            gsl_vector_set(b, i, gsl_vector_get(b, i) + gsl_matrix_get(ch->el, i, j) * gsl_vector_get(y, j));
        }
    }
}

void cholesky_elsolve(const Cholesky *ch, const gsl_vector *b, gsl_vector *y) {
    int i, j;
    double sum;
    if (b->size != ch->n || y->size != ch->n) {
        fprintf(stderr, "bad lengths in Cholesky elsolve\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < ch->n; i++) {
        sum = gsl_vector_get(b, i);
        for (j = 0; j < i; j++) {
            sum -= gsl_matrix_get(ch->el, i, j) * gsl_vector_get(y, j);
        }
        gsl_vector_set(y, i, sum / gsl_matrix_get(ch->el, i, i));
    }
}

void cholesky_inverse(const Cholesky *ch, gsl_matrix *ainv) {
    int i, j, k;
    double sum;
    gsl_matrix_set_zero(ainv);
    for (i = 0; i < ch->n; i++) {
        for (j = 0; j <= i; j++) {
            sum = (i == j) ? 1.0 : 0.0;
            for (k = i - 1; k >= j; k--) {
                sum -= gsl_matrix_get(ch->el, i, k) * gsl_matrix_get(ainv, j, k);
            }
            gsl_matrix_set(ainv, j, i, sum / gsl_matrix_get(ch->el, i, i));
        }
    }
    for (i = ch->n - 1; i >= 0; i--) {
        for (j = 0; j <= i; j++) {
            sum = (i < j) ? 0.0 : gsl_matrix_get(ainv, j, i);
            for (k = i + 1; k < ch->n; k++) {
                sum -= gsl_matrix_get(ch->el, k, i) * gsl_matrix_get(ainv, j, k);
            }
            double value = sum / gsl_matrix_get(ch->el, i, i);
            gsl_matrix_set(ainv, i, j, value);
            gsl_matrix_set(ainv, j, i, value);
        }
    }
}

double cholesky_logdet(const Cholesky *ch) {
    double sum = 0.0;
    for (int i = 0; i < ch->n; i++) {
        sum += log(gsl_matrix_get(ch->el, i, i));
    }
    return 2.0 * sum;
}
