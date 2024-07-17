//
// Created by jorda on 7/16/2024.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "qrdcmp.h"

#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SQR(x) ((x) * (x))

QRdcmp* qrdcmp_alloc(const gsl_matrix *a) {
    int i, j, k;
    int n = a->size1;
    double scale, sigma, sum, tau;
    gsl_vector *c = gsl_vector_alloc(n);
    gsl_vector *d = gsl_vector_alloc(n);

    QRdcmp *qr = (QRdcmp*)malloc(sizeof(QRdcmp));
    qr->n = n;
    qr->qt = gsl_matrix_alloc(n, n);
    qr->r = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(qr->r, a);
    qr->sing = 0;

    for (k = 0; k < n - 1; k++) {
        scale = 0.0;
        for (i = k; i < n; i++) scale = fmax(scale, fabs(gsl_matrix_get(qr->r, i, k)));
        if (scale == 0.0) {
            qr->sing = 1;
            gsl_vector_set(c, k, 0.0);
            gsl_vector_set(d, k, 0.0);
        } else {
            for (i = k; i < n; i++) gsl_matrix_set(qr->r, i, k, gsl_matrix_get(qr->r, i, k) / scale);
            for (sum = 0.0, i = k; i < n; i++) sum += SQR(gsl_matrix_get(qr->r, i, k));
            sigma = SIGN(sqrt(sum), gsl_matrix_get(qr->r, k, k));
            gsl_matrix_set(qr->r, k, k, gsl_matrix_get(qr->r, k, k) + sigma);
            gsl_vector_set(c, k, sigma * gsl_matrix_get(qr->r, k, k));
            gsl_vector_set(d, k, -scale * sigma);
            for (j = k + 1; j < n; j++) {
                for (sum = 0.0, i = k; i < n; i++) sum += gsl_matrix_get(qr->r, i, k) * gsl_matrix_get(qr->r, i, j);
                tau = sum / gsl_vector_get(c, k);
                for (i = k; i < n; i++) gsl_matrix_set(qr->r, i, j, gsl_matrix_get(qr->r, i, j) - tau * gsl_matrix_get(qr->r, i, k));
            }
        }
    }
    gsl_vector_set(d, n - 1, gsl_matrix_get(qr->r, n - 1, n - 1));
    if (gsl_vector_get(d, n - 1) == 0.0) qr->sing = 1;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) gsl_matrix_set(qr->qt, i, j, 0.0);
        gsl_matrix_set(qr->qt, i, i, 1.0);
    }
    for (k = 0; k < n - 1; k++) {
        if (gsl_vector_get(c, k) != 0.0) {
            for (j = 0; j < n; j++) {
                sum = 0.0;
                for (i = k; i < n; i++) sum += gsl_matrix_get(qr->r, i, k) * gsl_matrix_get(qr->qt, i, j);
                sum /= gsl_vector_get(c, k);
                for (i = k; i < n; i++) gsl_matrix_set(qr->qt, i, j, gsl_matrix_get(qr->qt, i, j) - sum * gsl_matrix_get(qr->r, i, k));
            }
        }
    }
    for (i = 0; i < n; i++) {
        gsl_matrix_set(qr->r, i, i, gsl_vector_get(d, i));
        for (j = 0; j < i; j++) gsl_matrix_set(qr->r, i, j, 0.0);
    }

    gsl_vector_free(c);
    gsl_vector_free(d);

    return qr;
}

void qrdcmp_free(QRdcmp *qr) {
    gsl_matrix_free(qr->qt);
    gsl_matrix_free(qr->r);
    free(qr);
}

void qrdcmp_solve(const QRdcmp *qr, const gsl_vector *b, gsl_vector *x) {
    qrdcmp_qtmult(qr, b, x);
    qrdcmp_rsolve(qr, x, x);
}

void qrdcmp_qtmult(const QRdcmp *qr, const gsl_vector *b, gsl_vector *x) {
    int i, j;
    double sum;
    for (i = 0; i < qr->n; i++) {
        sum = 0.0;
        for (j = 0; j < qr->n; j++) sum += gsl_matrix_get(qr->qt, i, j) * gsl_vector_get(b, j);
        gsl_vector_set(x, i, sum);
    }
}

void qrdcmp_rsolve(const QRdcmp *qr, const gsl_vector *b, gsl_vector *x) {
    int i, j;
    double sum;
    if (qr->sing) {
        fprintf(stderr, "Attempting to solve in a singular QR\n");
        exit(EXIT_FAILURE);
    }
    for (i = qr->n - 1; i >= 0; i--) {
        sum = gsl_vector_get(b, i);
        for (j = i + 1; j < qr->n; j++) sum -= gsl_matrix_get(qr->r, i, j) * gsl_vector_get(x, j);
        gsl_vector_set(x, i, sum / gsl_matrix_get(qr->r, i, i));
    }
}

void qrdcmp_update(QRdcmp *qr, const gsl_vector *u, const gsl_vector *v) {
    int i, k;
    gsl_vector *w = gsl_vector_alloc(qr->n);
    gsl_vector_memcpy(w, u);

    for (k = qr->n - 1; k >= 0; k--) {
        if (gsl_vector_get(w, k) != 0.0) break;
    }
    if (k < 0) k = 0;
    for (i = k - 1; i >= 0; i--) {
        qrdcmp_rotate(qr, i, gsl_vector_get(w, i), -gsl_vector_get(w, i + 1));
        if (gsl_vector_get(w, i) == 0.0)
            gsl_vector_set(w, i, fabs(gsl_vector_get(w, i + 1)));
        else if (fabs(gsl_vector_get(w, i)) > fabs(gsl_vector_get(w, i + 1)))
            gsl_vector_set(w, i, fabs(gsl_vector_get(w, i)) * sqrt(1.0 + SQR(gsl_vector_get(w, i + 1) / gsl_vector_get(w, i))));
        else gsl_vector_set(w, i, fabs(gsl_vector_get(w, i + 1)) * sqrt(1.0 + SQR(gsl_vector_get(w, i) / gsl_vector_get(w, i + 1))));
    }
    for (i = 0; i < qr->n; i++) gsl_matrix_set(qr->r, 0, i, gsl_matrix_get(qr->r, 0, i) + gsl_vector_get(w, 0) * gsl_vector_get(v, i));
    for (i = 0; i < k; i++) qrdcmp_rotate(qr, i, gsl_matrix_get(qr->r, i, i), -gsl_matrix_get(qr->r, i + 1, i));
    for (i = 0; i < qr->n; i++) if (gsl_matrix_get(qr->r, i, i) == 0.0) qr->sing = 1;

    gsl_vector_free(w);
}

void qrdcmp_rotate(QRdcmp *qr, const int i, const double a, const double b) {
    int j;
    double c, fact, s, w, y;
    if (a == 0.0) {
        c = 0.0;
        s = (b >= 0.0 ? 1.0 : -1.0);
    } else if (fabs(a) > fabs(b)) {
        fact = b / a;
        c = SIGN(1.0 / sqrt(1.0 + (fact * fact)), a);
        s = fact * c;
    } else {
        fact = a / b;
        s = SIGN(1.0 / sqrt(1.0 + (fact * fact)), b);
        c = fact * s;
    }
    for (j = i; j < qr->n; j++) {
        y = gsl_matrix_get(qr->r, i, j);
        w = gsl_matrix_get(qr->r, i + 1, j);
        gsl_matrix_set(qr->r, i, j, c * y - s * w);
        gsl_matrix_set(qr->r, i + 1, j, s * y + c * w);
    }
    for (j = 0; j < qr->n; j++) {
        y = gsl_matrix_get(qr->qt, i, j);
        w = gsl_matrix_get(qr->qt, i + 1, j);
        gsl_matrix_set(qr->qt, i, j, c * y - s * w);
        gsl_matrix_set(qr->qt, i + 1, j, s * y + c * w);
    }
}
