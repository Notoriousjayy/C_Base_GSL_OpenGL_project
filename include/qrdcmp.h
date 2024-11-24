//
// Created by jorda on 7/16/2024.
//

#ifndef QRDCMP_H
#define QRDCMP_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct {
    int n;
    gsl_matrix *qt;
    gsl_matrix *r;
    int sing;
} QRdcmp;

QRdcmp* qrdcmp_alloc(const gsl_matrix *a);
void qrdcmp_free(QRdcmp *qr);
void qrdcmp_solve(const QRdcmp *qr, const gsl_vector *b, gsl_vector *x);
void qrdcmp_qtmult(const QRdcmp *qr, const gsl_vector *b, gsl_vector *x);
void qrdcmp_rsolve(const QRdcmp *qr, const gsl_vector *b, gsl_vector *x);
void qrdcmp_update(QRdcmp *qr, const gsl_vector *u, const gsl_vector *v);
void qrdcmp_rotate(QRdcmp *qr, const int i, const double a, const double b);

#endif // QRDCMP_H
