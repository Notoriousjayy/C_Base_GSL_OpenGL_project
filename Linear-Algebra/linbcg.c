//
// Created by jorda on 7/15/2024.
//

#include "linbcg.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define SQR(a) ((a)*(a))
#define EPS 1.0e-14

void linbcg_solve(Linbcg *self, const gsl_vector *b, gsl_vector *x, const int itol, const double tol, const int itmax, int *iter, double *err) {
    double ak, akden, bk, bkden = 1.0, bknum, bnrm, dxnrm, xnrm, zm1nrm, znrm;
    int j, n = b->size;
    gsl_vector *p = gsl_vector_alloc(n);
    gsl_vector *pp = gsl_vector_alloc(n);
    gsl_vector *r = gsl_vector_alloc(n);
    gsl_vector *rr = gsl_vector_alloc(n);
    gsl_vector *z = gsl_vector_alloc(n);
    gsl_vector *zz = gsl_vector_alloc(n);

    *iter = 0;
    self->atimes(x, r, 0, self);
    for (j = 0; j < n; j++) {
        gsl_vector_set(r, j, gsl_vector_get(b, j) - gsl_vector_get(r, j));
        gsl_vector_set(rr, j, gsl_vector_get(r, j));
    }

    if (itol == 1) {
        bnrm = self->snrm(b, itol);
        self->asolve(r, z, 0, self);
    } else if (itol == 2) {
        self->asolve(b, z, 0, self);
        bnrm = self->snrm(z, itol);
        self->asolve(r, z, 0, self);
    } else if (itol == 3 || itol == 4) {
        self->asolve(b, z, 0, self);
        bnrm = self->snrm(z, itol);
        self->asolve(r, z, 0, self);
        znrm = self->snrm(z, itol);
    } else {
        fprintf(stderr, "illegal itol in linbcg\n");
        exit(EXIT_FAILURE);
    }

    while (*iter < itmax) {
        (*iter)++;
        self->asolve(rr, zz, 1, self);
        for (bknum = 0.0, j = 0; j < n; j++) {
            bknum += gsl_vector_get(z, j) * gsl_vector_get(rr, j);
        }
        if (*iter == 1) {
            for (j = 0; j < n; j++) {
                gsl_vector_set(p, j, gsl_vector_get(z, j));
                gsl_vector_set(pp, j, gsl_vector_get(zz, j));
            }
        } else {
            bk = bknum / bkden;
            for (j = 0; j < n; j++) {
                gsl_vector_set(p, j, bk * gsl_vector_get(p, j) + gsl_vector_get(z, j));
                gsl_vector_set(pp, j, bk * gsl_vector_get(pp, j) + gsl_vector_get(zz, j));
            }
        }
        bkden = bknum;
        self->atimes(p, z, 0, self);
        for (akden = 0.0, j = 0; j < n; j++) {
            akden += gsl_vector_get(z, j) * gsl_vector_get(pp, j);
        }
        ak = bknum / akden;
        self->atimes(pp, zz, 1, self);
        for (j = 0; j < n; j++) {
            gsl_vector_set(x, j, gsl_vector_get(x, j) + ak * gsl_vector_get(p, j));
            gsl_vector_set(r, j, gsl_vector_get(r, j) - ak * gsl_vector_get(z, j));
            gsl_vector_set(rr, j, gsl_vector_get(rr, j) - ak * gsl_vector_get(zz, j));
        }
        self->asolve(r, z, 0, self);
        if (itol == 1) {
            *err = self->snrm(r, itol) / bnrm;
        } else if (itol == 2) {
            *err = self->snrm(z, itol) / bnrm;
        } else if (itol == 3 || itol == 4) {
            zm1nrm = znrm;
            znrm = self->snrm(z, itol);
            if (fabs(zm1nrm - znrm) > EPS * znrm) {
                dxnrm = fabs(ak) * self->snrm(p, itol);
                *err = znrm / fabs(zm1nrm - znrm) * dxnrm;
            } else {
                *err = znrm / bnrm;
                continue;
            }
            xnrm = self->snrm(x, itol);
            if (*err <= 0.5 * xnrm) *err /= xnrm;
            else {
                *err = znrm / bnrm;
                continue;
            }
        }
        if (*err <= tol) break;
    }

    gsl_vector_free(p);
    gsl_vector_free(pp);
    gsl_vector_free(r);
    gsl_vector_free(rr);
    gsl_vector_free(z);
    gsl_vector_free(zz);
}

double linbcg_snrm(const gsl_vector *sx, const int itol) {
    int i, n = sx->size;
    double ans;
    if (itol <= 3) {
        ans = 0.0;
        for (i = 0; i < n; i++) ans += SQR(gsl_vector_get(sx, i));
        return sqrt(ans);
    } else {
        int isamax = 0;
        for (i = 0; i < n; i++) {
            if (fabs(gsl_vector_get(sx, i)) > fabs(gsl_vector_get(sx, isamax))) isamax = i;
        }
        return fabs(gsl_vector_get(sx, isamax));
    }
}
