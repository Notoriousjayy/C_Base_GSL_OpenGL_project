//
// Created by jorda on 7/15/2024.
//

#ifndef LINBCG_H
#define LINBCG_H

#include <gsl/gsl_vector.h>

typedef struct Linbcg {
    void (*asolve)(const gsl_vector *b, gsl_vector *x, const int itrnsp, void *solver);
    void (*atimes)(const gsl_vector *x, gsl_vector *r, const int itrnsp, void *solver);
    void (*solve)(struct Linbcg *self, const gsl_vector *b, gsl_vector *x, const int itol, const double tol, const int itmax, int *iter, double *err);
    double (*snrm)(const gsl_vector *sx, const int itol);
} Linbcg;

void linbcg_solve(Linbcg *self, const gsl_vector *b, gsl_vector *x, const int itol, const double tol, const int itmax, int *iter, double *err);
double linbcg_snrm(const gsl_vector *sx, const int itol);

#endif // LINBCG_H
