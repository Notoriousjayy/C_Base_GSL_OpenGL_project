//
// Created by jorda on 7/15/2024.
//

#ifndef BANDEC_H
#define BANDEC_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void banmul(const gsl_matrix *a, int m1, int m2, const gsl_vector *x, gsl_vector *b);
void solve_band_system(const gsl_matrix *a, const gsl_vector *b, gsl_vector *x, int m1, int m2);
double determinant_band_matrix(const gsl_matrix *a, int m1, int m2);

#endif
