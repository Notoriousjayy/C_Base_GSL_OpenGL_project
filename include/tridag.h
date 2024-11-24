//
// Created by jorda on 7/15/2024.
//

#ifndef TRIDAG_H
#define TRIDAG_H

#include <gsl/gsl_vector.h>

void tridag(const gsl_vector* a, const gsl_vector* b, const gsl_vector* c, const gsl_vector* r, gsl_vector* u);
void cyclic(const gsl_vector* a, const gsl_vector* b, const gsl_vector* c, double alpha, double beta, const gsl_vector* r, gsl_vector* x);

#endif // TRIDAG_H
