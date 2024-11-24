//
// Created by jorda on 7/15/2024.
//

#ifndef GAUSSJ_H
#define GAUSSJ_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void gaussj(gsl_matrix* a, gsl_matrix* b);
void gaussj_single(gsl_matrix* a);

#endif // GAUSSJ_H
