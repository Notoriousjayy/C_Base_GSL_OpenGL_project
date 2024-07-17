//
// Created by jorda on 7/15/2024.
//

#ifndef UTILS_H
#define UTILS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void print_matrix(const gsl_matrix *m, const char *name);
void print_vector(const gsl_vector *v, const char *name);

#endif // UTILS_H
