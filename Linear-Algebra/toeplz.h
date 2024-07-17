//
// Created by jorda on 7/16/2024.
//

#ifndef TOEPLZ_H
#define TOEPLZ_H

#include <gsl/gsl_vector.h>

// Function to solve Toeplitz system of linear equations
void toeplz(const gsl_vector *r, gsl_vector *x, const gsl_vector *y);

#endif // TOEPLZ_H
