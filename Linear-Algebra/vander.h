//
// Created by jorda on 7/16/2024.
//

#ifndef VANDER_H
#define VANDER_H

#include <gsl/gsl_vector.h>

// Function prototype for the vander function
void vander(const gsl_vector *x, gsl_vector *w, const gsl_vector *q);

#endif // VANDER_H
