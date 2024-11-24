//
// Created by jorda on 7/18/2024.
//

#ifndef CURRENT_H
#define CURRENT_H


#include <gsl/gsl_integration.h>

// Function to compute the current
double compute_current(double V, double t, double L, double R);

double compute_average_current(double a, double b);

int current_solver_func(double t, const double y[], double f[], void *params);
#endif // CURRENT_H
