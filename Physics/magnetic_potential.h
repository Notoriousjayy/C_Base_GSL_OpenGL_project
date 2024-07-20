//
// Created by jorda on 7/19/2024.
//

#ifndef MAGNETIC_POTENTIAL_H
#define MAGNETIC_POTENTIAL_H

#include <gsl/gsl_integration.h>

// Function to calculate the integrand
double integrand(double x, void *params);

// Function to calculate the magnetic potential P
double calculate_magnetic_potential(double N, double I, double r, double k, double c);

#endif // MAGNETIC_POTENTIAL_H
