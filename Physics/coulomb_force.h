//
// Created by jorda on 7/19/2024.
//

#ifndef COULOMB_FORCE_H
#define COULOMB_FORCE_H

#include <gsl/gsl_vector.h>

// Define the Coulomb constant (N m^2/C^2)
#define COULOMB_CONSTANT 8.9875517873681764e9

// Function to calculate the electric force vector
void calculate_electric_force(double q1, double q2, gsl_vector *r, gsl_vector *F);

#endif // COULOMB_FORCE_H
