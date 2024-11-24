//
// Created by jorda on 7/19/2024.
//

#include "../include/coulomb_force.h"
#include <gsl/gsl_blas.h>
#include <math.h>

void calculate_electric_force(double q1, double q2, gsl_vector *r, gsl_vector *F) {
    double r_norm = gsl_blas_dnrm2(r); // Calculate the norm of the vector r
    double scale = COULOMB_CONSTANT * q1 * q2 / (r_norm * r_norm * r_norm); // Calculate the scale factor

    // Scale the vector r to get the force vector F
    gsl_vector_memcpy(F, r);
    gsl_vector_scale(F, scale);
}
