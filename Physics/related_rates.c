#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include "../include/related_rates.h"

void calculate_current_rate(double V, double R, double dV_dt, double dR_dt) {
    // Calculate current I
    double I = V / R; // Amperes

    // Setting up the matrix equation A * x = b
    gsl_matrix *A = gsl_matrix_alloc(2, 2);
    gsl_vector *x = gsl_vector_alloc(2);
    gsl_vector *b = gsl_vector_alloc(2);

    // Matrix A
    gsl_matrix_set(A, 0, 0, R);
    gsl_matrix_set(A, 0, 1, 1);
    gsl_matrix_set(A, 1, 0, I);
    gsl_matrix_set(A, 1, 1, 0);

    // Vector b
    gsl_vector_set(b, 0, dV_dt);
    gsl_vector_set(b, 1, dR_dt);

    // Solve for x
    int s;
    gsl_permutation *p = gsl_permutation_alloc(2);

    gsl_linalg_LU_decomp(A, p, &s);
    gsl_linalg_LU_solve(A, p, b, x);

    double dI_dt = gsl_vector_get(x, 1);

    // Output the result
    printf("The rate at which the current I is changing is: %f amperes per second\n", dI_dt);

    // Free the allocated memory
    gsl_matrix_free(A);
    gsl_vector_free(x);
    gsl_vector_free(b);
    gsl_permutation_free(p);
}
