//
// Created by jorda on 7/19/2024.
//

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include "../include/magnetic_potential.h"

// Define the integrand function
double integrand(double x, void *params) {
    double r = *(double *)params;
    return 1.0 / pow(r * r + x * x, 1.5);
}

// Function to calculate the magnetic potential P
double calculate_magnetic_potential(double N, double I, double r, double k, double c) {
    // Set up the integration workspace
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

    // Define the integral parameters
    double result, error;
    gsl_function F;
    F.function = &integrand;
    F.params = &r;

    // Perform the integration
    int status = gsl_integration_qagi(&F, 0, 1e-7, 1000, workspace, &result, &error);

    // Check for errors
    if (status != GSL_SUCCESS) {
        fprintf(stderr, "GSL integration error: %s\n", gsl_strerror(status));
        gsl_integration_workspace_free(workspace);
        return NAN; // Return NaN to indicate an error
    }

    // Free the integration workspace
    gsl_integration_workspace_free(workspace);

    // Calculate the magnetic potential P
    double P = (2 * M_PI * N * I * r / k) * result;

    return P;
}
