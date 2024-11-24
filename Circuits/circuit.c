//
// Created by jorda on 7/18/2024.
//

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include "../include/circuit.h"

// Define the system of differential equations
int func(double t, const double y[], double f[], void *params) {
    double R = *(double *)params;
    double L = *((double *)params + 1);
    double E = *((double *)params + 2);
    f[0] = (E - R * y[0]) / L;
    return GSL_SUCCESS;
}

// Define the applied voltage function E(t)
double E(double t, void *params) {
    // Example voltage function: E(t) = sin(t)
    return sin(t);
}

// Function to integrate E(t) for the RL circuit
double integrate_E(double t) {
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &E;
    F.params = NULL;

    double result, error;
    gsl_integration_qag(&F, 0, t, 0, 1e-7, 1000, 6, workspace, &result, &error);
    gsl_integration_workspace_free(workspace);

    return result;
}

// Function to differentiate E(t) for the RC circuit
double differentiate_E(double t) {
    gsl_function F;
    F.function = &E;
    F.params = NULL;

    double result, abserr;
    gsl_deriv_central(&F, t, 1e-8, &result, &abserr);

    return result;
}