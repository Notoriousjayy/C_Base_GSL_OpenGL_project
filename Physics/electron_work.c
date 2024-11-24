//
// Created by jorda on 7/19/2024.
//

// electron_work.c
#include <gsl/gsl_integration.h>
#include "../include/electron_work.h"

#define K_CONSTANT 1.0  // Define the constant k

double force_function(double x, void *params) {
    double u = 2.0 - x;
    return K_CONSTANT / (u * u);
}

double compute_work(double lower_limit, double upper_limit) {
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

    double result, error;

    gsl_function F;
    F.function = &force_function;
    F.params = NULL;

    gsl_integration_qag(&F, lower_limit, upper_limit, 0, 1e-7, 1000, 6, workspace, &result, &error);

    gsl_integration_workspace_free(workspace);
    return result;
}
