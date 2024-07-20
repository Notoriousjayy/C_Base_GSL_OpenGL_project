//
// Created by jorda on 7/18/2024.
//

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "current.h"

// Define the parameters
const double E0 = 120.0;    // Electromotive force in volts
const double R = 600.0;     // Resistance in ohms
const double L = 4.0;       // Inductance in henrys

// Function to compute the current
double compute_current(double V, double t, double L, double R) {
    // Current I given by the derived formula
    return (V * t) / L;
}

double current(double t, void *params) {
    return 2 * sin(60 * M_PI * t) + cos(120 * M_PI * t);
}

double compute_average_current(double a, double b) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double result, error;

    gsl_function F;
    F.function = &current;
    F.params = NULL;

    // Increase the error tolerance and the maximum number of subdivisions
    double epsabs = 1e-7; // Absolute error tolerance
    double epsrel = 1e-7; // Relative error tolerance
    size_t limit = 1000;  // Maximum number of intervals

    // Perform the integration using a different rule, e.g., GSL_INTEG_GAUSS61
    int status = gsl_integration_qag(&F, a, b, epsabs, epsrel, limit, GSL_INTEG_GAUSS61, w, &result, &error);

    if (status != GSL_SUCCESS) {
        printf("GSL integration error: %s\n", gsl_strerror(status));
        gsl_integration_workspace_free(w);
        return NAN; // Return NaN to indicate an error
    }

    gsl_integration_workspace_free(w);

    return result / (b - a);

}

// Define the differential equation dI/dt = (E0/L) - (R/L) * I
int current_solver_func(double t, const double y[], double f[], void *params) {
    (void)(t); // Avoid unused parameter warning
    f[0] = (E0 / L) - (R / L) * y[0];
    return GSL_SUCCESS;
}