#include "electric_circuit.h"
#include <math.h>
#include <gsl/gsl_errno.h>

// Define the differential equation system
int electric_circuit_func(double t, const double y[], double f[], void *params) {
    double *params_array = (double *)params;
    double L = params_array[0];
    double R = params_array[1];
    f[0] = (1.0 / L) * (sin(2 * t) - R * y[0]);
    return GSL_SUCCESS;
}
