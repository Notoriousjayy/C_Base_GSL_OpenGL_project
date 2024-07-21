//
// Created by jorda on 7/19/2024.
//

#include "rc_circuit.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <gsl/gsl_sf_log.h>

int rc_ode(double t, const double y[], double f[], void *params) {
    double R = ((double *)params)[0];
    double C = ((double *)params)[1];
    double V = ((double *)params)[2];

    f[0] = (V - y[0] / C) / R;  // dy/dt = (V - Vc) / R

    return GSL_SUCCESS;
}

double calculate_time_to_reach_voltage(double R, double C, double V_final, double V_desired) {
    double ratio = V_desired / V_final;
    double ln_ratio = gsl_sf_log(ratio);
    double RC = R * C;
    double time = -RC * ln_ratio;
    return time * 1e12; // Return time in picoseconds
}

#define R 1.0
#define C 0.000001
#define OMEGA 100.0

// Define the system of ODEs
int rc_circuit_func(double t, const double y[], double f[], void *params) {
    f[0] = (sin(OMEGA * t) - y[0]) / (R * C);
    return GSL_SUCCESS;
}

// Jacobian matrix is not needed for this simple case
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
    return GSL_SUCCESS;
}


