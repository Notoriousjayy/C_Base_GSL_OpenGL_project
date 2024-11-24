//
// Created by jorda on 7/19/2024.
//

#include <math.h>
#include "../include/rl_circuit.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_log.h>

int rl_circuit_func(double t, const double y[], double f[], void *params) {
    (void)(t); // avoid unused parameter warning
    double R = 5.0;
    double L = 0.05;
    double E = 5 * cos(120 * t);
    f[0] = y[1];
    f[1] = (E - R * y[0]) / L;
    return GSL_SUCCESS;
}

int rl_ode(double t, const double y[], double f[], void *params) {
    double R = ((double *)params)[0];
    double L = ((double *)params)[1];
    double V = ((double *)params)[2];

    f[0] = (V - R * y[0]) / L;  // dy/dt = (V - IR) / L

    return GSL_SUCCESS;
}

double calculate_time_constant(double L, double R) {
    return L / R;
}

double time_to_reach_90_percent(double L, double R) {
    double tau = calculate_time_constant(L, R);
    double ln_0_1 = gsl_sf_log(0.1);
    return -tau * ln_0_1;
}