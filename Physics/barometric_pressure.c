//
// Created by jorda on 7/17/2024.
//


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "barometric_pressure.h"

#define SEA_LEVEL_PRESSURE 30.0  // Pressure at sea level in inches
#define PRESSURE_AT_18000  15.0  // Pressure at 18000 feet in inches
#define HEIGHT_18000       18000 // Height at 18000 feet in feet
#define HEIGHT_35000       35000 // Height at 35000 feet in feet

// Function to compute the right-hand side of the ODE
int barometric_fuc(double t, const double y[], double f[], void *params) {
    double k = *(double *) params;
    f[0] = -k * y[0];
    return GSL_SUCCESS;
}

// Function to find the proportionality constant k
double find_k() {
    double k;
    k = log(2) / HEIGHT_18000;
    return k;
}

// Function to compute the pressure at a given height
double compute_pressure(double k, double height) {
    gsl_odeiv2_system sys = {barometric_fuc, NULL, 1, &k};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    double t = 0.0, y[1] = {SEA_LEVEL_PRESSURE};
    gsl_odeiv2_driver_apply(d, &t, height, y);
    gsl_odeiv2_driver_free(d);

    return y[0];
}