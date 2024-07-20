//
// Created by jorda on 7/19/2024.
//

#include "rc_circuit.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

int rc_ode(double t, const double y[], double f[], void *params) {
    double R = ((double *)params)[0];
    double C = ((double *)params)[1];
    double V = ((double *)params)[2];

    f[0] = (V - y[0] / C) / R;  // dy/dt = (V - Vc) / R

    return GSL_SUCCESS;
}
