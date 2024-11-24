//
// Created by jorda on 7/19/2024.
//

#include "../include/inductance.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>

double calculate_inductance(double r, double h) {
    double constant = 0.00021;
    double log_value = gsl_sf_log(2 * h / r);
    double L = constant * (log_value - 0.75);
    return L;
}
