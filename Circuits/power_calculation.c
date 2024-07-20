//
// Created by jorda on 7/19/2024.
//

#include "power_calculation.h"
#include <stdio.h>

double calculate_power(double E, double R) {
    return (E * E) / R;
}

double calculate_relative_error(double E_error_percent, double R_error_percent) {
    return 2 * (E_error_percent / 100.0) + (R_error_percent / 100.0);
}

// Function to calculate power dissipated by a resistor
double power_resistor(double I, double R) {
    return I * I * R;
}

// Function to calculate power associated with an inductor
double energy_inductor(double L, double I) {
    return 0.5 * L * I * I;
}

double power_inductor(double L, double I, double dI_dt) {
    return L * I * dI_dt;
}

// Function to calculate power associated with a capacitor
double energy_capacitor(double C, double V) {
    return 0.5 * C * V * V;
}

double power_capacitor(double C, double V, double dV_dt) {
    return C * V * dV_dt;
}

// Derivative function for GSL
double derivative_inductor(const gsl_function *F, double x) {
    double result, abserr;
    gsl_deriv_central(F, x, 1e-8, &result, &abserr);
    return result;
}

double derivative_capacitor(const gsl_function *F, double x) {
    double result, abserr;
    gsl_deriv_central(F, x, 1e-8, &result, &abserr);
    return result;
}
