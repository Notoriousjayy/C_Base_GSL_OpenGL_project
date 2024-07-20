//
// Created by jorda on 7/19/2024.
//

#ifndef POWER_CALCULATION_H
#define POWER_CALCULATION_H

#include <gsl/gsl_deriv.h>

double calculate_power(double E, double R);
double calculate_relative_error(double E_error_percent, double R_error_percent);

// Function to calculate power dissipated by a resistor
double power_resistor(double I, double R);

// Function to calculate power associated with an inductor
double energy_inductor(double L, double I);
double power_inductor(double L, double I, double dI_dt);

// Function to calculate power associated with a capacitor
double energy_capacitor(double C, double V);
double power_capacitor(double C, double V, double dV_dt);

// Derivative function for GSL
double derivative_inductor(const gsl_function *F, double x);
double derivative_capacitor(const gsl_function *F, double x);

#endif // POWER_CALCULATION_H
