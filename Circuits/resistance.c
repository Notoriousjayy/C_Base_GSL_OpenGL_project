//
// Created by jorda on 7/18/2024.
//

#include <gsl/gsl_deriv.h>
#include "resistance.h"

double calculate_parallel_resistance(double R1, double R2) {
    return 1.0 / (1.0 / R1 + 1.0 / R2);
}

double calculate_rate_of_change(double R1, double R2, double dR1_dt, double dR2_dt) {
    double R = calculate_parallel_resistance(R1, R2);
    return R * R * (dR1_dt / (R1 * R1) + dR2_dt / (R2 * R2));
}

double calculate_resistance(double R1, double R2) {
    return 1.0 / (1.0 / R1 + 1.0 / R2);
}

double resistance_function(double R1, void *params) {
    double R2 = *(double *)params;
    return 1.0 / R1 + 1.0 / R2;
}

double calculate_dR_dt(double R1, double R2, double dR1_dt, double dR2_dt) {
    // Calculate R using the given R1 and R2
    double R = 1.0 / (1.0 / R1 + 1.0 / R2);

    // Calculate the derivatives directly
    double dR1_inv_dt = -1.0 / (R1 * R1) * dR1_dt;
    double dR2_inv_dt = -1.0 / (R2 * R2) * dR2_dt;

    // Calculate d(1/R)/dt
    double dR_inv_dt = dR1_inv_dt + dR2_inv_dt;

    // Calculate dR/dt
    double dR_dt = -R * R * dR_inv_dt;

    return dR_dt;
}

// Function to compute R(T)
double resistance_at_temperature(double T, void *params) {
    (void)(params); // Unused parameter
    return sqrt(0.001 * pow(T, 4) - 4 * T + 100);
}

// Function to compute the derivative dR/dT
double resistance_derivative(double T, void *params) {
    (void)(params); // Unused parameter
    double numerator = 0.004 * pow(T, 3) - 4;
    double denominator = 2 * sqrt(0.001 * pow(T, 4) - 4 * T + 100);
    return numerator / denominator;
}