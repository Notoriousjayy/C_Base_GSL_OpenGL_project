//
// Created by jorda on 7/19/2024.
//

#ifndef RC_CIRCUIT_H
#define RC_CIRCUIT_H

// Function declarations
int rc_circuit_func(double t, const double y[], double f[], void *params);
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params);
int rc_ode(double t, const double y[], double f[], void *params);
double calculate_time_to_reach_voltage(double R, double C, double V_final, double V_desired);

#endif // RC_CIRCUIT_H

