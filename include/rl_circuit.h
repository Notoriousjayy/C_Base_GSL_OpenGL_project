//
// Created by jorda on 7/19/2024.
//

#ifndef RL_CIRCUIT_H
#define RL_CIRCUIT_H

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

// Function to define the differential equations
int rl_circuit_func(double t, const double y[], double f[], void *params);

int rl_ode(double t, const double y[], double f[], void *params);

// Function to calculate the time constant of an RL circuit
double calculate_time_constant(double L, double R);

// Function to calculate the time to reach 90% of the final current value
double time_to_reach_90_percent(double L, double R);

#endif // RL_CIRCUIT_H
