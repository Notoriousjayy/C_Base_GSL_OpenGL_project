//
// Created by jorda on 7/18/2024.
//

#ifndef CIRCUIT_H
#define CIRCUIT_H

// Function declaration for the ODE system
int electric_circuit_func(double t, const double y[], double f[], void *params);
double E(double t, void *params);
double integrate_E(double t);
double differentiate_E(double t);

#endif // CIRCUIT_H
