//
// Created by jorda on 7/18/2024.
//

#ifndef RESISTANCE_H
#define RESISTANCE_H

double calculate_parallel_resistance(double R1, double R2);
double calculate_rate_of_change(double R1, double R2, double dR1_dt, double dR2_dt);

// Function to calculate 1/R from R1 and R2
double resistance_function(double R1, void *params);

// Function to calculate the rate of change of the combined resistance R
double calculate_dR_dt(double R1, double R2, double dR1_dt, double dR2_dt);

#endif // RESISTANCE_H
