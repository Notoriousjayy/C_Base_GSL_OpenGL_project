//
// Created by jorda on 7/17/2024.
//

#ifndef BAROMETRIC_PRESSURE_H
#define BAROMETRIC_PRESSURE_H

double find_k();
double compute_pressure(double k, double height);
int barometric_fuc(double t, const double y[], double f[], void *params);
#endif // BAROMETRIC_PRESSURE_H
