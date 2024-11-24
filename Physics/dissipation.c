//
// Created by jorda on 7/19/2024.
//

#include <math.h>
#include "../include/dissipation.h"

#define LN2 0.6931471805599453  // Natural logarithm of 2

double calculate_dissipation_time(double R, double C) {
    return R * C * LN2;
}
