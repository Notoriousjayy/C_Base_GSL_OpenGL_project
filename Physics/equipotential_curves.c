//
// Created by jorda on 7/19/2024.
//

#include <stdio.h>
#include <math.h>
#include "equipotential_curves.h"

void compute_circle(double radius, double x[], double y[]) {
    for (int i = 0; i < NUM_POINTS; i++) {
        double theta = 2 * M_PI * i / NUM_POINTS;
        x[i] = radius * cos(theta);
        y[i] = radius * sin(theta);
    }
}

void write_data_to_file(const char *filename, double x[], double y[]) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file %s!\n", filename);
        return;
    }
    for (int i = 0; i < NUM_POINTS; i++) {
        fprintf(file, "%f %f\n", x[i], y[i]);
    }
    fclose(file);
}
