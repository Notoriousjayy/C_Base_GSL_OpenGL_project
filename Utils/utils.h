//
// Created by jorda on 7/15/2024.
//

#ifndef UTILS_H
#define UTILS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "../Geometry/geometry.h"

void print_matrix(const gsl_matrix *m, const char *name);
void print_vector(const gsl_vector *v, const char *name);
int compareCircles(const Circle* c1, const Circle* c2);
void printCircle(const Circle* c);
void printMat2(const mat2* m);
void printMat3(const mat3* m);
void printMat4(const mat4* m);

#endif // UTILS_H
