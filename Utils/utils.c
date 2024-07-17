//
// Created by jorda on 7/15/2024.
//
#include "utils.h"
#include <stdio.h>

void print_matrix(const gsl_matrix *m, const char *name) {
    printf("Matrix %s:\n", name);
    for (size_t i = 0; i < m->size1; ++i) {
        for (size_t j = 0; j < m->size2; ++j) {
            printf("%f ", gsl_matrix_get(m, i, j));
        }
        printf("\n");
    }
}

void print_vector(const gsl_vector *v, const char *name) {
    printf("Vector %s:\n", name);
    for (size_t i = 0; i < v->size; ++i) {
        printf("%f ", gsl_vector_get(v, i));
    }
    printf("\n");
}
