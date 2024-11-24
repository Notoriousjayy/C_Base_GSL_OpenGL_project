//
// Created by jorda on 7/15/2024.
//
#include "../include/utils.h"
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

void printMat2(const mat2* m) {
    printf("mat2:\n");
    printf("[%f, %f]\n", m->_11, m->_12);
    printf("[%f, %f]\n", m->_21, m->_22);
}

void printMat3(const mat3* m) {
    printf("mat3:\n");
    printf("[%f, %f, %f]\n", m->_11, m->_12, m->_13);
    printf("[%f, %f, %f]\n", m->_21, m->_22, m->_23);
    printf("[%f, %f, %f]\n", m->_31, m->_32, m->_33);
}

void printMat4(const mat4* m) {
    printf("mat4:\n");
    printf("[%f, %f, %f, %f]\n", m->_11, m->_12, m->_13, m->_14);
    printf("[%f, %f, %f, %f]\n", m->_21, m->_22, m->_23, m->_24);
    printf("[%f, %f, %f, %f]\n", m->_31, m->_32, m->_33, m->_34);
    printf("[%f, %f, %f, %f]\n", m->_41, m->_42, m->_43, m->_44);
}