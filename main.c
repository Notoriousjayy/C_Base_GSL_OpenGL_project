#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include "Linear-Algebra/gaussj.h"

int main() {
    // Define matrices 'a' and 'b'
    gsl_matrix* a = gsl_matrix_alloc(3, 3);
    gsl_matrix* b = gsl_matrix_alloc(3, 1);

    // Initialize matrices with example values
    // a = [3 2 -4; 2 3 3; 5 -3 1]
    gsl_matrix_set(a, 0, 0, 3); gsl_matrix_set(a, 0, 1, 2); gsl_matrix_set(a, 0, 2, -4);
    gsl_matrix_set(a, 1, 0, 2); gsl_matrix_set(a, 1, 1, 3); gsl_matrix_set(a, 1, 2, 3);
    gsl_matrix_set(a, 2, 0, 5); gsl_matrix_set(a, 2, 1, -3); gsl_matrix_set(a, 2, 2, 1);

    // b = [3; 15; 14]
    gsl_matrix_set(b, 0, 0, 3);
    gsl_matrix_set(b, 1, 0, 15);
    gsl_matrix_set(b, 2, 0, 14);

    // Perform Gauss-Jordan elimination
    gaussj(a, b);

    // Print the result
    printf("Matrix A:\n");
    for (size_t i = 0; i < a->size1; ++i) {
        for (size_t j = 0; j < a->size2; ++j) {
            printf("%f ", gsl_matrix_get(a, i, j));
        }
        printf("\n");
    }

    printf("Matrix B:\n");
    for (size_t i = 0; i < b->size1; ++i) {
        for (size_t j = 0; j < b->size2; ++j) {
            printf("%f ", gsl_matrix_get(b, i, j));
        }
        printf("\n");
    }

    gsl_matrix_free(a);
    gsl_matrix_free(b);

    return 0;
}
