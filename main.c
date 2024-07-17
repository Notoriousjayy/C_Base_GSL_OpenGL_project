#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "Linear-Algebra/qrdcmp.h"

int main() {
    // Define the matrix A (3x3 matrix for example)
    double data[] = {
            12, -51, 4,
            6, 167, -68,
            -4, 24, -41
    };

    // Define the vector b
    double b_data[] = {1, 2, 3};

    // Define the vectors u and v for update (example vectors)
    double u_data[] = {1, 0, 0};
    double v_data[] = {0, 1, 0};

    // Create GSL matrices and vectors
    gsl_matrix_view A = gsl_matrix_view_array(data, 3, 3);
    gsl_vector_view b = gsl_vector_view_array(b_data, 3);
    gsl_vector_view u = gsl_vector_view_array(u_data, 3);
    gsl_vector_view v = gsl_vector_view_array(v_data, 3);
    gsl_vector *x = gsl_vector_alloc(3);

    // Perform QR decomposition
    QRdcmp *qr = qrdcmp_alloc(&A.matrix);

    // Solve the system Ax = b
    qrdcmp_solve(qr, &b.vector, x);

    // Print the solution x
    printf("Solution x:\n");
    gsl_vector_fprintf(stdout, x, "%g");

    // Perform an update of the QR decomposition
    qrdcmp_update(qr, &u.vector, &v.vector);

    // Solve the updated system
    qrdcmp_solve(qr, &b.vector, x);

    // Print the updated solution x
    printf("Updated solution x:\n");
    gsl_vector_fprintf(stdout, x, "%g");

    // Free the allocated memory
    qrdcmp_free(qr);
    gsl_vector_free(x);

    return 0;
}
