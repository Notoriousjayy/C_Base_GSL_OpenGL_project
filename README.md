# C_Base_GSL_OpenGL_project

## Usage/Examples

### Graph

#### TODO
```c
#include <string.h>
#include "Graph/graph.h"

int main() {
    // Initialize area for memory allocation
    Area area;
    init_area(area);

    // Create a new graph with 3 vertices
    Graph *g = new_graph(3);

    // Allocate memory for vertex names
    Vertex *v1 = &g->vertices[0];
    Vertex *v2 = &g->vertices[1];
    Vertex *v3 = &g->vertices[2];

    // Assign names to the vertices
    v1->name = save_string("Vertex 1");
    v2->name = save_string("Vertex 2");
    v3->name = save_string("Vertex 3");

    // Add edges between vertices
    new_edge(v1, v2, 10);  // Edge from Vertex 1 to Vertex 2 with length 10
    new_edge(v2, v3, 20);  // Edge from Vertex 2 to Vertex 3 with length 20
    new_edge(v3, v1, 30);  // Edge from Vertex 3 to Vertex 1 with length 30

    // Optionally, set graph ID and utility types
    snprintf(g->id, ID_FIELD_SIZE, "MyGraph");
    strncpy(g->util_types, "example", sizeof(g->util_types) - 1);

    // Perform operations on the graph
    switch_to_graph(g);
    // Example: look up a vertex by name
    Vertex *lookup_vertex = hash_lookup("Vertex 2", g);
    if (lookup_vertex) {
        printf("Found vertex: %s\n", lookup_vertex->name);
    } else {
        printf("Vertex not found\n");
    }

    // Clean up
    recycle(g);
    graph_free(area);

    return 0;
}
```

#### TODO
```c
#include <stdio.h>
#include <stdlib.h>
#include "Graph/graph_basic.h"
#include "Graph/graph.h"

void print_graph(Graph* g, const char* name) {
    if (g == NULL) {
        printf("%s: Graph creation failed.\n", name);
        return;
    }
    printf("%s: Graph with %ld vertices and %ld edges.\n", name, g->n, g->m);
    // Add more detailed printing or visualization of the graph if needed
    // Example: printing the adjacency matrix or list
}

int main() {
    Graph* g;

    // Testing macros for creating specific types of graphs
    g = complete(5);
    print_graph(g, "Complete Graph with 5 vertices");

    g = transitive(5);
    print_graph(g, "Transitive Graph with 5 vertices");

    g = empty(5);
    print_graph(g, "Empty Graph with 5 vertices");

    g = circuit(5);
    print_graph(g, "Circuit Graph with 5 vertices");

    g = cycle(5);
    print_graph(g, "Cycle Graph with 5 vertices");

    // Testing macros for creating graphs using subsets function
    g = disjoint_subsets(5, 2);
    print_graph(g, "Disjoint Subsets Graph (5, 2)");

    g = petersen();
    print_graph(g, "Petersen Graph");

    // Testing macros for creating graphs using perms function
    g = all_perms(3, 1);
    print_graph(g, "All Permutations Graph with 3 vertices (directed)");

    // Testing macros for creating graphs using parts function
    g = all_parts(3, 1);
    print_graph(g, "All Partitions Graph with 3 vertices (directed)");

    // Testing macros for creating graphs using binary function
    g = all_trees(3, 1);
    print_graph(g, "All Trees Graph with 3 vertices (directed)");

    // Testing additional graph creation functions
    g = bi_complete();
    print_graph(g, "Bipartite Complete Graph");

    g = wheel();
    print_graph(g, "Wheel Graph");

    // Free the graph memory if needed (depends on your graph implementation)
    // Example: free_graph(g);

    return 0;
}
```

#### TODO
```c
#include "Graph/graph_io.h"

int main() {
    // Test imap_chr and imap_ord
    char test_char = imap_chr(65);
    printf("imap_chr(65) = %c\n", test_char);

    long test_ord = imap_ord('A');
    printf("imap_ord('A') = %ld\n", test_ord);

    // Test graph_newline
    graph_newline();
    printf("graph_newline() called\n");

    // Test new_checksum
    char* test_str = "Hello, World!";
    long checksum = new_checksum(test_str, 0);
    printf("new_checksum(\"%s\", 0) = %ld\n", test_str, checksum);

    // Test graph_eof
    long eof = graph_eof();
    printf("graph_eof() = %ld\n", eof);

    // Test graph_char
    char next_char = graph_char();
    printf("graph_char() = %c\n", next_char);

    // Test graph_backup
    graph_backup();
    printf("graph_backup() called\n");

    // Test graph_digit
    long digit = graph_digit('5');
    printf("graph_digit('5') = %ld\n", digit);

    // Test graph_number
    unsigned long number = graph_number('5');
    printf("graph_number('5') = %lu\n", number);

    // Test graph_string
    char* result_str = graph_string(str_buf, '!');
    printf("graph_string(str_buf, '!') = %s\n", result_str);

    // Test graph_raw_open and graph_open
    graph_raw_open("testfile.txt");
    printf("graph_raw_open(\"testfile.txt\") called\n");

    long open_result = graph_open("testfile.txt");
    printf("graph_open(\"testfile.txt\") = %ld\n", open_result);

    // Test graph_close
    long close_result = graph_close();
    printf("graph_close() = %ld\n", close_result);

    return 0;
}

```

### Linear Algebra

#### Gaussian elimination (row reduction)
```c
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "Linear-Algebra/lu_decomp.h"

int main() {
    gsl_matrix* a = gsl_matrix_alloc(3, 3);
    gsl_matrix* b = gsl_matrix_alloc(3, 1);
    gsl_vector* x = gsl_vector_alloc(3);

    // Initialize matrices with example values
    gsl_matrix_set(a, 0, 0, 3); gsl_matrix_set(a, 0, 1, 2); gsl_matrix_set(a, 0, 2, -4);
    gsl_matrix_set(a, 1, 0, 2); gsl_matrix_set(a, 1, 1, 3); gsl_matrix_set(a, 1, 2, 3);
    gsl_matrix_set(a, 2, 0, 5); gsl_matrix_set(a, 2, 1, -3); gsl_matrix_set(a, 2, 2, 1);

    gsl_matrix_set(b, 0, 0, 3);
    gsl_matrix_set(b, 1, 0, 15);
    gsl_matrix_set(b, 2, 0, 14);

    LUdcmp* ludcmp = lu_decomp_create(a);

    lu_decomp_solve_matrix(ludcmp, b, b);

    printf("Matrix A after decomposition:\n");
    for (size_t i = 0; i < a->size1; ++i) {
        for (size_t j = 0; j < a->size2; ++j) {
            printf("%f ", gsl_matrix_get(ludcmp->lu, i, j));
        }
        printf("\n");
    }

    printf("Matrix B (solution):\n");
    for (size_t i = 0; i < b->size1; ++i) {
        for (size_t j = 0; j < b->size2; ++j) {
            printf("%f ", gsl_matrix_get(b, i, j));
        }
        printf("\n");
    }

    gsl_matrix_free(a);
    gsl_matrix_free(b);
    gsl_vector_free(x);
    lu_decomp_free(ludcmp);

    return 0;
}
```

#### Lower–Upper (LU) decomposition (Factorization)
```c
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "Linear-Algebra/lu_decomp.h"

int main() {
    gsl_matrix* a = gsl_matrix_alloc(3, 3);
    gsl_matrix* b = gsl_matrix_alloc(3, 1);
    gsl_vector* x = gsl_vector_alloc(3);

    // Initialize matrices with example values
    gsl_matrix_set(a, 0, 0, 3); gsl_matrix_set(a, 0, 1, 2); gsl_matrix_set(a, 0, 2, -4);
    gsl_matrix_set(a, 1, 0, 2); gsl_matrix_set(a, 1, 1, 3); gsl_matrix_set(a, 1, 2, 3);
    gsl_matrix_set(a, 2, 0, 5); gsl_matrix_set(a, 2, 1, -3); gsl_matrix_set(a, 2, 2, 1);

    gsl_matrix_set(b, 0, 0, 3);
    gsl_matrix_set(b, 1, 0, 15);
    gsl_matrix_set(b, 2, 0, 14);

    LUdcmp* ludcmp = lu_decomp_create(a);

    lu_decomp_solve_matrix(ludcmp, b, b);

    printf("Matrix A after decomposition:\n");
    for (size_t i = 0; i < a->size1; ++i) {
        for (size_t j = 0; j < a->size2; ++j) {
            printf("%f ", gsl_matrix_get(ludcmp->lu, i, j));
        }
        printf("\n");
    }

    printf("Matrix B (solution):\n");
    for (size_t i = 0; i < b->size1; ++i) {
        for (size_t j = 0; j < b->size2; ++j) {
            printf("%f ", gsl_matrix_get(b, i, j));
        }
        printf("\n");
    }

    gsl_matrix_free(a);
    gsl_matrix_free(b);
    gsl_vector_free(x);
    lu_decomp_free(ludcmp);

    return 0;
}
```

#### Triangular matrix
```c
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include "Linear-Algebra/tridag.h"

int main() {
    int n = 5;
    gsl_vector* a = gsl_vector_alloc(n);
    gsl_vector* b = gsl_vector_alloc(n);
    gsl_vector* c = gsl_vector_alloc(n);
    gsl_vector* r = gsl_vector_alloc(n);
    gsl_vector* u = gsl_vector_alloc(n);
    gsl_vector* x = gsl_vector_alloc(n);

    // Initialize vectors with example values
    for (int i = 0; i < n; i++) {
        gsl_vector_set(a, i, -1.0);
        gsl_vector_set(b, i, 2.0);
        gsl_vector_set(c, i, -1.0);
        gsl_vector_set(r, i, 1.0);
    }

    // Solve the tridiagonal system
    tridag(a, b, c, r, u);

    printf("Solution vector u:\n");
    for (int i = 0; i < n; i++) {
        printf("%f ", gsl_vector_get(u, i));
    }
    printf("\n");

    // Solve the cyclic system
    cyclic(a, b, c, 1.0, 1.0, r, x);

    printf("Solution vector x:\n");
    for (int i = 0; i < n; i++) {
        printf("%f ", gsl_vector_get(x, i));
    }
    printf("\n");

    gsl_vector_free(a);
    gsl_vector_free(b);
    gsl_vector_free(c);
    gsl_vector_free(r);
    gsl_vector_free(u);
    gsl_vector_free(x);

    return 0;
}
```

#### Banded matrix
```c
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "Linear-Algebra/bandec.h"
#include "Utils/utils.h"

int main() {
    int n = 5;
    int m1 = 1;
    int m2 = 1;
    gsl_matrix *a = gsl_matrix_alloc(n, n);  // Changed to n x n
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector *b = gsl_vector_alloc(n);

    if (a == NULL || x == NULL || b == NULL) {
        fprintf(stderr, "Failed to allocate memory for matrices or vectors.\n");
        if (a != NULL) gsl_matrix_free(a);
        if (x != NULL) gsl_vector_free(x);
        if (b != NULL) gsl_vector_free(b);
        return -1;
    }

    // Initialize the band matrix 'a' with example values
    gsl_matrix_set_zero(a);

    // Band matrix with 1 sub-diagonal and 1 super-diagonal
    // Elements on main diagonal
    gsl_matrix_set(a, 0, 0, 4.0);
    gsl_matrix_set(a, 1, 1, 4.0);
    gsl_matrix_set(a, 2, 2, 4.0);
    gsl_matrix_set(a, 3, 3, 4.0);
    gsl_matrix_set(a, 4, 4, 4.0);

    // Elements on super-diagonal
    gsl_matrix_set(a, 0, 1, 1.0);
    gsl_matrix_set(a, 1, 2, 1.0);
    gsl_matrix_set(a, 2, 3, 1.0);
    gsl_matrix_set(a, 3, 4, 1.0);

    // Elements on sub-diagonal
    gsl_matrix_set(a, 1, 0, 1.0);
    gsl_matrix_set(a, 2, 1, 1.0);
    gsl_matrix_set(a, 3, 2, 1.0);
    gsl_matrix_set(a, 4, 3, 1.0);

    // Initialize the vector 'b' with example values
    for (int i = 0; i < n; i++) {
        gsl_vector_set(b, i, 1.0);
    }

    // Ensure the dimensions of the matrix and vectors are correct
    printf("Matrix dimensions: %lu x %lu\n", a->size1, a->size2);
    printf("Vector dimensions: %lu\n", b->size);

    print_matrix(a, "A");
    print_vector(b, "B");

    // Solve the system
    gsl_vector_memcpy(x, b);
    solve_band_system(a, b, x, m1, m2);
    printf("Solution vector x:\n");
    print_vector(x, "x");

    // Calculate the determinant
    double det = determinant_band_matrix(a, m1, m2);
    printf("Determinant: %f\n", det);

    gsl_matrix_free(a);
    gsl_vector_free(x);
    gsl_vector_free(b);

    return 0;
}
```

#### Singular Value Decomposition (SVD) )
```c
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

int main() {
    const size_t n = 5; // Example size
    gsl_spmatrix *A = gsl_spmatrix_alloc(n, n); // triplet format

    // Populate the matrix A
    gsl_spmatrix_set(A, 0, 0, 10.0);
    gsl_spmatrix_set(A, 1, 1, 20.0);
    gsl_spmatrix_set(A, 2, 2, 30.0);
    gsl_spmatrix_set(A, 3, 3, 40.0);
    gsl_spmatrix_set(A, 4, 4, 50.0);

    // Allocate the destination matrix for the compressed column format
    gsl_spmatrix *C = gsl_spmatrix_alloc_nzmax(n, n, gsl_spmatrix_nnz(A), GSL_SPMATRIX_CSC);

    // Convert to compressed column format
    gsl_spmatrix_csc(C, A);

    // Right hand side vector b
    gsl_vector *b = gsl_vector_alloc(n);
    gsl_vector_set_all(b, 1.0); // Example input vector b

    // Solution vector x
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set_zero(x); // Initial guess for x

    const double tol = 1e-6;  // solution relative tolerance
    const size_t max_iter = 1000; // maximum iterations
    const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, n, 0);

    size_t iter = 0;
    int status;

    do {
        status = gsl_splinalg_itersolve_iterate(C, b, tol, x, work);

        double residual = gsl_splinalg_itersolve_normr(work);
        fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);

        if (status == GSL_SUCCESS) {
            fprintf(stderr, "Converged\n");
            break;
        }
    } while (status == GSL_CONTINUE && ++iter < max_iter);

    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Failed to converge\n");
    }

    printf("Iterations: %zu\n", iter);
    gsl_vector_fprintf(stdout, x, "%g");

    gsl_splinalg_itersolve_free(work);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_spmatrix_free(A);
    gsl_spmatrix_free(C);

    return 0;
}
```

#### Bi-Conjugate Gradient
```c
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

int main() {
    const size_t n = 5; // Example size
    gsl_spmatrix *A = gsl_spmatrix_alloc(n, n); // triplet format

    // Populate the matrix A
    gsl_spmatrix_set(A, 0, 0, 10.0);
    gsl_spmatrix_set(A, 1, 1, 20.0);
    gsl_spmatrix_set(A, 2, 2, 30.0);
    gsl_spmatrix_set(A, 3, 3, 40.0);
    gsl_spmatrix_set(A, 4, 4, 50.0);

    // Allocate the destination matrix for the compressed column format
    gsl_spmatrix *C = gsl_spmatrix_alloc_nzmax(n, n, gsl_spmatrix_nnz(A), GSL_SPMATRIX_CSC);

    // Convert to compressed column format
    gsl_spmatrix_csc(C, A);

    // Right hand side vector b
    gsl_vector *b = gsl_vector_alloc(n);
    gsl_vector_set_all(b, 1.0); // Example input vector b

    // Solution vector x
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector_set_zero(x); // Initial guess for x

    const double tol = 1e-6;  // solution relative tolerance
    const size_t max_iter = 1000; // maximum iterations
    const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, n, 0);

    size_t iter = 0;
    int status;

    do {
        status = gsl_splinalg_itersolve_iterate(C, b, tol, x, work);

        double residual = gsl_splinalg_itersolve_normr(work);
        fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);

        if (status == GSL_SUCCESS) {
            fprintf(stderr, "Converged\n");
            break;
        }
    } while (status == GSL_CONTINUE && ++iter < max_iter);

    if (status != GSL_SUCCESS) {
        fprintf(stderr, "Failed to converge\n");
    }

    printf("Iterations: %zu\n", iter);
    gsl_vector_fprintf(stdout, x, "%g");

    gsl_splinalg_itersolve_free(work);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_spmatrix_free(A);
    gsl_spmatrix_free(C);

    return 0;
}
```

#### Vandermonde matrix
```c
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include "Linear-Algebra/toeplz.h"
int main() {
    size_t n = 4;
    gsl_vector *r = gsl_vector_alloc(2 * n - 1);
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector *y = gsl_vector_alloc(n);

    // Initialize vectors
    gsl_vector_set(r, 0, 4.0);
    gsl_vector_set(r, 1, 3.0);
    gsl_vector_set(r, 2, 2.0);
    gsl_vector_set(r, 3, 1.0);
    gsl_vector_set(r, 4, 1.0);
    gsl_vector_set(r, 5, 2.0);
    gsl_vector_set(r, 6, 3.0);

    gsl_vector_set(y, 0, 1.0);
    gsl_vector_set(y, 1, 2.0);
    gsl_vector_set(y, 2, 3.0);
    gsl_vector_set(y, 3, 4.0);

    toeplz(r, x, y);

    // Print result
    for (size_t i = 0; i < n; i++) {
        printf("x[%zu] = %g\n", i, gsl_vector_get(x, i));
    }

    gsl_vector_free(r);
    gsl_vector_free(x);
    gsl_vector_free(y);

    return 0;
}
```

#### Toeplitz matrix
```c
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include "Linear-Algebra/toeplz.h"
int main() {
    size_t n = 4;
    gsl_vector *r = gsl_vector_alloc(2 * n - 1);
    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector *y = gsl_vector_alloc(n);

    // Initialize vectors
    gsl_vector_set(r, 0, 4.0);
    gsl_vector_set(r, 1, 3.0);
    gsl_vector_set(r, 2, 2.0);
    gsl_vector_set(r, 3, 1.0);
    gsl_vector_set(r, 4, 1.0);
    gsl_vector_set(r, 5, 2.0);
    gsl_vector_set(r, 6, 3.0);

    gsl_vector_set(y, 0, 1.0);
    gsl_vector_set(y, 1, 2.0);
    gsl_vector_set(y, 2, 3.0);
    gsl_vector_set(y, 3, 4.0);

    toeplz(r, x, y);

    // Print result
    for (size_t i = 0; i < n; i++) {
        printf("x[%zu] = %g\n", i, gsl_vector_get(x, i));
    }

    gsl_vector_free(r);
    gsl_vector_free(x);
    gsl_vector_free(y);

    return 0;
}
```

#### Cholesky decomposition
```c
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "Linear-Algebra/cholesky.h"
int main() {
    size_t n = 3;
    gsl_matrix *a = gsl_matrix_alloc(n, n);
    gsl_vector *b = gsl_vector_alloc(n);
    gsl_vector *x = gsl_vector_alloc(n);

    gsl_matrix_set(a, 0, 0, 25.0);
    gsl_matrix_set(a, 0, 1, 15.0);
    gsl_matrix_set(a, 0, 2, -5.0);
    gsl_matrix_set(a, 1, 0, 15.0);
    gsl_matrix_set(a, 1, 1, 18.0);
    gsl_matrix_set(a, 1, 2, 0.0);
    gsl_matrix_set(a, 2, 0, -5.0);
    gsl_matrix_set(a, 2, 1, 0.0);
    gsl_matrix_set(a, 2, 2, 11.0);

    gsl_vector_set(b, 0, 1.0);
    gsl_vector_set(b, 1, 1.0);
    gsl_vector_set(b, 2, 1.0);

    Cholesky ch;
    cholesky_decomp(&ch, a);
    cholesky_solve(&ch, b, x);

    for (size_t i = 0; i < n; i++) {
        printf("x[%zu] = %g\n", i, gsl_vector_get(x, i));
    }

    gsl_matrix_free(a);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_matrix_free(ch.el);

    return 0;
}
```

#### QR decomposition
```c
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
```

### Electronics

#### RL and RC Circuit Proportions
If the resistance in the RL circuit is zero, show that the current I(t) is directly proportional to the integral of the applied voltage E(t). Similarly show that if the resistance in the RC circuit is zero, the current is directly proportional to the derivative of the applied voltage. (In engineering applications, it is often necessary to generate a voltage, rather than a current, which is the integral or derivative of another voltage.)
```c
#include <stdio.h>
#include "circuit.h"

int main() {
    double L = 1.0;  // Inductance
    double C = 1.0;  // Capacitance

    for (double t = 0; t <= 10; t += 0.1) {
        double integral = integrate_E(t);
        double derivative = differentiate_E(t);

        double I_RL = integral / L;    // Current in RL circuit
        double I_RC = C * derivative;  // Current in RC circuit

        printf("t = %4.1f: I_RL = %10.7f, I_RC = %10.7f\n", t, I_RL, I_RC);
    }

    return 0;
}
```
#### Power Laws for Circuits
The power generated or dissipated by a circuit element equals the voltage across the element times the current through the element. Show that the power dissipated by a resistor equals **I^2** R, the power associated with an inductor equals the derivative of **(1/2) LI^2**, and the power associated with a capacitor equals the derivative of (1/2)C(E^2)_E.
```c
#include <stdio.h>
#include "power_calculation.h"

int main() {
    double I = 2.0; // Current in Amperes
    double R = 5.0; // Resistance in Ohms
    double L = 0.1; // Inductance in Henry
    double C = 0.01; // Capacitance in Farads
    double V = 10.0; // Voltage in Volts

    // Calculate power dissipated by a resistor
    double P_resistor = power_resistor(I, R);
    printf("Power dissipated by the resistor: %f W\n", P_resistor);

    // Calculate power associated with an inductor
    double dI_dt = 1.0; // Change in current per second (A/s)
    double P_inductor = power_inductor(L, I, dI_dt);
    printf("Power associated with the inductor: %f W\n", P_inductor);

    // Calculate power associated with a capacitor
    double dV_dt = 2.0; // Change in voltage per second (V/s)
    double P_capacitor = power_capacitor(C, V, dV_dt);
    printf("Power associated with the capacitor: %f W\n", P_capacitor);

    // Using GSL to compute derivatives (optional)
    gsl_function F;
    F.function = &energy_inductor;
    F.params = &L;
    double dE_L_dt = derivative_inductor(&F, I);
    printf("Derivative of inductor energy: %f\n", dE_L_dt);

    F.function = &energy_capacitor;
    F.params = &C;
    double dE_C_dt = derivative_capacitor(&F, V);
    printf("Derivative of capacitor energy: %f\n", dE_C_dt);

    return 0;
}
```

#### Power Balance in Circuits (RL)
Derive a power balance equation for the RL and RC circuits. Discuss the significance of the signs of the three power terms.
### Significance of the Signs of the Power Terms

- **Source Power (\(P(t)\)):**
    - Positive when power is supplied to the circuit.

- **Resistive Power (\(P_R(t)\)):**
    - Always positive since it represents the power dissipated as heat in the resistor.

- **Inductive Power (\(P_L(t)\)):**
    - Can be positive or negative. Positive when energy is being stored in the magnetic field of the inductor and negative when energy is being released from the inductor back to the circuit.

- **Capacitive Power (\(P_C(t)\)):**
    - Can be positive or negative. Positive when energy is being stored in the electric field of the capacitor and negative when energy is being released from the capacitor back to the circuit.

In summary, the signs of the power terms indicate the direction of energy flow: resistors always dissipate energy, while inductors and capacitors can either store or release energy depending on the current and voltage relationships at any given time.
```c
#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "rl_circuit.h"

int main() {
    double R = 1.0;  // Resistance in ohms
    double L = 1.0;  // Inductance in henrys
    double V = 1.0;  // Voltage in volts

    double params[3] = {R, L, V};

    gsl_odeiv2_system sys = {rl_ode, NULL, 1, params};

    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    double t = 0.0;
    double y[1] = {0.0};  // Initial current
    double t1 = 10.0;     // Time to integrate to
    double h = 0.1;       // Step size

    for (double ti = t; ti < t1; ti += h) {
        int status = gsl_odeiv2_driver_apply(d, &ti, ti + h, y);
        if (status != GSL_SUCCESS) {
            printf("Error: %d\n", status);
            break;
        }
        printf("t = %.2f, I = %.5f\n", ti, y[0]);
    }

    gsl_odeiv2_driver_free(d);

    return 0;
}
```

#### Power Balance in Circuits (RC)
Derive a power balance equation for the RL and RC circuits. Discuss the significance of the signs of the three power terms.
### Significance of the Signs of the Power Terms

- **Source Power (\(P(t)\)):**
    - Positive when power is supplied to the circuit.

- **Resistive Power (\(P_R(t)\)):**
    - Always positive since it represents the power dissipated as heat in the resistor.

- **Inductive Power (\(P_L(t)\)):**
    - Can be positive or negative. Positive when energy is being stored in the magnetic field of the inductor and negative when energy is being released from the inductor back to the circuit.

- **Capacitive Power (\(P_C(t)\)):**
    - Can be positive or negative. Positive when energy is being stored in the electric field of the capacitor and negative when energy is being released from the capacitor back to the circuit.

In summary, the signs of the power terms indicate the direction of energy flow: resistors always dissipate energy, while inductors and capacitors can either store or release energy depending on the current and voltage relationships at any given time.
```c
#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "rc_circuit.h"

int main() {
    double R = 1.0;  // Resistance in ohms
    double C = 1.0;  // Capacitance in farads
    double V = 1.0;  // Voltage in volts

    double params[3] = {R, C, V};

    gsl_odeiv2_system sys = {rc_ode, NULL, 1, params};

    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    double t = 0.0;
    double y[1] = {0.0};  // Initial voltage across capacitor
    double t1 = 10.0;     // Time to integrate to
    double h = 0.1;       // Step size

    for (double ti = t; ti < t1; ti += h) {
        int status = gsl_odeiv2_driver_apply(d, &ti, ti + h, y);
        if (status != GSL_SUCCESS) {
            printf("Error: %d\n", status);
            break;
        }
        printf("t = %.2f, Vc = %.5f\n", ti, y[0]);
    }

    gsl_odeiv2_driver_free(d);

    return 0;
}
```

#### Time for RL Circuit
An industrial electromagnet can be modeled as an RL circuit, while it is being energized with a voltage source. If the inductance is 10 H and the wire windings contain 3 Ω of resistance, how long does it take a constant applied voltage to energize the electromagnet to within 90% of its final value (that is, the current equals 90% of its asymptotic value)?
```c
#include <stdio.h>
#include "rl_circuit.h"

int main() {
    double L = 10.0;  // Inductance in Henries
    double R = 3.0;   // Resistance in Ohms

    // Calculate the time to reach 90% of final value
    double time_to_90_percent = time_to_reach_90_percent(L, R);

    // Print the result
    printf("Time to reach 90%% of final value: %f seconds\n", time_to_90_percent);

    return 0;
}
```

#### Capacitor Discharge Time
A 10⁻⁸-F capacitor (10 nanofarads) is charged to 50 V and then disconnected. One can model the charge leakage of the capacitor with an RC circuit with no voltage source and the resistance of the air between the capacitor plates. On a cold dry day, the resistance of the air gap is 5 × 10¹³ Ω; on a humid day, the resistance is 7 × 10⁶ Ω. How long will it take the capacitor voltage to dissipate to half its original value on each day?
```c
#include <stdio.h>
#include "dissipation.h"

int main() {
    double R_cold = 5e13;   // Resistance on a cold dry day (in ohms)
    double R_humid = 7e6;   // Resistance on a humid day (in ohms)
    double C = 10e-9;       // Capacitance (in farads)

    double t_cold = calculate_dissipation_time(R_cold, C);
    double t_humid = calculate_dissipation_time(R_humid, C);

    printf("Time to dissipate to half the original value:\n");
    printf("On a cold dry day: %f seconds (approximately %f hours)\n", t_cold, t_cold / 3600.0);
    printf("On a humid day: %f seconds\n", t_humid);

    return 0;
}
```

#### Rate Change R Calculation
**Electricity** The combined electrical resistance **R** of two resistors **R_1** and **R_2**, connected in parallel, is given by **1/R = 1/R_1 + 1/R_2** where R, R_1, and R_2 are measured in ohms. R_1 and R_2 are increasing at rates of 1 and 1.5 ohms per second, respectively. At what rate is R changing when R_1 = 50 ohms and R_2 = 75 ohms?
```c
#include <stdio.h>
#include "resistance.h"

int main() {
    // Given values
    double R1 = 50.0;  // ohms
    double R2 = 75.0;  // ohms
    double dR1_dt = 1.0;  // ohms per second
    double dR2_dt = 1.5;  // ohms per second

    // Calculate the rate of change of the equivalent resistance
    double dR_dt = calculate_rate_of_change(R1, R2, dR1_dt, dR2_dt);

    // Output the result
    printf("The rate at which R is changing: %f ohms per second\n", dR_dt);

    return 0;
}
```

#### Current Rate Given Voltages
**Electrical Circuit** The voltage V in volts of an electrical circuit is V = IR, where R is the resistance in ohms and I is the current in amperes. R is increasing at a rate of 2 ohms per second, and V is increasing at a rate of 3 volts per second. At what rate is I changing when V = 12 volts and R = 4 ohms?
```c
#include "related_rates.h"

int main() {
    // Given values
    double V = 12.0; // Volts
    double R = 4.0; // Ohms
    double dV_dt = 3.0; // Volts per second
    double dR_dt = 2.0; // Ohms per second

    calculate_current_rate(V, R, dV_dt, dR_dt);

    return 0;
}
```

#### Rate of Parallel Resistance
**Electricity** The combined electrical resistance R of two resistors R_1 and R_2, connected in parallel, is given by

1/R = 1/R_1 + 1/R_2

where R, R_1, and R_2 are measured in ohms. R_1 and R_2 are increasing at rates of 1 and 1.5 ohms per second, respectively. At what rate is R changing when R_1 = 50 ohms and R_2 = 75 ohms?
```c
#include <stdio.h>
#include "resistance.h"

int main() {
    // Given values
    double R1 = 50.0;
    double R2 = 75.0;
    double dR1_dt = 1.0;
    double dR2_dt = 1.5;

    // Calculate R using the given R1 and R2
    double R = 1.0 / (1.0 / R1 + 1.0 / R2);

    // Print R
    printf("R = %f ohms\n", R);

    // Calculate the rate of change of R
    double dR_dt = calculate_dR_dt(R1, R2, dR1_dt, dR2_dt);

    // Print the rate of change of R
    printf("dR/dt = %f ohms/second\n", dR_dt);

    return 0;
}
```

#### Current Calculation L'Hôpital
**Electric Circuit** If voltage V is first applied at time t = 0, then the current I flowing through the circuit at time t is given by I =V/R(1 - e^{-Rt/L) where L is the inductance and R is the resistance. Use L'Hôpital's Rule to find the formula for the current by fixing V and L and letting R approach 0 from the right.
```c
#include <stdio.h>
#include "current.h"

int main() {
    // Define the values of V, t, L, and R
    double V = 5.0;   // Voltage in volts
    double t = 2.0;   // Time in seconds
    double L = 0.01;  // Inductance in henries
    double R = 0.0;   // Resistance approaching 0 (not used in the final formula)

    // Calculate the current using the derived formula
    double current = compute_current(V, t, L, R);

    // Print the result
    printf("The current I as R approaches 0 is: %f A\n", current);

    return 0;
}
```

#### Circuit Current and Limit
**Electric Circuit** A model of the current I, in amperes (A), at time t is given by the first-order differential equation: L dI/dt + RI = E(t) where E(t) is the voltage (V) produced by the power source, R is the resistance, in ohms Ω, and L is the inductance, in henrys (H). Suppose the electric circuit consists of a 24-V power source, a 12-Ω resistor, and a 4-H inductor.

**a** Sketch a slope field for the differential equation.

**b** What is the limiting value of the current?

```c
#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "circuit.h"

int main(void) {
    // Constants for the problem
    double R = 12.0; // Ohms
    double L = 4.0;  // Henrys
    double E = 24.0; // Volts

    // Initial conditions
    double I0 = 0.0; // Initial current
    double t = 0.0;  // Initial time
    double t1 = 10.0; // End time
    double h = 1e-6;  // Initial step size

    // Parameters to pass to the system function
    double params[3] = {R, L, E};

    // Initialize the ODE system
    gsl_odeiv2_system sys = {electric_circuit_func, NULL, 1, &params};

    // Define the step type and control variables
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, h, 1e-8, 1e-8);

    double y[1] = {I0}; // Array to hold the current

    // Print initial condition
    printf("t = %.5f, I(t) = %.5f\n", t, y[0]);

    // Integrate the ODE
    for (int i = 1; i <= 100; i++) {
        double ti = i * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

        if (status != GSL_SUCCESS) {
            printf("Error: %s\n", gsl_strerror(status));
            break;
        }

        printf("t = %.5f, I(t) = %.5f\n", t, y[0]);
    }

    // Free the driver
    gsl_odeiv2_driver_free(d);

    return 0;
}
```

#### Electric Circuit Differential Equation
**Electric Circuit** Find the current I as a function of time t (in seconds), given that I satisfies the differential equation L dI/dt + RI = sin(2t) where R and L are nonzero constants.
```c
#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "electric_circuit.h"

// Define the parameters for the differential equation
double L = 1.0; // Inductance
double R = 1.0; // Resistance

int main() {
    // Define the initial conditions
    double y[1] = {0.0}; // Initial current

    // Combine parameters into an array
    double params[2] = {L, R};

    // Define the ODE system
    gsl_odeiv2_system sys = {electric_circuit_func, NULL, 1, params};

    // Create a stepper and define the stepping function
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    // Solve the ODE system for a range of t values
    double t = 0.0;
    double t1 = 10.0;
    double dt = 0.1;

    // Print the header
    printf("t\t\tI(t)\n");

    for (int i = 0; t < t1; i++) {
        int status = gsl_odeiv2_driver_apply(driver, &t, t + dt, y);

        if (status != GSL_SUCCESS) {
            printf("Error: %s\n", gsl_strerror(status));
            break;
        }

        printf("%lf\t%lf\n", t, y[0]);
    }

    // Free the driver
    gsl_odeiv2_driver_free(driver);

    return 0;
}
```

#### Work Done Moving Electrons
**Electric Force** Two electrons repel each other with a force that varies inversely as the square of the distance between them. One electron is fixed at the point (2, 4). Find the work done in moving the second electron from (-2, 4) to (1, 4).
```c
#include <stdio.h>
#include "electron_work.h"

int main() {
double lower_limit = -2.0;
double upper_limit = 1.0;

double work_done = compute_work(lower_limit, upper_limit);

printf("The work done in moving the second electron is: %f\n", work_done);

return 0;
}
```

#### Electric force fields summary
**Electric force fields are defined by Coulomb's Law, which states that the force exerted on a particle with electric charge q_1 located at (x, y, z) by a particle with electric charge q_2 located at (0, 0, 0) is** F(x, y, z) = (cq_1 q_2 / ||r||^2)u where r = xi + yj + zk, u = r/ ||r||, and c is a constant that depends on the choice of units for ||r||, q_1, and q_2.**

**Note that an electric force field has the same form as a gravitational field. That is,**

F(x, y, z) = (k/ ||r||)^2)u

**Such a force field is called an inverse square field.**

```c
#include <stdio.h>
#include "coulomb_force.h"

int main() {
    // Charges in coulombs
    double q1 = 1.0; // Charge q1
    double q2 = -1.0; // Charge q2

    // Position vector r (x, y, z)
    gsl_vector *r = gsl_vector_alloc(3);
    gsl_vector_set(r, 0, 1.0); // x
    gsl_vector_set(r, 1, 2.0); // y
    gsl_vector_set(r, 2, 3.0); // z

    // Force vector F (output)
    gsl_vector *F = gsl_vector_alloc(3);

    // Calculate the electric force
    calculate_electric_force(q1, q2, r, F);

    // Print the result
    printf("Electric force vector F: (%.5e, %.5e, %.5e)\n",
           gsl_vector_get(F, 0), gsl_vector_get(F, 1), gsl_vector_get(F, 2));

    // Free the allocated memory
    gsl_vector_free(r);
    gsl_vector_free(F);

    return 0;
}
```

#### Equipotential Curves for Voltages
**Electric Potential** The electric potential V at any point (x, y) is V(x, y) = 5/sqrt(25 + x^2 + y^2)
Sketch the equipotential curves for V = 1/2, V = 1/3, and V = 1/4.

```c
#include <stdio.h>
#include <math.h>
#include "equipotential_curves.h"

int main() {
    double radii[] = {5 * sqrt(3), 10 * sqrt(2), 5 * sqrt(15)};
    const char *filenames[] = {"circle_1_2.dat", "circle_1_3.dat", "circle_1_4.dat"};

    double x[NUM_POINTS], y[NUM_POINTS];

    for (int i = 0; i < 3; i++) {
        compute_circle(radii[i], x, y);
        write_data_to_file(filenames[i], x, y);
    }

    // Plotting with Gnuplot
    FILE *gnuplot = popen("gnuplot", "w");
    if (gnuplot == NULL) {
        printf("Error opening Gnuplot!\n");
        return 1;
    }

    fprintf(gnuplot, "set title 'Equipotential Curves'\n");
    fprintf(gnuplot, "set xlabel 'x'\n");
    fprintf(gnuplot, "set ylabel 'y'\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "plot 'circle_1_2.dat' with lines title 'V = 1/2', "
                     "'circle_1_3.dat' with lines title 'V = 1/3', "
                     "'circle_1_4.dat' with lines title 'V = 1/4'\n");
    fflush(gnuplot);

    printf("Press Enter to exit...");
    getchar();

    pclose(gnuplot);

    return 0;
}
```

#### Max Percent Error Calculation
**Power** Electrical power P is given by P = E^2/R where E is voltage and R is resistance. Approximate the maximum percent error in calculating power when 120 volts is applied to a 2000-ohm resistor and the possible percent errors in measuring E and R are 3% and 4%, respectively.

```c
#include <stdio.h>
#include "power_calculation.h"

int main() {
    // Input values
    double E = 120.0;
    double R = 2000.0;
    double E_error_percent = 3.0;
    double R_error_percent = 4.0;

    // Calculate power
    double power = calculate_power(E, R);
    printf("Calculated Power: %f W\n", power);

    // Calculate relative error
    double relative_error = calculate_relative_error(E_error_percent, R_error_percent);
    double percent_error = relative_error * 100.0;
    printf("Maximum Percent Error in Power: %f%%\n", percent_error);

    return 0;
}
```

#### Approximate Inductance with Given Values
**Inductance** The inductance L (in microhenrys) of a straight nonmagnetic wire in free space is  L = 0.00021(ln( 2h/r - 0.75)) where h is the length of the wire in millimeters and r is the radius of a circular cross section. Approximate L when r = 2 +- 1/16 millimeters} and h = 100 +- 1/100 millimeters.

```c
#include <stdio.h>
#include "inductance.h"

int main() {
    // Given values
    double r = 2.0; // millimeters
    double h = 100.0; // millimeters

    // Compute the inductance L
    double L = calculate_inductance(r, h);

    // Print the result
    printf("The inductance L is approximately: %lf microhenrys\n", L);

    return 0;
}
```

#### Average Current Calculations
**Electricity** The oscillating current in an electrical circuit is I = 2 sin(60πt) + cos(120πt) where I is measured in amperes and t is measured in seconds. Find the average current for each time interval.

(a) 0 <= t <= 1/60

(b) 0 <= t <= 1/240

(c) 0 <= t <= 1/30
```c
#include <stdio.h>
#include "current.h"

int main() {
    double a = 0.0;
    double b1 = 1.0 / 60.0;
    double b2 = 1.0 / 240.0;
    double b3 = 1.0 / 30.0;

    double avg_current1 = compute_average_current(a, b1);
    double avg_current2 = compute_average_current(a, b2);
    double avg_current3 = compute_average_current(a, b3);

    if (!isnan(avg_current1)) {
        printf("The average current over the interval [0, 1/60] is: %f\n", avg_current1);
    } else {
        printf("Failed to compute the average current over the interval [0, 1/60]\n");
    }

    if (!isnan(avg_current2)) {
        printf("The average current over the interval [0, 1/240] is: %f\n", avg_current2);
    } else {
        printf("Failed to compute the average current over the interval [0, 1/240]\n");
    }

    if (!isnan(avg_current3)) {
        printf("The average current over the interval [0, 1/30] is: %f\n", avg_current3);
    } else {
        printf("Failed to compute the average current over the interval [0, 1/30]\n");
    }

    return 0;
}
```

#### Magnetic Potential Calculation
**Electromagnetic Theory** The magnetic potential P at a point on the axis of a circular coil is given by P = 2π NIr/k ∫_c^∞ 1/((r^2 + x^2)^3/2) dx where N, I, r, k, and c are constants. Find P.

```c
#include <stdio.h>
#include <math.h>
#include "magnetic_potential.h"

int main() {
    // Define constants
    double N = 100;     // Number of turns
    double I = 5;       // Current in amperes
    double r = 0.1;     // Radius of the coil in meters
    double k = 1;       // Constant (e.g., permeability)
    double c = 0.05;    // Lower limit of the integral

    // Calculate the magnetic potential P
    double P = calculate_magnetic_potential(N, I, r, k, c);

    // Check for NaN indicating an error
    if (isnan(P)) {
        fprintf(stderr, "Failed to calculate the magnetic potential P due to an integration error.\n");
        return 1; // Return an error code
    }

    // Print the result
    printf("The magnetic potential P is: %g\n", P);

    return 0;
}
```

#### Current Solution and Analysis
**Electric Circuits** use the differential equation for electric circuits given by

L dI/dt + RI = E

In this equation, I is the current, R is the resistance, L is the inductance, and E is the electromotive force (voltage).

Solve the differential equation for the current given a constant voltage E_0.

Use the result to find the equation for the current when I(0) = 0, E_0 = 120 volts, R = 600 ohms, and L = 4 henrys. When does the current reach 90% of its limiting value?
```c
#include <stdio.h>
#include "current.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

int main() {
    const double E0 = 120.0;    // Electromotive force in volts
    const double R = 600.0;     // Resistance in ohms

    gsl_odeiv2_system sys = {current_solver_func, NULL, 1, NULL};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    double t = 0.0, t1 = 0.1;
    double y[1] = {0.0}; // Initial condition I(0) = 0

    // Find the limiting value of the current
    double I_lim = E0 / R;

    // Find 90% of the limiting value
    double I_target = 0.9 * I_lim;

    // Integrate until the current reaches 90% of its limiting value
    while (y[0] < I_target) {
        int status = gsl_odeiv2_driver_apply(d, &t, t1, y);
        if (status != GSL_SUCCESS) {
            printf("Error: %s\n", gsl_strerror(status));
            break;
        }
    }

    printf("The current reaches 90%% of its limiting value at t = %.5f seconds\n", t);

    gsl_odeiv2_driver_free(d);
    return 0;
}
```

#### RL Circuit Response
An RL circuit with a 5-Ω resistor and a 0.05-H inductor carries a current of 1 A at t = 0, at which time a voltage source E(t) = 5 cos 120t V is added. Determine the subsequent inductor current and voltage.
```c
#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include "rl_circuit.h"

int main(void) {
    const double L = 0.05;
    const double initial_current = 1.0;
    double t = 0.0, t1 = 0.05;
    double y[2] = { initial_current, 0.0 }; // y[0] = current, y[1] = di/dt

    gsl_odeiv2_system sys = { rl_circuit_func, NULL, 2, NULL };
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

    FILE *fp = fopen("output.dat", "w");
    if (!fp) {
        perror("Unable to open output.dat");
        return 1;
    }

    for (int i = 1; i <= 1000; i++) {
        double ti = i * t1 / 1000.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS) {
            printf("error, return value=%d\n", status);
            break;
        }
        fprintf(fp, "%.5e %.5e %.5e\n", t, y[0], L * y[1]);
    }
    fclose(fp);
    gsl_odeiv2_driver_free(d);

    // Plot using Gnuplot
    FILE *gp = popen("gnuplot -persistent", "w");
    if (gp) {
        fprintf(gp, "set title 'RL Circuit'\n");
        fprintf(gp, "set xlabel 'Time (s)'\n");
        fprintf(gp, "set ylabel 'Current (A) and Voltage (V)'\n");
        fprintf(gp, "plot 'output.dat' using 1:2 with lines title 'Inductor Current', 'output.dat' using 1:3 with lines title 'Inductor Voltage'\n");
        pclose(gp);
    } else {
        perror("Unable to open Gnuplot");
    }

    return 0;
}
```

#### RC Circuit Solution
An RC circuit with a 1-Ω resistor and a 0.000001-F capacitor is driven by a voltage **E(t) = sin(100t) V**. If the initial capacitor voltage is zero, determine the subsequent resistor and capacitor voltages and the current.

```c
#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "Circuits/rc_circuit.h"

#define R 1.0
#define C 0.000001
#define OMEGA 100.0

int main() {
    // Define the ODE system
    gsl_odeiv2_system sys = {rc_circuit_func, jac, 1, NULL};

    // Define the step type and step object
    gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, 1);
    gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(1e-6, 0.0);
    gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(1);

    double t = 0.0, t1 = 1.0; // Time range
    double h = 1e-6; // Initial step size
    double y[1] = {0.0}; // Initial condition for V_C(t)

    // Open a file to save results
    FILE *fp = fopen("rc_circuit_output.txt", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }

    // Time stepping loop
    while (t < t1) {
        int status = gsl_odeiv2_evolve_apply(evolve, control, step, &sys, &t, t1, &h, y);
        if (status != GSL_SUCCESS) {
            fprintf(stderr, "Error: %s\n", gsl_strerror(status));
            break;
        }

        // Calculate current and resistor voltage
        double i_t = y[0] / (R * C);
        double v_r = i_t * R;

        // Write the results to file
        fprintf(fp, "%g %g %g %g\n", t, y[0], v_r, i_t);
    }

    // Free memory
    gsl_odeiv2_evolve_free(evolve);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);

    fclose(fp);

    return 0;
}
```

#### RC Circuit Time Calculation
The voltage source models the transmitting gate, and the capacitor models the receiving gate. Typically, the resistance is 100 Ω, and the capacitance is very small, **10^-12 F** (1 picofarad, pF). If the capacitor is initially uncharged and the transmitting gate changes instantaneously from 0 to 5 V, how long will it take for the voltage at the receiving gate to reach 3 V?
```c
#include <stdio.h>
#include "Circuits/rc_circuit.h"

int main() {
    // Given values
    double R = 100; // Resistance in ohms
    double C = 1e-12; // Capacitance in farads (1 picofarad)
    double V_final = 5; // Final voltage in volts
    double V_desired = 3; // Desired voltage in volts

    // Calculate the time to reach the desired voltage
    double time_in_ps = calculate_time_to_reach_voltage(R, C, V_final, V_desired);

    // Print the result
    printf("Time to reach %.2fV: %.2f picoseconds\n", V_desired, time_in_ps);

    return 0;
}
```

#### Pressure Calculation at Altitude
**Air Pressure** Under ideal conditions, air pressure decreases continuously with the height above sea level at a rate proportional to the pressure at that height. The barometer reads 30 inches at sea level and 15 inches at 18,000 feet. Find the barometric pressure at 35,000 feet.
```c
#include <stdio.h>
#include "Physics/barometric_pressure.h"
#define HEIGHT_35000       35000 // Height at 35000 feet in feet

int main() {
    double k = find_k();
    double pressure_at_35000 = compute_pressure(k, HEIGHT_35000);
    printf("The barometric pressure at 35,000 feet is approximately %.2f inches\n", pressure_at_35000);
    return 0;
}
```

#### Minimum Resistance at T
**Electrical Resistance** The resistance **R** of a certain type of resistor is given by the formula: **R = sqrt(0.001T^4 - 4T + 100) where **R** is measured in ohms and the temperature **T** is measured in degrees Celsius.

**a.** Find **dR/dT** and the critical number of the function. Determine the minimum resistance for this type of resistor.

**b.** Graph the function **R** and use the graph to approximate the minimum resistance for this type of resistor.

```c
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include "Circuits/resistance.h"

int main() {
    gsl_function F;
    double result, abserr;

    // Using GSL to find the derivative at a point
    F.function = &resistance_derivative;
    F.params = NULL;

    double T_critical = 10.0; // Critical point we found analytically
    double R_min = resistance_at_temperature(T_critical, NULL); // Minimum resistance at T = 10

    printf("Critical Temperature (T): %g\n", T_critical);
    printf("Minimum Resistance (R) at T = %g: %g ohms\n", T_critical, R_min);

    // Using GSL to numerically differentiate at the critical point
    gsl_deriv_central(&F, T_critical, 1e-8, &result, &abserr);

    printf("Numerical derivative dR/dT at T = %g: %g\n", T_critical, result);
    printf("Estimated error: %g\n", abserr);

    return 0;
}
```

#### Change in resistance
**Resistance** The total resistance **R** (in ohms) of two resistors connected in parallel is given by **1/R= 1/R1+ 1/R_2**. Approximate the change in **R** as **R_1** is increased from 10 ohms to 10.5 ohms and **R_2** is decreased from 15 ohms to 13 ohms.
```c
#include <stdio.h>
#include "Circuits/resistance.h"

int main() {
    // Initial resistances
    double R1_initial = 10.0;
    double R2_initial = 15.0;

    // New resistances
    double R1_new = 10.5;
    double R2_new = 13.0;

    // Calculate initial total resistance
    double R_initial = calculate_resistance(R1_initial, R2_initial);
    printf("Initial total resistance: %.2f ohms\n", R_initial);

    // Calculate new total resistance
    double R_new = calculate_resistance(R1_new, R2_new);
    printf("New total resistance: %.2f ohms\n", R_new);

    // Calculate change in total resistance
    double delta_R = R_new - R_initial;
    printf("Change in total resistance: %.2f ohms\n", delta_R);

    return 0;
}

```

### Graphics

#### 2D Graphics
```c
#include "Graphics/graphics2d.h"
#include <GL/glut.h>

Canvas canvas;
RGBpixmap pixmap;

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT);

    // Set window and viewport
    setCanvasWindow(&canvas, 0, 500, 0, 500);
    setCanvasViewport(&canvas, 0, 500, 0, 500);

    // Draw a triangle
    canvasMoveTo(&canvas, 100, 100);
    canvasLineTo(&canvas, 200, 100);
    canvasLineTo(&canvas, 150, 200);
    canvasLineTo(&canvas, 100, 100);

    // Draw a BMP image if it was successfully read
    if (pixmap.pixel != NULL) {
        drawRGBpixmap(&pixmap);
    }

    glFlush();
}

int main(int argc, char** argv) {
    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Graphics2D Test");

    // Set up the display function
    glutDisplayFunc(display);

    // Initialize the canvas
    initCanvas(&canvas, 500, 500, "Canvas");

    // Read a BMP file into the pixmap
    if (!readBMPFile(&pixmap, "example.bmp")) {
        printf("Failed to read BMP file\n");
    }

    // Start the GLUT main loop
    glutMainLoop();

    // Free the pixmap memory
    freeRGBpixmap(&pixmap);

    return 0;
}
```

#### 2D Graphics (Main Game Loop)
```c
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "Graphics/graphics2d.h"
#include "Graphics/graphics_app2D.h"


int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        fprintf(stderr, "Failed to initialize GLFW\n");
        return -1;
    }

    // Create a GLFW window
    GLFWwindow* window = glfwCreateWindow(800, 600, "2D Graphics Window", NULL, NULL);
    if (!window) {
        fprintf(stderr, "Failed to create GLFW window\n");
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        return -1;
    }

    // Set the viewport
    glViewport(0, 0, 800, 600);

    // Set the clear color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    // Initialize graphics
    initGraphics();

    // Main game loop
    while (!glfwWindowShouldClose(window)) {
        // Poll for and process events
        glfwPollEvents();

        // Render scene
        render();
    }

    // Cleanup and exit
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
```

#### TODO: NAME
```c
#include <stdio.h>
#include "geometry.h"

int main()
{
    Point P0 = {0, 0, 0};
    Point P1 = {1, 1, 0};
    Point P2 = {2, 2, 0};

    printf("isLeft: %d\n", isLeft(P0, P1, P2));

    Point V0 = {0, 0, 0};
    Point V1 = {1, 0, 0};
    Point V2 = {0, 1, 0};

    printf("orientation2D_Triangle: %d\n", orientation2D_Triangle(V0, V1, V2));
    printf("area2D_Triangle: %f\n", area2D_Triangle(V0, V1, V2));

    Point polygon[] = {
            {0, 0, 0},
            {4, 0, 0},
            {4, 4, 0},
            {0, 4, 0},
            {0, 0, 0}
    };

    printf("orientation2D_Polygon: %d\n", orientation2D_Polygon(4, polygon));
    printf("area2D_Polygon: %f\n", area2D_Polygon(4, polygon));

    Point N = {0, 0, 1};
    printf("area3D_Polygon: %f\n", area3D_Polygon(4, polygon, N));

    return 0;
}
```

#### TODO: NAME
```c
#include <stdio.h>
#include "geometry.h"

int main() {
    Point points[] = { {1.0, 2.0, 0}, {2.0, 3.0, 0}, {3.0, 4.0, 0} };
    Line line = { {0.0, 0.0, 0}, {1.0, 1.0, 0} };
    Segment segment = { {0.0, 0.0, 0}, {1.0, 1.0, 0} };
    int n = sizeof(points) / sizeof(points[0]);

    int closestIndex = closest2D_Point_to_Line(points, n, line);
    printf("Closest point to line is at index: %d\n", closestIndex);

    Point P = {2.0, 2.0, 0};
    float distanceToLine = dist_Point_to_Line(P, line);
    printf("Distance from point to line: %f\n", distanceToLine);

    float distanceToSegment = dist_Point_to_Segment(P, segment);
    printf("Distance from point to segment: %f\n", distanceToSegment);

    return 0;
}
```

#### TODO: NAME
```c
#include <stdio.h>
#include "geometry.h"

int main() {
    // Example usage
    Point polygon[] = {{0, 0}, {5, 0}, {5, 5}, {0, 5}, {0, 0}};
    Point testPoint = {3, 3};

    int cn_result = cn_PnPoly(testPoint, polygon, 4);
    int wn_result = wn_PnPoly(testPoint, polygon, 4);

    printf("Crossing Number Test: %d\n", cn_result);
    printf("Winding Number Test: %d\n", wn_result);

    return 0;
}
```
#### TODO: NAME
```c
#include <stdio.h>
#include "geometry.h"

// Example usage
int main() {
    Point P = {1.0, 2.0, 3.0};
    Plane PL = {{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};
    Point B;

    float distance = dist_Point_to_Plane(P, PL, &B);
    printf("Distance: %f\n", distance);
    printf("Base Point: (%f, %f, %f)\n", B.x, B.y, B.z);

    return 0;
}
```

#### TODO: NAME
```c
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "Graphics/graphics2d.h"

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window);

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        fprintf(stderr, "Failed to initialize GLFW\n");
        return -1;
    }

    // Create a GLFW window
    GLFWwindow* window = glfwCreateWindow(800, 600, "OpenGL Window", NULL, NULL);
    if (!window) {
        fprintf(stderr, "Failed to create GLFW window\n");
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        return -1;
    }

    // Set the viewport
    glViewport(0, 0, 800, 600);

    // Set the clear color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    // Initialize Canvas
    Canvas canvas;
    initCanvas(&canvas, 800, 600, "OpenGL Window");

    // Set the viewport and projection
    setCanvasViewport(&canvas, 0, 800, 0, 600);
    setCanvasWindow(&canvas, -1.0f, 1.0f, -1.0f, 1.0f);

    // Set the framebuffer size callback
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Main game loop
    while (!glfwWindowShouldClose(window)) {
        // Process input
        processInput(window);

        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT);

        // Drawing commands using the Canvas functions
        drawLine(0.0, 0.0, 1.0, 1.0);
        drawCircle(0.5, 0.5, 0.3);
        drawFilledCircle(-0.5, -0.5, 0.3);
        drawRectangle(-0.5, 0.5, 0.4, 0.2);
        drawFilledRectangle(0.5, -0.5, 0.4, 0.2);

        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Cleanup and exit
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}

void processInput(GLFWwindow *window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, 1);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}
```

#### Circle
```c
#include <stdio.h>
#include "geometry.h"
#include "Utils/utils.h"

int main() {
    // Create two circles
    Circle c1 = { {{0.0f, 0.0f}, {1.0f, 1.0f}}, 5.0f };
    Circle c2 = { {{0.0f, 0.0f}, {1.0f, 1.0f}}, 5.0f };

    // Print circles
    printCircle(&c1);
    printCircle(&c2);

    // Compare circles
    if (compareCircles(&c1, &c2)) {
        printf("The circles are equal.\n");
    } else {
        printf("The circles are not equal.\n");
    }

    return 0;
}
```

#### Rectangle
```c
#include "geometry.h"
#include <stdio.h>

void printRectangle2D(const Rectangle2D* rect) {
    printf("Rectangle Origin: (%.2f, %.2f), Size: (%.2f, %.2f)\n",
           rect->origin.x, rect->origin.y, rect->size.x, rect->size.y);
}

int main() {
    Rectangle2D rect1;
    initRectangle2D(&rect1);

    Point origin = {2.0f, 3.0f};
    Vector size = {4.0f, 5.0f};
    Rectangle2D rect2;
    initRectangle2DWithValues(&rect2, &origin, &size);

    printRectangle2D(&rect1);
    printRectangle2D(&rect2);

    return 0;
}
```
#### Oriented Rectangle
```c
#include "geometry.h"
#include <stdio.h>

void printOrientedRectangle(const OrientedRectangle* rect) {
printf("Oriented Rectangle Position: (%.2f, %.2f), Half Extents: (%.2f, %.2f), Rotation: %.2f\n",
rect->position.x, rect->position.y, rect->halfExtents.x, rect->halfExtents.y, rect->rotation);
}

int main() {
OrientedRectangle rect1;
initOrientedRectangle(&rect1);

Point position = {2.0f, 3.0f};
Vector halfExtents = {4.0f, 5.0f};
float rotation = 45.0f;
OrientedRectangle rect2;
initOrientedRectangleWithParams(&rect2, &position, &halfExtents, rotation);

printOrientedRectangle(&rect1);
printOrientedRectangle(&rect2);

return 0;
}
```

#### Point in Circle
```c
#include "geometry.h"
#include <stdio.h>

int main() {
Point point = {1.0f, 1.0f};
Circle circle = {{0.0f, 0.0f}, 2.0f};
Rectangle2D rectangle = {{0.0f, 0.0f}, {3.0f, 3.0f}};
OrientedRectangle orientedRectangle = {{0.0f, 0.0f}, {3.0f, 3.0f}, 0.0f};

if (PointInCircle(&point, &circle)) {
printf("Point is inside the circle.\n");
} else {
printf("Point is outside the circle.\n");
}

if (PointInRectangle(&point, &rectangle)) {
printf("Point is inside the rectangle.\n");
} else {
printf("Point is outside the rectangle.\n");
}

if (PointInOrientedRectangle(&point, &orientedRectangle)) {
printf("Point is inside the oriented rectangle.\n");
} else {
printf("Point is outside the oriented rectangle.\n");
}

return 0;
}
```

#### Circle to circle collision
```c
#include "geometry.h"
#include <stdio.h>

int main() {
Circle circle1 = {{0.0f, 0.0f}, 2.0f};
Circle circle2 = {{3.0f, 0.0f}, 2.0f};

if (CircleCircle(&circle1, &circle2)) {
printf("The circles intersect.\n");
} else {
printf("The circles do not intersect.\n");
}

return 0;
}
```

#### Circle to Rectangle collision
```c
#include "geometry.h"
#include <stdio.h>

int main() {
Circle circle = {{{2.0f, 3.0f, 0.0f}, {2.0f, 3.0f, 0.0f}}, 1.0f};
Rectangle2D rect = {{0.0f, 0.0f, 0.0f}, {4.0f, 4.0f, 0.0f}};

if (CircleRectangle(&circle, &rect)) {
printf("The circle intersects with the rectangle.\n");
} else {
printf("The circle does not intersect with the rectangle.\n");
}

return 0;
}
```

#### Circle to Oriented Rectangle collision
```c
#include "geometry.h"
#include <stdio.h>

int main() {
Circle circle = {{{2.0f, 3.0f, 0.0f}, {2.0f, 3.0f, 0.0f}}, 1.0f};
OrientedRectangle rect = {{0.0f, 0.0f, 0.0f}, {4.0f, 4.0f, 0.0f}, 45.0f};

if (CircleOrientedRectangle(&circle, &rect)) {
printf("The circle intersects with the oriented rectangle.\n");
} else {
printf("The circle does not intersect with the oriented rectangle.\n");
}

return 0;
}
```

#### Rectangle to Rectangle collision
```c
#include "geometry.h"
#include <stdio.h>

int main() {
    Rectangle2D rect1 = {{0.0f, 0.0f, 0.0f}, {4.0f, 4.0f, 0.0f}};
    Rectangle2D rect2 = {{2.0f, 2.0f, 0.0f}, {4.0f, 4.0f, 0.0f}};

    if (RectangleRectangle(&rect1, &rect2)) {
        printf("The rectangles overlap.\n");
    } else {
        printf("The rectangles do not overlap.\n");
    }

    return 0;
}
```

#### Rectangle to Oriented Rectangle collision
```c
#include "geometry.h"
#include <stdio.h>

int main() {
Rectangle2D rect1 = {{0.0f, 0.0f, 0.0f}, {4.0f, 4.0f, 0.0f}};
OrientedRectangle rect2 = {{2.0f, 2.0f, 0.0f}, {2.0f, 2.0f, 0.0f}, 45.0f};

if (RectangleOrientedRectangle(&rect1, &rect2)) {
printf("The rectangles overlap.\n");
} else {
printf("The rectangles do not overlap.\n");
}

return 0;
}
```

#### Rectangle to Oriented Rectangle collision
```c
#include "geometry.h"
#include <stdio.h>

int main() {
OrientedRectangle rect1 = {{0.0f, 0.0f, 0.0f}, {2.0f, 2.0f, 0.0f}, 30.0f};
OrientedRectangle rect2 = {{3.0f, 3.0f, 0.0f}, {2.0f, 2.0f, 0.0f}, 45.0f};

if (OrientedRectangleOrientedRectangle(&rect1, &rect2)) {
printf("The oriented rectangles overlap.\n");
} else {
printf("The oriented rectangles do not overlap.\n");
}

return 0;
}
```

#### TODO
```c
#include "Geometry/geometry.h"
#include <stdio.h>

int main() {
    Vector v1 = {1.0f, 2.0f};
    Vector v2 = {3.0f, 4.0f};
    Vector v3 = {1.0f, 2.0f, 3.0f};
    Vector v4 = {4.0f, 5.0f, 6.0f};

    Vector v2_add = Vector_add(&v1, &v2);
    Vector v3_add = Vector_add(&v3, &v4);

    printf("Vector add: (%f, %f)\n", v2_add.x, v2_add.y);
    printf("Vector add: (%f, %f, %f)\n", v3_add.x, v3_add.y, v3_add.z);

    Vector v2_sub = Vector_sub(&v1, &v2);
    Vector v3_sub = Vector_sub(&v3, &v4);

    printf("Vector sub: (%f, %f)\n", v2_sub.x, v2_sub.y);
    printf("Vector sub: (%f, %f, %f)\n", v3_sub.x, v3_sub.y, v3_sub.z);

    Vector v2_mul = Vector_mul(&v1, &v2);
    Vector v3_mul = Vector_mul(&v3, &v4);

    printf("Vector mul: (%f, %f)\n", v2_mul.x, v2_mul.y);
    printf("Vector mul: (%f, %f, %f)\n", v3_mul.x, v3_mul.y, v3_mul.z);

    Vector v2_mul_scalar = Vector_mul_scalar(&v1, 2.0f);
    Vector v3_mul_scalar = Vector_mul_scalar(&v3, 2.0f);

    printf("Vector mul scalar: (%f, %f)\n", v2_mul_scalar.x, v2_mul_scalar.y);
    printf("Vector mul scalar: (%f, %f, %f)\n", v3_mul_scalar.x, v3_mul_scalar.y, v3_mul_scalar.z);

    printf("Vector eq: %d\n", Vector_eq(&v1, &v2));
    printf("Vector eq: %d\n", Vector_eq(&v3, &v4));

    printf("Vector neq: %d\n", Vector_neq(&v1, &v2));
    printf("Vector neq: %d\n", Vector_neq(&v3, &v4));

    return 0;
}
```

#### TODO
```c
#include "Geometry/geometry.h"
#include <stdio.h>

int main() {
Vector v1 = {1.0f, 2.0f};
Vector v2 = {3.0f, 4.0f};
Vector v3 = {1.0f, 2.0f, 3.0f};
Vector v4 = {4.0f, 5.0f, 6.0f};

Vector v2_add = Vector_add(&v1, &v2);
Vector v3_add = Vector_add(&v3, &v4);

printf("Vector add: (%f, %f)\n", v2_add.x, v2_add.y);
printf("Vector add: (%f, %f, %f)\n", v3_add.x, v3_add.y, v3_add.z);

Vector v2_sub = Vector_sub(&v1, &v2);
Vector v3_sub = Vector_sub(&v3, &v4);

printf("Vector sub: (%f, %f)\n", v2_sub.x, v2_sub.y);
printf("Vector sub: (%f, %f, %f)\n", v3_sub.x, v3_sub.y, v3_sub.z);

Vector v2_mul = Vector_mul(&v1, &v2);
Vector v3_mul = Vector_mul(&v3, &v4);

printf("Vector mul: (%f, %f)\n", v2_mul.x, v2_mul.y);
printf("Vector mul: (%f, %f, %f)\n", v3_mul.x, v3_mul.y, v3_mul.z);

Vector v2_mul_scalar = Vector_mul_scalar(&v1, 2.0f);
Vector v3_mul_scalar = Vector_mul_scalar(&v3, 2.0f);

printf("Vector mul scalar: (%f, %f)\n", v2_mul_scalar.x, v2_mul_scalar.y);
printf("Vector mul scalar: (%f, %f, %f)\n", v3_mul_scalar.x, v3_mul_scalar.y, v3_mul_scalar.z);

printf("Vector eq: %d\n", Vector_eq(&v1, &v2));
printf("Vector eq: %d\n", Vector_eq(&v3, &v4));

printf("Vector neq: %d\n", Vector_neq(&v1, &v2));
printf("Vector neq: %d\n", Vector_neq(&v3, &v4));

float dot_v2 = dot_Vector(&v1, &v2);
float dot_v3 = dot_Vector(&v3, &v4);

printf("dot Vector: %f\n", dot_v2);
printf("dot Vector: %f\n", dot_v3);

float mag_v2 = magnitude_Vector(&v1);
float mag_v3 = magnitude_Vector(&v3);

printf("magnitude Vector: %f\n", mag_v2);
printf("magnitude Vector: %f\n", mag_v3);

float dist_v3 = distance_Vector(&v3, &v4);

printf("distance Vector: %f\n", dist_v3);

Vector norm_v2 = normalized_Vector(&v1);
Vector norm_v3 = normalized_Vector(&v3);

printf("normalized Vector: (%f, %f)\n", norm_v2.x, norm_v2.y);
printf("normalized Vector: (%f, %f, %f)\n", norm_v3.x, norm_v3.y, norm_v3.z);

Vector cross_v3 = cross_Vector(&v3, &v4);

printf("cross Vector: (%f, %f, %f)\n", cross_v3.x, cross_v3.y, cross_v3.z);

float angle_v2 = angle_Vector(&v1, &v2);
float angle_v3 = angle_Vector(&v3, &v4);

printf("angle Vector: %f\n", angle_v2);
printf("angle Vector: %f\n", angle_v3);

Vector proj_v2 = project_Vector(&v1, &v2);
Vector proj_v3 = project_Vector(&v3, &v4);

printf("project Vector: (%f, %f)\n", proj_v2.x, proj_v2.y);
printf("project Vector: (%f, %f, %f)\n", proj_v3.x, proj_v3.y, proj_v3.z);

Vector perp_v2 = perpendicular_Vector(&v1, &v2);
Vector perp_v3 = perpendicular_Vector(&v3, &v4);

printf("perpendicular Vector: (%f, %f)\n", perp_v2.x, perp_v2.y);
printf("perpendicular Vector: (%f, %f, %f)\n", perp_v3.x, perp_v3.y, perp_v3.z);

Vector refl_v2 = reflection_Vector(&v1, &v2);
Vector refl_v3 = reflection_Vector(&v3, &v4);

printf("reflection Vector: (%f, %f)\n", refl_v2.x, refl_v2.y);
printf("reflection Vector: (%f, %f, %f)\n", refl_v3.x, refl_v3.y, refl_v3.z);

return 0;
}
```

#### TODO
```c
#include <stdio.h>
#include "Geometry/geometry.h"
#include "Utils/utils.h"


int main() {
    // Define test matrices
    mat2 matA2 = {1.0f, 2.0f, 3.0f, 4.0f};
    mat2 matB2 = {5.0f, 6.0f, 7.0f, 8.0f};
    mat3 matA3 = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f};
    mat3 matB3 = {9.0f, 8.0f, 7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f};
    mat4 matA4 = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 16.0f};
    mat4 matB4 = {16.0f, 15.0f, 14.0f, 13.0f, 12.0f, 11.0f, 10.0f, 9.0f, 8.0f, 7.0f, 6.0f, 5.0f, 4.0f, 3.0f, 2.0f, 1.0f};

    // Print original matrices
    printf("Original mat2 A:\n");
    printMat2(&matA2);
    printf("Original mat3 A:\n");
    printMat3(&matA3);
    printf("Original mat4 A:\n");
    printMat4(&matA4);

    // Transpose matrices
    mat2 transposed2 = Transpose2(&matA2);
    mat3 transposed3 = Transpose3(&matA3);
    mat4 transposed4 = Transpose4(&matA4);

    printf("\nTransposed mat2 A:\n");
    printMat2(&transposed2);
    printf("\nTransposed mat3 A:\n");
    printMat3(&transposed3);
    printf("\nTransposed mat4 A:\n");
    printMat4(&transposed4);

    // Multiply matrices
    mat2 multiplied2 = mat2_mult(&matA2, &matB2);
    mat3 multiplied3 = mat3_mult(&matA3, &matB3);
    mat4 multiplied4 = mat4_mult(&matA4, &matB4);

    printf("\nMultiplied mat2 A * B:\n");
    printMat2(&multiplied2);
    printf("\nMultiplied mat3 A * B:\n");
    printMat3(&multiplied3);
    printf("\nMultiplied mat4 A * B:\n");
    printMat4(&multiplied4);

    // Determinant
    float det2 = Determinant_mat2(&matA2);
    float det3 = Determinant_mat3(&matA3);
    float det4 = Determinant_mat4(&matA4);

    printf("\nDeterminant of mat2 A: %f\n", det2);
    printf("Determinant of mat3 A: %f\n", det3);
    printf("Determinant of mat4 A: %f\n", det4);

    // Inverse
    mat2 inverse2 = Inverse_mat2(&matA2);
    mat3 inverse3 = Inverse_mat3(&matA3);
    mat4 inverse4 = Inverse_mat4(&matA4);

    printf("\nInverse of mat2 A:\n");
    printMat2(&inverse2);
    printf("\nInverse of mat3 A:\n");
    printMat3(&inverse3);
    printf("\nInverse of mat4 A:\n");
    printMat4(&inverse4);

    return 0;
}

```