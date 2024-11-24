//
// Created by jorda on 7/15/2024.
//
#include "../include/tridag.h"
#include <stdio.h>
#include <stdlib.h>

void tridag(const gsl_vector* a, const gsl_vector* b, const gsl_vector* c, const gsl_vector* r, gsl_vector* u) {
    int n = a->size;
    gsl_vector* gam = gsl_vector_alloc(n);
    double bet;

    if (gsl_vector_get(b, 0) == 0.0) {
        fprintf(stderr, "Error 1 in tridag\n");
        exit(EXIT_FAILURE);
    }

    gsl_vector_set(u, 0, gsl_vector_get(r, 0) / (bet = gsl_vector_get(b, 0)));

    for (int j = 1; j < n; j++) {
        gsl_vector_set(gam, j, gsl_vector_get(c, j - 1) / bet);
        bet = gsl_vector_get(b, j) - gsl_vector_get(a, j) * gsl_vector_get(gam, j);
        if (bet == 0.0) {
            fprintf(stderr, "Error 2 in tridag\n");
            exit(EXIT_FAILURE);
        }
        gsl_vector_set(u, j, (gsl_vector_get(r, j) - gsl_vector_get(a, j) * gsl_vector_get(u, j - 1)) / bet);
    }

    for (int j = n - 2; j >= 0; j--) {
        gsl_vector_set(u, j, gsl_vector_get(u, j) - gsl_vector_get(gam, j + 1) * gsl_vector_get(u, j + 1));
    }

    gsl_vector_free(gam);
}

void cyclic(const gsl_vector* a, const gsl_vector* b, const gsl_vector* c, double alpha, double beta, const gsl_vector* r, gsl_vector* x) {
    int n = a->size;
    gsl_vector* bb = gsl_vector_alloc(n);
    gsl_vector* u = gsl_vector_alloc(n);
    gsl_vector* z = gsl_vector_alloc(n);
    double gamma = -gsl_vector_get(b, 0);
    double fact;

    if (n <= 2) {
        fprintf(stderr, "n too small in cyclic\n");
        exit(EXIT_FAILURE);
    }

    gsl_vector_memcpy(bb, b);
    gsl_vector_set(bb, 0, gsl_vector_get(b, 0) - gamma);
    gsl_vector_set(bb, n - 1, gsl_vector_get(b, n - 1) - alpha * beta / gamma);

    tridag(a, bb, c, r, x);

    gsl_vector_set(u, 0, gamma);
    gsl_vector_set(u, n - 1, alpha);
    for (int i = 1; i < n - 1; i++) {
        gsl_vector_set(u, i, 0.0);
    }

    tridag(a, bb, c, u, z);

    fact = (gsl_vector_get(x, 0) + beta * gsl_vector_get(x, n - 1) / gamma) / (1.0 + gsl_vector_get(z, 0) + beta * gsl_vector_get(z, n - 1) / gamma);

    for (int i = 0; i < n; i++) {
        gsl_vector_set(x, i, gsl_vector_get(x, i) - fact * gsl_vector_get(z, i));
    }

    gsl_vector_free(bb);
    gsl_vector_free(u);
    gsl_vector_free(z);
}
