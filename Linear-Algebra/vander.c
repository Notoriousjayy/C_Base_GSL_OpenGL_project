//
// Created by jorda on 7/16/2024.
//

#include <stdio.h>
#include "vander.h"

void vander(const gsl_vector *x, gsl_vector *w, const gsl_vector *q) {
    size_t i, j, k, n = q->size;
    double b, s, t, xx;
    gsl_vector *c = gsl_vector_calloc(n);

    if (n == 1) {
        gsl_vector_set(w, 0, gsl_vector_get(q, 0));
    } else {
        gsl_vector_set(c, n-1, -gsl_vector_get(x, 0));
        for (i = 1; i < n; i++) {
            xx = -gsl_vector_get(x, i);
            for (j = n-1-i; j < n-1; j++) {
                gsl_vector_set(c, j, gsl_vector_get(c, j) + xx * gsl_vector_get(c, j+1));
            }
            gsl_vector_set(c, n-1, gsl_vector_get(c, n-1) + xx);
        }

        for (i = 0; i < n; i++) {
            xx = gsl_vector_get(x, i);
            t = b = 1.0;
            s = gsl_vector_get(q, n-1);
            for (k = n-1; k > 0; k--) {
                b = gsl_vector_get(c, k) + xx * b;
                s += gsl_vector_get(q, k-1) * b;
                t = xx * t + b;
            }
            gsl_vector_set(w, i, s / t);
        }
    }

    gsl_vector_free(c);
}