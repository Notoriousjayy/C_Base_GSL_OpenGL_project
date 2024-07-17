//
// Created by jorda on 7/16/2024.
//

#include <stdio.h>
#include <stdlib.h>
#include "toeplz.h"

void toeplz(const gsl_vector *r, gsl_vector *x, const gsl_vector *y) {
    int j, k, m, m1, m2, n1, n = y->size;
    double pp, pt1, pt2, qq, qt1, qt2, sd, sgd, sgn, shn, sxn;
    n1 = n - 1;

    if (gsl_vector_get(r, n1) == 0.0) {
        fprintf(stderr, "toeplz-1 singular principal minor\n");
        exit(EXIT_FAILURE);
    }

    gsl_vector_set(x, 0, gsl_vector_get(y, 0) / gsl_vector_get(r, n1));
    if (n1 == 0) return;

    gsl_vector *g = gsl_vector_calloc(n1);
    gsl_vector *h = gsl_vector_calloc(n1);

    gsl_vector_set(g, 0, gsl_vector_get(r, n1 - 1) / gsl_vector_get(r, n1));
    gsl_vector_set(h, 0, gsl_vector_get(r, n1 + 1) / gsl_vector_get(r, n1));

    for (m = 0; m < n; m++) {
        m1 = m + 1;
        sxn = -gsl_vector_get(y, m1);
        sd = -gsl_vector_get(r, n1);
        for (j = 0; j < m + 1; j++) {
            sxn += gsl_vector_get(r, n1 + m1 - j) * gsl_vector_get(x, j);
            sd += gsl_vector_get(r, n1 + m1 - j) * gsl_vector_get(g, m - j);
        }

        if (sd == 0.0) {
            fprintf(stderr, "toeplz-2 singular principal minor\n");
            exit(EXIT_FAILURE);
        }

        gsl_vector_set(x, m1, sxn / sd);

        for (j = 0; j < m + 1; j++) {
            gsl_vector_set(x, j, gsl_vector_get(x, j) - gsl_vector_get(x, m1) * gsl_vector_get(g, m - j));
        }

        if (m1 == n1) {
            gsl_vector_free(g);
            gsl_vector_free(h);
            return;
        }

        sgn = -gsl_vector_get(r, n1 - m1 - 1);
        shn = -gsl_vector_get(r, n1 + m1 + 1);
        sgd = -gsl_vector_get(r, n1);
        for (j = 0; j < m + 1; j++) {
            sgn += gsl_vector_get(r, n1 + j - m1) * gsl_vector_get(g, j);
            shn += gsl_vector_get(r, n1 + m1 - j) * gsl_vector_get(h, j);
            sgd += gsl_vector_get(r, n1 + j - m1) * gsl_vector_get(h, m - j);
        }

        if (sgd == 0.0) {
            fprintf(stderr, "toeplz-3 singular principal minor\n");
            exit(EXIT_FAILURE);
        }

        gsl_vector_set(g, m1, sgn / sgd);
        gsl_vector_set(h, m1, shn / sd);
        k = m;
        m2 = (m + 2) >> 1;
        pp = gsl_vector_get(g, m1);
        qq = gsl_vector_get(h, m1);
        for (j = 0; j < m2; j++) {
            pt1 = gsl_vector_get(g, j);
            pt2 = gsl_vector_get(g, k);
            qt1 = gsl_vector_get(h, j);
            qt2 = gsl_vector_get(h, k);
            gsl_vector_set(g, j, pt1 - pp * qt2);
            gsl_vector_set(g, k, pt2 - pp * qt1);
            gsl_vector_set(h, j, qt1 - qq * pt2);
            gsl_vector_set(h, k--, qt2 - qq * pt1);
        }
    }

    fprintf(stderr, "toeplz - should not arrive here!\n");
    exit(EXIT_FAILURE);
}