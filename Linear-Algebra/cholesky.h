//
// Created by jorda on 7/16/2024.
//

#ifndef CHOLESKY_H
#define CHOLESKY_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct {
    int n;
    gsl_matrix *el;
} Cholesky;

/**
 * @brief Performs Cholesky decomposition on the given matrix.
 *
 * @param ch Pointer to the Cholesky structure to be initialized.
 * @param a Input matrix to be decomposed.
 */
void cholesky_decomp(Cholesky *ch, const gsl_matrix *a);

/**
 * @brief Solves the linear system Ax = b using the Cholesky decomposition.
 *
 * @param ch Pointer to the Cholesky structure.
 * @param b Right-hand side vector.
 * @param x Solution vector to be computed.
 */
void cholesky_solve(const Cholesky *ch, const gsl_vector *b, gsl_vector *x);

/**
 * @brief Performs element-wise multiplication using the Cholesky decomposition matrix.
 *
 * @param ch Pointer to the Cholesky structure.
 * @param y Input vector for multiplication.
 * @param b Resultant vector after multiplication.
 */
void cholesky_elmult(const Cholesky *ch, const gsl_vector *y, gsl_vector *b);

/**
 * @brief Solves the lower triangular system using the Cholesky decomposition matrix.
 *
 * @param ch Pointer to the Cholesky structure.
 * @param b Right-hand side vector.
 * @param y Solution vector to be computed.
 */
void cholesky_elsolve(const Cholesky *ch, const gsl_vector *b, gsl_vector *y);

/**
 * @brief Computes the inverse of the matrix using the Cholesky decomposition.
 *
 * @param ch Pointer to the Cholesky structure.
 * @param ainv Output matrix for the inverse.
 */
void cholesky_inverse(const Cholesky *ch, gsl_matrix *ainv);

/**
 * @brief Computes the log determinant of the matrix using the Cholesky decomposition.
 *
 * @param ch Pointer to the Cholesky structure.
 * @return The log determinant of the matrix.
 */
double cholesky_logdet(const Cholesky *ch);

#endif // CHOLESKY_H
