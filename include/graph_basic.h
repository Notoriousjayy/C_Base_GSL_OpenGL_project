//
// Created by jorda on 7/17/2024.
//

#ifndef GRAPHBASIC_H
#define GRAPHBASIC_H

#include "graph.h"

/* Function declarations for creating various types of graphs */
extern Graph* board();
extern Graph* simplex();
extern Graph* subsets();
extern Graph* perms();
extern Graph* parts();
extern Graph* binary();

extern Graph* complement();
extern Graph* gunion();
extern Graph* intersection();
extern Graph* lines();
extern Graph* product();
extern Graph* induced();

/* Macros for creating specific types of graphs using the board function */
#define complete(n)      board((long)(n), 0L, 0L, 0L, -1L, 0L, 0L)
#define transitive(n)    board((long)(n), 0L, 0L, 0L, -1L, 0L, 1L)
#define empty(n)         board((long)(n), 0L, 0L, 0L, 2L, 0L, 0L)
#define circuit(n)       board((long)(n), 0L, 0L, 0L, 1L, 1L, 0L)
#define cycle(n)         board((long)(n), 0L, 0L, 0L, 1L, 1L, 1L)

/* Macros for creating specific types of graphs using the subsets function */
#define disjoint_subsets(n, k) subsets((long)(k), 1L, (long)(1 - (n)), 0L, 0L, 0L, 1L, 0L)
#define petersen()              disjoint_subsets(5, 2)

/* Macros for creating specific types of graphs using the perms function */
#define all_perms(n, directed)  perms((long)(1 - (n)), 0L, 0L, 0L, 0L, 0L, (long)(directed))

/* Macros for creating specific types of graphs using the parts function */
#define all_parts(n, directed)  parts((long)(n), 0L, 0L, (long)(directed))

/* Macros for creating specific types of graphs using the binary function */
#define all_trees(n, directed)  binary((long)(n), 0L, (long)(directed))

/* Constants for product types */
#define cartesian 0
#define direct    1
#define strong    2

/* Macros for accessing utility fields in vertex structures */
#define ind       z.I    // Index
#define IND_GRAPH 1000000000
#define subst     y.G    // Subgraph

/* Additional function declarations */
extern Graph* bi_complete();
extern Graph* wheel();

#endif // GRAPHBASIC_H
