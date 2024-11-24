//
// Created by jorda on 7/23/2024.
//

#include "../include/graph_flip.h"

#define graph_next_rand() (*graph_fptr>=0?*graph_fptr--:graph_flip_cycle() )  \

#define mod_diff(x,y) (((x) -(y) ) &0x7fffffff)  \

#define two_to_the_31 ((unsigned long) 0x80000000)  \
                                                    \
// Static array initialization
static long A[56] = {-1};

// Pointer to the random number array
long* graph_fptr = A;

// Function definitions

long graph_flip_cycle(void) {
    long* ii;
    long* jj;

    for (ii = &A[1], jj = &A[32]; jj <= &A[55]; ii++, jj++) {
        *ii = mod_diff(*ii, *jj);
    }
    for (jj = &A[1]; ii <= &A[55]; ii++, jj++) {
        *ii = mod_diff(*ii, *jj);
    }
    graph_fptr = &A[54];
    return A[55];
}

void graph_init_rand(long seed) {
    long i;
    long prev = seed;
    long next = 1;
    seed = prev = mod_diff(prev, 0);
    A[55] = prev;

    for (i = 21; i; i = (i + 21) % 55) {
        A[i] = next;
        next = mod_diff(prev, next);
        if (seed & 1) {
            seed = 0x40000000 + (seed >> 1);
        } else {
            seed >>= 1;
        }
        next = mod_diff(next, seed);
        prev = A[i];
    }

    for (i = 0; i < 5; i++) {
        (void)graph_flip_cycle();
    }
}

long graph_unif_rand(long m) {
    unsigned long t = two_to_the_31 - (two_to_the_31 % m);
    long r;

    do {
        r = graph_next_rand();
    } while (t <= (unsigned long)r);

    return r % m;
}
