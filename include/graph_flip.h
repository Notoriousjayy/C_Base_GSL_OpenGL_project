//
// Created by jorda on 7/23/2024.
//

#ifndef GB_FLIP_H
#define GB_FLIP_H

// Macro to get the next random number
#define graph_next_rand() ((*graph_fptr >= 0) ? (*graph_fptr--) : graph_flip_cycle())

// External declarations
extern long* graph_fptr;
extern long graph_flip_cycle(void);

// Function declarations
void graph_init_rand(long seed);
long graph_unif_rand(long m);

#endif // GB_FLIP_H
