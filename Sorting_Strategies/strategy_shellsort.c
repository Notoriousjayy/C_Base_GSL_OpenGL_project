//
// Created by jorda on 9/20/2024.
//

// strategy_shellsort.c
#include "../include/strategy_shellsort.h"
#include "../include/shellsort.h"

// Shellsort strategy function
void strategy_shellsort(int arr[], int aux[], int n) {
    // Shellsort needs an increment array, so we define it here.
    // Example: Using increments based on Shell's original sequence
    int inc[] = {5, 3, 1};
    int t = sizeof(inc) / sizeof(inc[0]);

    // Perform Shellsort using the provided array
    shellsort(arr, n, inc, t);
}
