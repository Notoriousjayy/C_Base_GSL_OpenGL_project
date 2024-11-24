//
// Created by jorda on 9/20/2024.
//

// strategy_radix_sort.c
#include "../include/strategy_radix_sort.h"
#include "../include/radix_sort.h"

// Radix Sort strategy function
void strategy_radix_sort(int arr[], int aux[], int n) {
    radix_sort(arr, n, sizeof(int));  // Perform Radix Sort for 32-bit integers
}
