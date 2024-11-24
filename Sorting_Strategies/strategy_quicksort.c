//
// Created by jorda on 9/20/2024.
//

// strategy_quicksort.c
#include "../include/strategy_quicksort.h"
#include "../include/quicksort.h"

// Quick Sort strategy function
void strategy_quicksort(int arr[], int aux[], int n) {
    // Perform Quick Sort using the provided array
    quicksort(arr, 0, n - 1);
}
