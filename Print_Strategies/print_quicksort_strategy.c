//
// Created by jorda on 9/20/2024.
//

// print_quicksort_strategy.c
#include "../include/print_quicksort_strategy.h"
#include "../include/quicksort.h"

// Wrapper for printing the array sorted by Quick Sort
void print_quicksort_strategy(void* data) {
    // Extract the array and size from the data
    int* arr = ((int**)data)[0];
    int n = *((int*)(((int**)data)[1]));

    // Call the core print function
    print_quicksort(arr, n);
}
