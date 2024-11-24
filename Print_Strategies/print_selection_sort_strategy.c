//
// Created by jorda on 9/20/2024.
//

// print_selection_sort_strategy.c
#include "../include/print_selection_sort_strategy.h"
#include "../include/selection_sort.h"

// Wrapper for printing the array sorted by Selection Sort
void print_selection_sort_strategy(void* data) {
    // Extract the array and size from the data
    int* arr = ((int**)data)[0];
    int n = *((int*)(((int**)data)[1]));

    // Call the core print function
    print_selection_sort(arr, n);
}
