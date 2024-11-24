//
// Created by jorda on 9/20/2024.
//

// print_radix_sort_strategy.c
#include "../include/print_radix_sort_strategy.h"
#include "../include/radix_sort.h"

// Wrapper for printing the array sorted by Radix Sort
void print_radix_sort_strategy(void* data) {
    // Extract the array and size from the data
    int* arr = ((int**)data)[0];
    int n = *((int*)(((int**)data)[1]));

    // Call the core print function
    print_radix_sort(arr, n);
}
