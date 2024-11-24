//
// Created by jorda on 9/19/2024.
//

// print_comparison_counting_sort_strategy.c
#include "../include/print_comparison_counting_sort_strategy.h"
#include <stdio.h>

// Wrapper for the comparison counting sort function
void print_comparison_counting_sort_strategy(void* data) {
    // Extract keys, count, and n from the passed data
    int* keys = ((int**)data)[0];
    int* count = ((int**)data)[1];
    int n = *((int*)(((int**)data)[2]));  // Number of elements

    // Perform the comparison counting sort
    comparison_counting_sort(keys, count, n);
    print_comparison_counting_sort(keys, count, n);
}
