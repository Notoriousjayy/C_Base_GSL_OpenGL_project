//
// Created by jorda on 9/20/2024.
//

// print_shellsort_strategy.c
#include "../include/print_shellsort_strategy.h"
#include <stdio.h>
#include "../include/shellsort.h"

// Wrapper for printing the array sorted by Shellsort
void print_shellsort_strategy(void* data) {
    // Extract the array and size from the data
    int* arr = ((int**)data)[0];
    int n = *((int*)(((int**)data)[1]));

    // Call the core print function
    print_array(arr, n);
}
