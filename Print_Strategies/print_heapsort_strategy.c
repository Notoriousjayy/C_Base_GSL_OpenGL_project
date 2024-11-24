//
// Created by jorda on 9/20/2024.
//

// print_heapsort_strategy.c
#include "../include/print_heapsort_strategy.h"
#include "../include/heapsort.h"

// Wrapper for printing the array sorted by Heap Sort
void print_heapsort_strategy(void* data) {
    // Extract the array and size from the data
    int* arr = ((int**)data)[0];
    int n = *((int*)(((int**)data)[1]));

    // Call the core print function
    print_heapsort(arr, n);
}
