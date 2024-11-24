//
// Created by jorda on 9/20/2024.
//

// strategy_uniform_binary_search.c
#include "../include/strategy_uniform_binary_search.h"
#include "../include/uniform_binary_search.h"

// Uniform binary search strategy
void strategy_uniform_binary_search(int arr[], int aux[], int n, int key, int* result) {
    // Assume aux contains the deltas array
    *result = uniform_binary_search(arr, n, key, aux);
}
