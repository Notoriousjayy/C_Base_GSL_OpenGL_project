//
// Created by jorda on 9/20/2024.
//

// strategy_binary_search.c
#include "../include/strategy_binary_search.h"
#include "../include/binary_search.h"

// Binary search strategy
void strategy_binary_search(int arr[], int aux[], int n, int key, int* result) {
    *result = binary_search(arr, n, key);
}
