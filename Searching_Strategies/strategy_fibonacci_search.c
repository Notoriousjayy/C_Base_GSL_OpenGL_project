//
// Created by jorda on 9/20/2024.
//

// strategy_fibonacci_search.c
#include "../include/strategy_fibonacci_search.h"
#include "../include/fibonacci_search.h"

// Fibonacci search strategy
void strategy_fibonacci_search(int arr[], int aux[], int n, int key, int* result) {
    *result = fibonacci_search(arr, n, key);
}
