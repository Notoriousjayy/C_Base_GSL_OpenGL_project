//
// Created by jorda on 9/20/2024.
//

// strategy_quicker_sequential_search.c
#include "../include/strategy_quicker_sequential_search.h"
#include "../include/quicker_sequential_search.h"


// Quicker sequential search strategy
void strategy_quicker_sequential_search(int arr[], int aux[], int n, int key, int* result) {
    *result = quicker_sequential_search(arr, n, key);
}
