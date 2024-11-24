//
// Created by jorda on 9/20/2024.
//

// strategy_sequential_search.c
#include "../include/strategy_sequential_search.h"
#include "../include/sequential_search.h"

// Sequential search strategy
void strategy_sequential_search(int arr[], int aux[], int n, int key) {
    // Perform sequential search
    int result = sequential_search(arr, n, key);

    // Print the result using the print function
    print_sequential_search_result(result, key);
}

