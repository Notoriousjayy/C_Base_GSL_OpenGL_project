//
// Created by jorda on 9/20/2024.
//

// strategy_radix_exchange_sort.c
#include "../include/strategy_radix_exchange_sort.h"
#include "../include/radix_exchange_sort.h"

// Radix Exchange Sort strategy function
void strategy_radix_exchange_sort(int arr[], int aux[], int n) {
    int max_bit = 31;  // Assuming 32-bit integers, we start sorting with the highest bit
    radix_exchange_sort(arr, 0, n - 1, max_bit);
}
