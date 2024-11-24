//
// Created by jorda on 9/20/2024.
//

// strategy_merge_sort.c
#include "../include/strategy_merge_sort.h"
#include "../include/merge_sort.h"

// Merge Sort strategy function
void strategy_merge_sort(int arr[], int aux[], int n) {
    merge_sort(arr, 0, n - 1);
}
