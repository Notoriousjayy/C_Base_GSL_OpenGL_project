//
// Created by jorda on 9/19/2024.
//

// strategy_multi_list_insertion_sort.c
#include "../include/strategy_multi_list_insertion_sort.h"
#include "../include/multi_list_insertion.h"

// Multi-list insertion sort strategy
void strategy_multi_list_insertion_sort(int keys[], int link[], int n) {
    int head[M];  // Array to hold heads of multiple lists

    // Initialize heads and links
    initialize_heads(head, link, M);

    // Call the multi-list insertion sort function
    multi_list_insertion_sort(keys, link, head, M, n);
}
