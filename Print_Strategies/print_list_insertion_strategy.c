//
// Created by jorda on 9/19/2024.
//

// print_list_insertion_strategy.c
#include "../include/print_list_insertion_strategy.h"
#include <stdio.h>

// Wrapper for the existing print_list_insertion function
void print_list_insertion_strategy(void* data) {
    int* keys = ((int**)data)[0];
    int* link = ((int**)data)[1];
    int n = *((int*)(((int**)data)[2]));  // Number of elements

    // Call the original print_list_insertion function
    print_list_insertion(keys, link, n);
}
