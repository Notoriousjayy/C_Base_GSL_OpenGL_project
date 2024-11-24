//
// Created by jorda on 9/19/2024.
//

// print_multi_list_strategy.c
#include "../include/print_multi_list_strategy.h"
#include <stdio.h>
#include "../include/multi_list_insertion.h"

// Wrapper for printing sorted lists to match the PrintListStrategy signature
void print_sorted_multi_lists(void* data) {
    // Extract the keys, link, and head arrays from data
    int* keys = ((int**)data)[0];
    int* link = ((int**)data)[1];
    int* head = ((int**)data)[2];

    // Call the existing function to print sorted lists
    print_sorted_lists(keys, link, head, M);
}
