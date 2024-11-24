//
// Created by jorda on 9/19/2024.
//

#include "../include/strategy_list_insertion.h"
#include "../include/list_insertion.h"

// Wrapper for the insertion sort strategy
void strategy_list_insertion_sort(int keys[], int link[], int n) {
    initialize_links(link, keys, n);
    list_insertion_sort(keys, link, n);
}
