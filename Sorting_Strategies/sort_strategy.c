//
// Created by jorda on 9/19/2024.
//

#include "../include/sort_strategy.h"

// Function to set the current sorting strategy
void set_sort_strategy(SortContext* context, SortStrategy sort_strategy) {
    context->sort_strategy = sort_strategy;
}

// Function to execute the current sorting strategy
void execute_sort_strategy(SortContext* context, int keys[], int aux[], int n) {
    if (context->sort_strategy) {
        context->sort_strategy(keys, aux, n);
    }
}
