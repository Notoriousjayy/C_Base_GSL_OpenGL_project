//
// Created by jorda on 9/20/2024.
//

#include "../include/search_strategy.h"
#include <stddef.h>  // Add this line to define NULL

// Function to set the current searching strategy
void set_search_strategy(SearchContext* context, SearchStrategy search_strategy) {
    context->search_strategy = search_strategy;
}

// Function to execute the current searching strategy
void execute_search_strategy(SearchContext* context, void* data, const char* key, int* result) {
    if (context->search_strategy != NULL) {
        context->search_strategy(data, key, result);
    }
}