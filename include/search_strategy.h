//
// Created by jorda on 9/20/2024.
//

// search_strategy.h
#ifndef SEARCH_STRATEGY_H
#define SEARCH_STRATEGY_H

// Define a function pointer type for searching strategies
typedef void (*SearchStrategy)(void* data, const char* key, int* result);

// Context to hold the current searching strategy
typedef struct {
    SearchStrategy search_strategy;
} SearchContext;

// Function to set the current searching strategy
void set_search_strategy(SearchContext* context, SearchStrategy search_strategy);

// Function to execute the current searching strategy
void execute_search_strategy(SearchContext* context, void* data, const char* key, int* result);

#endif // SEARCH_STRATEGY_H


