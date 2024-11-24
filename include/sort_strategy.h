//
// Created by jorda on 9/19/2024.
//

#ifndef SORT_STRATEGY_H
#define SORT_STRATEGY_H

// Define a function pointer type for sorting strategies
typedef void (*SortStrategy)(int keys[], int aux[], int n);

// Context to hold the current sorting strategy
typedef struct {
    SortStrategy sort_strategy;
} SortContext;

// Function to set the current sorting strategy
void set_sort_strategy(SortContext* context, SortStrategy sort_strategy);

// Function to execute the current sorting strategy
void execute_sort_strategy(SortContext* context, int keys[], int aux[], int n);

#endif // SORT_STRATEGY_H
