//
// Created by jorda on 9/19/2024.
//

// list_print_strategy.c
#include <stdio.h>
#include "../include/print_strategy.h"

// Function to set the current strategy
void set_print_strategy(PrintContext* context, PrintListStrategy strategy) {
    context->strategy = strategy;
}

// Function to execute the current strategy
void execute_print_strategy(PrintContext* context, void* data) {
    if (context->strategy != NULL) {
        context->strategy(data);
    } else {
        printf("No print strategy set.\n");
    }
}
