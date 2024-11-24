//
// Created by jorda on 9/19/2024.
//

// list_print_strategy.h
#ifndef LIST_PRINT_STRATEGY_H
#define LIST_PRINT_STRATEGY_H

// Define a function pointer type for print strategies
typedef void (*PrintListStrategy)(void* data);

// Context to hold the current strategy
typedef struct {
    PrintListStrategy strategy;
} PrintContext;

// Function to set the current strategy
void set_print_strategy(PrintContext* context, PrintListStrategy strategy);

// Function to execute the current strategy
void execute_print_strategy(PrintContext* context, void* data);

#endif // LIST_PRINT_STRATEGY_H
