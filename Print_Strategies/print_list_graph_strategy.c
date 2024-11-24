//
// Created by jorda on 9/19/2024.
//

// print_list_graph_strategy.c
#include "../include/print_list_graph_strategy.h"
#include "../include/graph_sort.h"
#include <stdio.h>

// Wrapper for the existing print_list_graph function
void print_list_graph_strategy(void* data) {
    node* head = (node*)data;  // Cast the data to the correct type

    // Call the original print_list_graph function
    print_list_graph(head);
}
