//
// Created by jorda on 7/23/2024.
//

#ifndef GB_SORT_H
#define GB_SORT_H

typedef struct node_struct {
    long key;
    struct node_struct* link;
} node;

// Function declarations
void graph_linksort(node *);
extern char* graph_sorted[];
void print_list_graph(node* head);

#endif // GB_SORT_H
