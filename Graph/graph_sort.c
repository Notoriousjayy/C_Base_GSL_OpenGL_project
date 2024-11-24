//
// Created by jorda on 7/23/2024.
//

//
// Created by jorda on 7/23/2024.
//

#include <stdio.h>
#include "../include/graph_flip.h"

typedef struct node_struct {
    long key;
    struct node_struct* link;
} node;

node* graph_sorted[256];
static node* alt_sorted[256];

void graph_linksort(node* l) {
    long k;
    node** pp;
    node* p;
    node* q;

    // Initialize alt_sorted array to NULL
    for (pp = alt_sorted + 255; pp >= alt_sorted; pp--) {
        *pp = NULL;
    }

    // First pass of distribution using graph_next_rand
    for (p = l; p; p = q) {
        k = graph_next_rand() >> 23;
        q = p->link;
        p->link = alt_sorted[k];
        alt_sorted[k] = p;
    }

    // Initialize graph_sorted array to NULL
    for (pp = graph_sorted + 255; pp >= graph_sorted; pp--) {
        *pp = NULL;
    }

    // Second pass of distribution using graph_next_rand
    for (pp = alt_sorted + 255; pp >= alt_sorted; pp--) {
        for (p = *pp; p; p = q) {
            k = graph_next_rand() >> 23;
            q = p->link;
            p->link = graph_sorted[k];
            graph_sorted[k] = p;
        }
    }

    // Initialize alt_sorted array to NULL
    for (pp = alt_sorted + 255; pp >= alt_sorted; pp--) {
        *pp = NULL;
    }

    // Third pass of distribution using the least significant byte of the key
    for (pp = graph_sorted + 255; pp >= graph_sorted; pp--) {
        for (p = *pp; p; p = q) {
            k = p->key & 0xff;
            q = p->link;
            p->link = alt_sorted[k];
            alt_sorted[k] = p;
        }
    }

    // Initialize graph_sorted array to NULL
    for (pp = graph_sorted + 255; pp >= graph_sorted; pp--) {
        *pp = NULL;
    }

    // Fourth pass of distribution using the next byte of the key
    for (pp = alt_sorted; pp < alt_sorted + 256; pp++) {
        for (p = *pp; p; p = q) {
            k = (p->key >> 8) & 0xff;
            q = p->link;
            p->link = graph_sorted[k];
            graph_sorted[k] = p;
        }
    }

    // Initialize alt_sorted array to NULL
    for (pp = alt_sorted + 255; pp >= alt_sorted; pp--) {
        *pp = NULL;
    }

    // Fifth pass of distribution using the next byte of the key
    for (pp = graph_sorted + 255; pp >= graph_sorted; pp--) {
        for (p = *pp; p; p = q) {
            k = (p->key >> 16) & 0xff;
            q = p->link;
            p->link = alt_sorted[k];
            alt_sorted[k] = p;
        }
    }

    // Initialize graph_sorted array to NULL
    for (pp = graph_sorted + 255; pp >= graph_sorted; pp--) {
        *pp = NULL;
    }

    // Sixth pass of distribution using the most significant byte of the key
    for (pp = alt_sorted; pp < alt_sorted + 256; pp++) {
        for (p = *pp; p; p = q) {
            k = (p->key >> 24) & 0xff;
            q = p->link;
            p->link = graph_sorted[k];
            graph_sorted[k] = p;
        }
    }
}
// Function to print the linked list

void print_list_graph(node* head) {
    node* current = head;
    while (current != NULL) {
        printf("%ld -> ", current->key);
        current = current->link;
    }
    printf("NULL\n");
}



