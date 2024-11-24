// list_insertion.c
#include "../include/list_insertion.h"
#include <stdio.h>

// Function to initialize the 'link' array and keys
void initialize_links(int link[], int keys[], int n) {
    for (int i = 0; i < n; i++) {
        link[i] = -1;  // Initialize links to -1 (null in this context)
    }
}

// Function for list insertion sorting
void list_insertion_sort(int keys[], int link[], int n) {
    int p, q, j, k, kp;

    for (j = 1; j < n; j++) {
        k = keys[j];
        q = -1; // Start with no previous node (null equivalent)
        p = link[0]; // Start from the first element in the linked list

        while (p != -1 && keys[p] < k) {
            q = p;
            p = link[p]; // Traverse to the next element
        }

        if (q == -1) {
            link[j] = link[0];
            link[0] = j;
        } else {
            link[j] = link[q];
            link[q] = j;
        }
    }
}

// Function to print the sorted list
void print_list_insertion(int keys[], int link[], int n) {
    int i = link[0]; // Start from the first element in the linked list
    printf("Sorted keys: ");
    while (i != -1) {
        printf("%d ", keys[i]);
        i = link[i];
    }
    printf("\n");
}
