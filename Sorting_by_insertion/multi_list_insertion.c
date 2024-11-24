//
// Created by jorda on 9/19/2024.
//

// multi_list_insertion.c
#include "../include/multi_list_insertion.h"
#include <stdio.h>

// Function to initialize heads and links
void initialize_heads(int head[], int link[], int m) {
    for (int i = 0; i < m; i++) {
        head[i] = -1;  // Initialize heads to -1 (null in this context)
    }
    for (int i = 0; i < N; i++) {
        link[i] = -1;  // Initialize links to -1 (null in this context)
    }
}

// Function for multiple list insertion sorting
void multi_list_insertion_sort(int keys[], int link[], int head[], int m, int n) {
    int p, q, j, k;
    int index;

    for (j = 0; j < n; j++) {
        k = keys[j];
        index = (k * m) >> (E - 3); // M * K_j / 2^E
        q = -1; // Start with no previous node
        p = head[index]; // Start from the first element in the corresponding list

        while (p != -1 && keys[p] < k) {
            q = p;
            p = link[p]; // Traverse to the next element
        }

        if (q == -1) {
            link[j] = head[index];
            head[index] = j;
        } else {
            link[j] = link[q];
            link[q] = j;
        }
    }
}

// Function to print sorted keys based on the lists
void print_sorted_lists(int keys[], int link[], int head[], int m) {
    for (int i = 0; i < m; i++) {
        int p = head[i];
        printf("List %d: ", i);
        while (p != -1) {
            printf("%d ", keys[p]);
            p = link[p];
        }
        printf("\n");
    }
}
