//
// Created by jorda on 9/19/2024.
//

#include <stdio.h>
#include "../include/comparison_counting_sort.h"
//comparison_counting_sort.c

// Function to perform the comparison counting sort
void comparison_counting_sort(int keys[], int count[], int n) {
    // Initialize the counts to 0
    for (int i = 0; i < n; i++) {
        count[i] = 0;
    }

    // Compare each key with every other key
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (keys[i] > keys[j]) {
                count[i]++; // Increment the count for keys[i] if it's larger
            }
        }
    }
}

// Function to print the sorted list using the count array
void print_comparison_counting_sort(int keys[], int count[], int n) {
    printf("Sorted keys (Comparison Counting Sort): ");

    // Print the keys in sorted order based on the count array
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (count[j] == i) {
                printf("%d ", keys[j]);
            }
        }
    }
    printf("\n");
}