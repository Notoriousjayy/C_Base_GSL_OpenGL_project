//
// Created by jorda on 9/20/2024.
//

// bubble_sort.c
#include "../include/bubble_sort.h"
#include <stdio.h>

// Bubble Sort implementation
void bubble_sort(int arr[], int n) {
    int bound, j, temp, exchanges;

    bound = n - 1;  // Initialize BOUND
    do {
        exchanges = 0;
        for (j = 0; j < bound; j++) {
            if (arr[j] > arr[j + 1]) {
                // Swap arr[j] and arr[j + 1]
                temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;

                exchanges = 1;  // Record that a swap happened
            }
        }
        bound--;  // Decrease the boundary after each pass
    } while (exchanges);  // Continue until no swaps are made
}

// Function to print the array
void print_bubble_sort(int arr[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}
