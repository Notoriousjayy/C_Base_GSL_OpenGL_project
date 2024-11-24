//
// Created by jorda on 9/20/2024.
//

// uniform_binary_search.c
#include "../include/uniform_binary_search.h"
#include <stdio.h>

// Uniform Binary Search function
int uniform_binary_search(int arr[], int n, int key, int deltas[]) {
    int j = 0;  // Starting index in the deltas array
    int i = deltas[j];  // Start at the midpoint defined by the delta

    while (deltas[j] != 0) {
        // Check if the current index is within bounds
        if (i < 0 || i >= n) {
            return -1;  // Return -1 if out of bounds
        }

        // Compare key with the element at index i
        if (arr[i] == key) {
            return i;  // Return index if the key is found
        }

        // If key is greater, move to the right
        if (arr[i] < key) {
            i += deltas[++j];  // Move right and go to the next delta
        }
            // If key is smaller, move to the left
        else {
            i -= deltas[++j];  // Move left and go to the next delta
        }
    }

    return -1;  // Return -1 if the key is not found
}

// Function to print the result of the search
void print_uniform_binary_search_result(int result, int key) {
    if (result == -1) {
        printf("Key %d not found in the array.\n", key);
    } else {
        printf("Key %d found at index %d.\n", key, result);
    }
}
