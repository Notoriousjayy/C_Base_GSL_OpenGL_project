//
// Created by jorda on 9/20/2024.
//

// quicker_sequential_search.c
#include "../include/quicker_sequential_search.h"
#include <stdio.h>

// Quicker Sequential Search function that searches for the key in the array
int quicker_sequential_search(int arr[], int n, int key) {
    for (int i = 0; i < n; i += 2) {
        // Compare the current and the next element in each iteration
        if (arr[i] == key) {
            return i;  // Return the index if key is found
        }
        if (i + 1 < n && arr[i + 1] == key) {
            return i + 1;  // Return the index if the key is found in the next element
        }
    }
    return -1;  // Return -1 if key is not found
}

// Function to print the result of the search
void print_quicker_sequential_search_result(int result, int key) {
    if (result == -1) {
        printf("Key %d not found in the array.\n", key);
    } else {
        printf("Key %d found at index %d.\n", key, result);
    }
}
