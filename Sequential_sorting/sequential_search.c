//
// Created by jorda on 9/20/2024.
//

// sequential_search.c
#include "../include/sequential_search.h"
#include <stdio.h>

// Sequential search function that searches for the key in the array
int sequential_search(int arr[], int n, int key) {
    for (int i = 0; i < n; i++) {
        if (arr[i] == key) {
            return i;  // Return the index if key is found
        }
    }
    return -1;  // Return -1 if key is not found
}

// Function to print the result of the search
void print_sequential_search_result(int result, int key) {
    if (result == -1) {
        printf("Key %d not found in the array.\n", key);
    } else {
        printf("Key %d found at index %d.\n", key, result);
    }
}
