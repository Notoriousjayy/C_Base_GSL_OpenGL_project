//
// Created by jorda on 9/20/2024.
//

// binary_search.c
#include "../include/binary_search.h"
#include <stdio.h>

// Binary Search function that searches for the key in the sorted array
int binary_search(int arr[], int n, int key) {
    int left = 0;
    int right = n - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;

        // Check if the key is at the midpoint
        if (arr[mid] == key) {
            return mid;
        }

        // If key is greater, ignore the left half
        if (arr[mid] < key) {
            left = mid + 1;
        }
            // If key is smaller, ignore the right half
        else {
            right = mid - 1;
        }
    }

    // Return -1 if key is not found
    return -1;
}

// Function to print the result of the search
void print_binary_search_result(int result, int key) {
    if (result == -1) {
        printf("Key %d not found in the array.\n", key);
    } else {
        printf("Key %d found at index %d.\n", key, result);
    }
}
