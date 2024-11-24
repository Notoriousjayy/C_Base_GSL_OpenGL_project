//
// Created by jorda on 9/20/2024.
//

// selection_sort.c
#include "../include/selection_sort.h"
#include <stdio.h>

// Function to perform selection sort
void selection_sort(int arr[], int n) {
    for (int j = n - 1; j > 0; j--) {
        int max_index = j;
        // Find the maximum element in arr[0...j]
        for (int i = 0; i < j; i++) {
            if (arr[i] > arr[max_index]) {
                max_index = i;
            }
        }
        // Swap the found maximum element with the element at index j
        int temp = arr[j];
        arr[j] = arr[max_index];
        arr[max_index] = temp;
    }
}

// Function to print an array
void print_selection_sort(int arr[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}
