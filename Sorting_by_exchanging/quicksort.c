//
// Created by jorda on 9/20/2024.
//

// quicksort.c
#include "../include/quicksort.h"
#include <stdio.h>

// Function to perform the partitioning for Quicksort
int partition(int arr[], int left, int right) {
    int pivot = arr[right];
    int i = left - 1;
    int temp;

    for (int j = left; j < right; j++) {
        if (arr[j] <= pivot) {
            i++;
            // Swap arr[i] and arr[j]
            temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }

    // Swap arr[i + 1] and arr[right] (or pivot)
    temp = arr[i + 1];
    arr[i + 1] = arr[right];
    arr[right] = temp;

    return i + 1;
}

// Recursive function to perform Quicksort
void quicksort(int arr[], int left, int right) {
    if (left < right) {
        // Partition the array
        int pivot_index = partition(arr, left, right);

        // Recursively sort the left and right partitions
        quicksort(arr, left, pivot_index - 1);
        quicksort(arr, pivot_index + 1, right);
    }
}

// Function to print the array
void print_quicksort(int arr[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}
