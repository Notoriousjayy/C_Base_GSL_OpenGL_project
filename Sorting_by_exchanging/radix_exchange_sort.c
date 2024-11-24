//
// Created by jorda on 9/20/2024.
//

// radix_exchange_sort.c
#include "../include/radix_exchange_sort.h"
#include <stdio.h>

// Function to perform the radix exchange sort recursively
void radix_exchange_sort(int arr[], int left, int right, int bit) {
    if (left >= right || bit < 0) {
        return;  // Base case: no more bits or no more range to sort
    }

    int i = left, j = right;
    int mask = 1 << bit;  // Create a mask for the current bit

    while (i <= j) {
        // Move 'i' to the right while arr[i] has a 0 in the current bit
        while (i <= j && (arr[i] & mask) == 0) {
            i++;
        }
        // Move 'j' to the left while arr[j] has a 1 in the current bit
        while (i <= j && (arr[j] & mask) != 0) {
            j--;
        }
        // Swap arr[i] and arr[j] if needed
        if (i < j) {
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
            i++;
            j--;
        }
    }

    // Recursively sort the left and right partitions
    radix_exchange_sort(arr, left, j, bit - 1);   // Sort elements with 0 in the current bit
    radix_exchange_sort(arr, i, right, bit - 1);  // Sort elements with 1 in the current bit
}

// Function to print an array
void print_radix_exchange_sort(int arr[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}
