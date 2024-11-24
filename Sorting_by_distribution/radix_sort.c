//
// Created by jorda on 9/20/2024.
//

// radix_sort.c
#include "../include/radix_sort.h"
#include <stdio.h>

// Get the kth digit of the number, assuming each digit is a byte
int get_digit(int num, int k) {
    return (num >> (8 * k)) & 0xFF;
}

// Radix List Sort implementation
void radix_sort(int arr[], int n, int p) {
    int output[n];  // Output array to store sorted numbers
    int count[MAX_DIGIT];  // Count array to store the count of each digit

    // Perform sorting on each digit (from least significant to most significant)
    for (int k = 0; k < p; k++) {
        // Initialize count array
        for (int i = 0; i < MAX_DIGIT; i++) {
            count[i] = 0;
        }

        // Count occurrences of each digit
        for (int i = 0; i < n; i++) {
            int digit = get_digit(arr[i], k);
            count[digit]++;
        }

        // Change count[i] so that count[i] contains the position of the digit in output[]
        for (int i = 1; i < MAX_DIGIT; i++) {
            count[i] += count[i - 1];
        }

        // Build the output array
        for (int i = n - 1; i >= 0; i--) {
            int digit = get_digit(arr[i], k);
            output[count[digit] - 1] = arr[i];
            count[digit]--;
        }

        // Copy the output array to arr[], so that arr[] contains the sorted numbers
        for (int i = 0; i < n; i++) {
            arr[i] = output[i];
        }
    }
}

// Function to print an array
void print_radix_sort(int arr[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}
