//
// Created by jorda on 9/20/2024.
//

// fibonacci_search.c
#include "../include/fibonacci_search.h"
#include <stdio.h>

// Function to generate Fibonacci numbers up to n
int fibonacci(int n) {
    if (n == 0) return 0;
    if (n == 1) return 1;
    return fibonacci(n - 1) + fibonacci(n - 2);
}

// Fibonacci Search function
int fibonacci_search(int arr[], int n, int key) {
    int fibM2 = 0;  // (m-2)'th Fibonacci number
    int fibM1 = 1;  // (m-1)'th Fibonacci number
    int fibM = fibM2 + fibM1;  // m'th Fibonacci number

    int offset = -1;

    // Initialize Fibonacci numbers to find the smallest Fibonacci number greater than or equal to n
    while (fibM < n) {
        fibM2 = fibM1;
        fibM1 = fibM;
        fibM = fibM2 + fibM1;
    }

    // While there are elements to be inspected
    while (fibM > 1) {
        int i = (offset + fibM2 < n - 1) ? (offset + fibM2) : (n - 1);

        // If key is greater than the value at index i, cut the subarray from offset to i
        if (arr[i] < key) {
            fibM = fibM1;
            fibM1 = fibM2;
            fibM2 = fibM - fibM1;
            offset = i;
        }
            // If key is less than the value at index i, cut the subarray after i
        else if (arr[i] > key) {
            fibM = fibM2;
            fibM1 = fibM1 - fibM2;
            fibM2 = fibM - fibM1;
        }
            // Key found
        else {
            return i;
        }
    }

    // Comparing the last element with key
    if (fibM1 && arr[offset + 1] == key) {
        return offset + 1;
    }

    return -1;  // Return -1 if key is not found
}

// Function to print the result of the search
void print_fibonacci_search_result(int result, int key) {
    if (result == -1) {
        printf("Key %d not found in the array.\n", key);
    } else {
        printf("Key %d found at index %d.\n", key, result);
    }
}
