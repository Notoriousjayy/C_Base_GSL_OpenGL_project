//
// Created by jorda on 9/20/2024.
//

// shellsort.c
#include "../include/shellsort.h"
#include <stdio.h>

// Function to perform Shellsort
void shellsort(int arr[], int n, int inc[], int t) {
    int h, i, j, temp;

    // Loop over each increment in reverse order
    for (int s = t - 1; s >= 0; s--) {
        h = inc[s];  // Get the increment value from the table

        // Perform insertion sort with gap h
        for (i = h; i < n; i++) {
            temp = arr[i];
            j = i;

            // Compare and swap elements separated by gap h
            while (j >= h && arr[j - h] > temp) {
                arr[j] = arr[j - h];
                j -= h;
            }

            // Place the element in its correct position
            arr[j] = temp;
        }
    }
}

// Function to print an array
void print_array(int arr[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}
