//
// Created by jorda on 9/19/2024.
//

// multi_list_insertion.h
#ifndef MULTI_LIST_INSERTION_H
#define MULTI_LIST_INSERTION_H

#define N 16   // Example size for N; can be changed
#define M 4    // Number of lists for the multi-list insertion
#define E 64   // Assuming keys are less than 2^e

// Function declarations
void initialize_heads(int head[], int link[], int m);
void multi_list_insertion_sort(int keys[], int link[], int head[], int m, int n);
void print_sorted_lists(int keys[], int link[], int head[], int m);

#endif
