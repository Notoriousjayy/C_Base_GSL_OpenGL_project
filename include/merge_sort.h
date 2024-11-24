//
// Created by jorda on 9/20/2024.
//
// merge_sort.h
#ifndef MERGE_SORT_H
#define MERGE_SORT_H

#define N 16  // Example size for N (array size)

void merge_sort(int arr[], int left, int right);
void merge(int arr[], int left, int mid, int right);
void print_merge_sort(int arr[], int n);

#endif
