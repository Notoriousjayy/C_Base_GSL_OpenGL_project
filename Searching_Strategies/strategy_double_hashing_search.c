//
// Created by jorda on 9/20/2024.
//

// strategy_double_hashing_search.c
#include "../include/strategy_double_hashing_search.h"
#include "../include/open_addressing_double_hashing.h"
#include <stdlib.h>
#include <stdio.h>

// Search strategy for double hashing
void strategy_double_hashing_search(void* data, const char* key_str, int* result) {
    // Convert key from string to int
    int key = atoi(key_str);

    // Cast the data to HashTable*
    HashTable* table = (HashTable*)data;

    // Perform the search using double hashing
    int value = search_open_addressing_double_hashing(table, key);

    *result = (value != EMPTY) ? value : -1;
}
