//
// Created by jorda on 9/20/2024.
//

// strategy_hash_table_search.c
#include <stdlib.h>
#include "../include/strategy_chained_hash_table_search.h"
#include "../include/chained_hash_table_search.h"

// Hash table search strategy
void strategy_hash_table_search(void* data, const char* key_str, int* result) {
    // Convert the string key to an integer
    int key = atoi(key_str);

    // Extract the hash table
    HashTable* ht = (HashTable*)data;

    // Perform the search
    int index = search_hash_table(ht, key);
    *result = (index != -1) ? index : -1; // Return index if found, otherwise -1
}
