//
// Created by jorda on 9/20/2024.
//

// strategy_linear_probing_hash_table_search.c
#include <stdlib.h>
#include "../include/strategy_linear_probing_hash_table_search.h"
#include "../include/linear_probing_hash_table_search.h"

// Linear probing hash table search strategy
void strategy_linear_probing_hash_table_search(void* data, const char* key_str, int* result) {
    // Convert the string key to an integer
    int key = atoi(key_str);

    // Extract the hash table
    LinearTableHashTable* ht = (LinearTableHashTable*)data;

    // Perform the search
    int index = search_linear_probing_hash_table(ht, key);
    *result = (index != -1) ? index : -1;  // Return the index if found, otherwise -1
}