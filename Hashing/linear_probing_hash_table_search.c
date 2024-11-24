// linear_probing_hash_table_search.c
#include "../include/linear_probing_hash_table_search.h"
#include <stdio.h>
#include <stdlib.h>

// Initialize the hash table
void init_hash_table(LinearTableHashTable* ht, int size) {
    ht->table = (int*) malloc(sizeof(int) * size);
    ht->vacancies = (int*) malloc(sizeof(int) * size);
    ht->M = size;
    ht->N = 0;

    for (int i = 0; i < size; i++) {
        ht->table[i] = -1;       // Mark as empty
        ht->vacancies[i] = 0;    // Mark as vacant
    }
}

// Hash function (division method)
int linear_probing_hash_function(int key, int M) {
    return key % M;
}

// Search for a key using linear probing
int search_linear_probing_hash_table(LinearTableHashTable* ht, int key) {
    int i = linear_probing_hash_function(key, ht->M);
    int original_i = i;

    while (ht->table[i] != -1) {
        if (ht->table[i] == key) {
            return i;  // Key found
        }
        i = (i + 1) % ht->M;

        // If we loop back to the start, the key is not in the table
        if (i == original_i) {
            return -1;  // Key not found
        }
    }
    return -1;  // Key not found
}

// Insert a key using linear probing
int insert_linear_probing_hash_table(LinearTableHashTable* ht, int key) {
    int i = linear_probing_hash_function(key, ht->M);
    int original_i = i;

    while (ht->table[i] != -1) {
        i = (i + 1) % ht->M;

        if (i == original_i) {
            return -1;  // Table is full
        }
    }

    ht->table[i] = key;
    ht->vacancies[i] = 1;
    ht->N++;
    return i;
}

// Free the memory used by the hash table
void free_linear_probing_hash_table(LinearTableHashTable* ht) {
    free(ht->table);
    free(ht->vacancies);
}
