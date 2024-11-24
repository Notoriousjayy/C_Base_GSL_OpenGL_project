//
// Created by jorda on 9/20/2024.
//

#include "../include/chained_hash_table_search.h"
#include <stdlib.h>
#include <stdio.h>

// Create a new hash table
HashTable* create_hash_table(int size) {
    HashTable* ht = (HashTable*)malloc(sizeof(HashTable));
    ht->table = (Node*)calloc(size, sizeof(Node));
    ht->M = size;
    return ht;
}

// Hash function (modulo division)
int hash_function(int key, int M) {
    return key % M;
}

// Insert a key into the hash table using chaining
void insert_hash_table(HashTable* ht, int key) {
    int i = hash_function(key, ht->M);
    int original_i = i; // Save original index
    while (ht->table[i].key != 0) { // Find an empty spot
        i = (i + 1) % ht->M;
        if (i == original_i) {
            printf("Hash table is full!\n");
            return;
        }
    }
    ht->table[i].key = key;
    ht->table[i].link = -1; // No linked node for simplicity
    printf("Inserted key %d at index %d\n", key, i);
}

// Search for a key in the hash table
int search_hash_table(HashTable* ht, int key) {
    int i = hash_function(key, ht->M);
    while (ht->table[i].key != 0) {
        if (ht->table[i].key == key) {
            return i; // Key found
        }
        i = (i + 1) % ht->M; // Linear probing
    }
    return -1; // Key not found
}

// Free the hash table memory
void free_hash_table(HashTable* ht) {
    free(ht->table);
    free(ht);
}
