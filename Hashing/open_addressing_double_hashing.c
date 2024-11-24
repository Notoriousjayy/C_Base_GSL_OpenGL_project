// open_addressing_double_hashing.c
#include "../include/open_addressing_double_hashing.h"
#include <stdio.h>
#include <stdlib.h>

// Create and initialize the hash table
HashTable *create_table(int size) {
    HashTable *new_table = (HashTable *)malloc(sizeof(HashTable));
    new_table->size = size;
    new_table->vacancies = size;
    new_table->table = (HashEntry *)malloc(size * sizeof(HashEntry));

    for (int i = 0; i < size; i++) {
        new_table->table[i].key = EMPTY;
        new_table->table[i].value = 0;
    }

    return new_table;
}

// Destroy the hash table and free memory
void destroy_table(HashTable *table) {
    free(table->table);
    free(table);
}

// First hash function (hash1)
int hash1(int key, int size) {
    return key % size;
}

// Second hash function (hash2) for double hashing
int hash2(int key, int size) {
    return 1 + (key % (size - 1));
}

// Insert a key-value pair into the hash table
int insert_open_addressing_double_hashing(HashTable *table, int key, int value) {
    if (table->vacancies == 0) {
        printf("Hash table is full.\n");
        return 0;
    }

    int index = hash1(key, table->size);
    int step = hash2(key, table->size);

    while (table->table[index].key != EMPTY && table->table[index].key != DELETED) {
        index = (index + step) % table->size;
    }

    table->table[index].key = key;
    table->table[index].value = value;
    table->vacancies--;
    return 1;
}

// Search for a key in the hash table
int search_open_addressing_double_hashing(HashTable *table, int key) {
    int index = hash1(key, table->size);
    int step = hash2(key, table->size);

    while (table->table[index].key != EMPTY) {
        if (table->table[index].key == key) {
            return table->table[index].value;
        }
        index = (index + step) % table->size;
    }

    return EMPTY;  // Key not found
}

// Delete an entry in the hash table by key
int delete_entry(HashTable *table, int key) {
    int index = hash1(key, table->size);
    int step = hash2(key, table->size);

    while (table->table[index].key != EMPTY) {
        if (table->table[index].key == key) {
            table->table[index].key = DELETED;
            table->vacancies++;
            return 1;  // Successfully deleted
        }
        index = (index + step) % table->size;
    }

    return 0;  // Key not found
}
