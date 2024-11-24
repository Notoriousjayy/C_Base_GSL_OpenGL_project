//
// Created by jorda on 9/20/2024.
//

#ifndef HASH_TABLE_H
#define HASH_TABLE_H

#define TABLE_SIZE 1000  // Adjust based on needs
#define EMPTY -1
#define DELETED -2

typedef struct {
    int key;
    int value;
} HashEntry;

typedef struct {
    HashEntry *table;
    int size;
    int vacancies;
} HashTable;

// Hash table operations
HashTable *create_table(int size);
void destroy_table(HashTable *table);
int hash1(int key, int size);
int hash2(int key, int size);
int insert_open_addressing_double_hashing(HashTable *table, int key, int value);
int search_open_addressing_double_hashing(HashTable *table, int key);
int delete_entry(HashTable *table, int key);

#endif