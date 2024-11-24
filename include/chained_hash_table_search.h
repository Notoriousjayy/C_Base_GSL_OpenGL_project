//
// Created by jorda on 9/20/2024.
//

#ifndef HASH_TABLE_H
#define HASH_TABLE_H

// Define constants and structures for the hash table
#define TABLE_SIZE 256 // Example size of the table

typedef struct {
    int key;
    int link; // Link field for chained hash table
} Node;

typedef struct {
    Node* table; // Pointer to the haash table array
    int M;       // The size of the hash table
} HashTable;

// Function prototypes
HashTable* create_hash_table(int size);
int hash_function(int key, int M);
void insert_hash_table(HashTable* ht, int key);
int search_hash_table(HashTable* ht, int key);
void free_hash_table(HashTable* ht);

#endif
