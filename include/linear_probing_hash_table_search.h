// linear_probing_hash_table_search.h
#ifndef LINEAR_PROBING_HASH_TABLE_SEARCH_H
#define LINEAR_PROBING_HASH_TABLE_SEARCH_H

// Define a hash table structure
typedef struct {
    int* table;       // Pointer to the array of keys
    int* vacancies;   // Pointer to the array of vacancies
    int M;            // Table size (assumed to be prime)
    int N;            // Number of keys inserted
} LinearTableHashTable;

// Initialize the hash table
void init_hash_table(LinearTableHashTable* ht, int size);

// Hash function (division method)
int linear_probing_hash_function(int key, int M);

// Search for a key in the hash table using linear probing
int search_linear_probing_hash_table(LinearTableHashTable* ht, int key);

// Insert a key into the hash table using linear probing
int insert_linear_probing_hash_table(LinearTableHashTable* ht, int key);

// Free the hash table memory
void free_linear_probing_hash_table(LinearTableHashTable* ht);

#endif
