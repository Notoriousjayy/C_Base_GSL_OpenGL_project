#include "../include/strategy_btree_search.h"
#include <stdio.h>

// B-Tree search strategy implementation using the SearchContext
void strategy_btree_search(void* data, void* key, int* result) {
    BTreeNode* root = (BTreeNode*)data;   // Cast data back to BTreeNode*
    int search_key = *(int*)key;          // Cast key back to int
    printf("Searching for key %d in the B-Tree...\n", search_key);

    // Perform the search and check if the node is found
    BTreeNode* found_node = btree_search(root, search_key);

    if (found_node != NULL) {
        *result = 1;  // Key found
        printf("Key %d found in the B-Tree.\n", search_key);
    } else {
        *result = -1;  // Key not found
        printf("Key %d not found in the B-Tree.\n", search_key);
    }
}

