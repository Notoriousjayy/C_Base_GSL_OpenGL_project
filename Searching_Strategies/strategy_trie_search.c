//
// Created by jorda on 9/20/2024.
//

// strategy_trie_search.c
#include "../include/strategy_trie_search.h"
#include "../include/trie.h"


// Trie search strategy
void strategy_trie_search(TrieNode* root, const char* key, int* result) {
    int found = search_trie(root, key);
    *result = found ? 1 : -1;  // Return 1 if the key is found, otherwise -1
}

