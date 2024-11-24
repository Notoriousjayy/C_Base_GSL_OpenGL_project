//
// Created by jorda on 9/20/2024.
//

// trie.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/trie.h"

// Helper function to create a new Trie node
TrieNode* createTrieNode(void) {
    TrieNode* node = (TrieNode*)malloc(sizeof(TrieNode));
    if (node) {
        node->isEndOfWord = 0;
        for (int i = 0; i < ALPHABET_SIZE; i++) {
            node->children[i] = NULL;
        }
    }
    return node;
}

// Function to insert a word into the Trie
void insert_trie(TrieNode* root, const char* key) {
    TrieNode* pCrawl = root;
    for (int level = 0; level < strlen(key); level++) {
        int index = key[level] - 'A';  // Convert character to index
        if (!pCrawl->children[index]) {
            pCrawl->children[index] = createTrieNode();
        }
        pCrawl = pCrawl->children[index];
    }
    pCrawl->isEndOfWord = 1;  // Mark end of word
}

// Function to search for a word in the Trie
int search_trie(TrieNode* root, const char* key) {
    TrieNode* pCrawl = root;
    for (int level = 0; level < strlen(key); level++) {
        int index = key[level] - 'A';  // Convert character to index
        if (!pCrawl->children[index]) {
            return 0;  // Word not found
        }
        pCrawl = pCrawl->children[index];
    }
    return (pCrawl != NULL && pCrawl->isEndOfWord);
}
