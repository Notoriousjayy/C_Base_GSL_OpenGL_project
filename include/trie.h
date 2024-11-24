//
// Created by jorda on 9/20/2024.
//

// trie.h
#ifndef TRIE_H
#define TRIE_H

#define ALPHABET_SIZE 26  // Since we are using only uppercase A-Z

// Structure to represent a node in the Trie
typedef struct TrieNode {
    struct TrieNode* children[ALPHABET_SIZE];
    int isEndOfWord;  // 1 if the node represents the end of a word
} TrieNode;

// Function prototypes
TrieNode* createTrieNode(void);
void insert_trie(TrieNode* root, const char* key);
int search_trie(TrieNode* root, const char* key);

#endif
