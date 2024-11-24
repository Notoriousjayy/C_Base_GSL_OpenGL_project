//
// Created by jorda on 9/20/2024.
//

// strategy_balanced_tree_search.c
#include <stddef.h>
#include "../include/strategy_balanced_tree_search.h"
#include "../include/balanced_tree.h"

// Balanced tree search strategy
void strategy_balanced_tree_search(TreeNode* root, int key, int* result) {
    TreeNode* node = search_balanced_tree(root, key);
    if (node != NULL) {
        *result = node->key;  // Return the key if found
    } else {
        *result = -1;  // Return -1 if the key is not found
    }
}
