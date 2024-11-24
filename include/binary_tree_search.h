//
// Created by jorda on 9/20/2024.
//

// binary_tree_search.h
#ifndef BINARY_TREE_SEARCH_H
#define BINARY_TREE_SEARCH_H

typedef struct TreeNode {
    int key;
    struct TreeNode* left;
    struct TreeNode* right;
} TreeNode;

TreeNode* insert(TreeNode* root, int key);
TreeNode* search(TreeNode* root, int key);
void inorder_print(TreeNode* root);

#endif
