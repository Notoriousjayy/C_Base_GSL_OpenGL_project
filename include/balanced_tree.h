//
// Created by jorda on 9/20/2024.
//

// balanced_tree.h
#ifndef BALANCED_TREE_H
#define BALANCED_TREE_H

typedef struct TreeNode {
    int key;
    int balance;  // Balance factor: -1, 0, or +1
    struct TreeNode* left;
    struct TreeNode* right;
} TreeNode;

TreeNode* insert_balanced_tree(TreeNode* root, int key);
TreeNode* search_balanced_tree(TreeNode* root, int key);
void inorder_print(TreeNode* root);

#endif
