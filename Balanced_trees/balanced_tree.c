//
// Created by jorda on 9/20/2024.
//

// balanced_tree.c
#include "../include/balanced_tree.h"
#include <stdio.h>
#include <stdlib.h>

// Helper function to create a new node
TreeNode* create_balanced_tree_node(int key) {
    TreeNode* new_node = (TreeNode*)malloc(sizeof(TreeNode));
    if (new_node == NULL) {
        printf("Memory allocation failed.\n");
        return NULL;
    }
    new_node->key = key;
    new_node->balance = 0;
    new_node->left = NULL;
    new_node->right = NULL;
    return new_node;
}

// Function to search for a key in the balanced tree
TreeNode* search_balanced_tree(TreeNode* root, int key) {
    if (root == NULL || root->key == key) {
        return root;  // Key found or tree is empty
    }

    if (key < root->key) {
        return search_balanced_tree(root->left, key);
    }

    return search_balanced_tree(root->right, key);
}

// Function to perform a right rotation
TreeNode* right_rotate(TreeNode* y) {
    TreeNode* x = y->left;
    TreeNode* T2 = x->right;

    // Perform rotation
    x->right = y;
    y->left = T2;

    return x;  // New root
}

// Function to perform a left rotation
TreeNode* left_rotate(TreeNode* x) {
    TreeNode* y = x->right;
    TreeNode* T2 = y->left;

    // Perform rotation
    y->left = x;
    x->right = T2;

    return y;  // New root
}

// Function to insert a key into the balanced tree and maintain balance
TreeNode* insert_balanced_tree(TreeNode* root, int key) {
    // Perform the normal BST insertion
    if (root == NULL) {
        return create_balanced_tree_node(key);
    }

    if (key < root->key) {
        root->left = insert_balanced_tree(root->left, key);
    } else if (key > root->key) {
        root->right = insert_balanced_tree(root->right, key);
    } else {
        return root;  // Duplicate keys are not allowed
    }

    // Balance the tree (this is a simplified example; full balance code omitted)
    // Check balance and perform rotations if necessary
    if (root->balance == -2) {
        if (key < root->left->key) {
            return right_rotate(root);  // Left-Left case
        } else {
            root->left = left_rotate(root->left);  // Left-Right case
            return right_rotate(root);
        }
    } else if (root->balance == 2) {
        if (key > root->right->key) {
            return left_rotate(root);  // Right-Right case
        } else {
            root->right = right_rotate(root->right);  // Right-Left case
            return left_rotate(root);
        }
    }

    return root;
}

// Inorder traversal to print the tree
void inorder_balanced_tree_print(TreeNode* root) {
    if (root != NULL) {
        inorder_print(root->left);
        printf("%d ", root->key);
        inorder_print(root->right);
    }
}
