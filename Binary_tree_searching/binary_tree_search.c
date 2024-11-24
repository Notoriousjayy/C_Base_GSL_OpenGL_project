//
// Created by jorda on 9/20/2024.
//

// binary_tree_search.c
#include "../include/binary_tree_search.h"
#include <stdio.h>
#include <stdlib.h>

// Function to create a new node in the tree
TreeNode* create_node(int key) {
    TreeNode* new_node = (TreeNode*)malloc(sizeof(TreeNode));
    if (new_node == NULL) {
        printf("Memory allocation failed.\n");
        return NULL;
    }
    new_node->key = key;
    new_node->left = NULL;
    new_node->right = NULL;
    return new_node;
}

// Function to insert a new key into the tree
TreeNode* insert(TreeNode* root, int key) {
    if (root == NULL) {
        return create_node(key);
    }

    if (key < root->key) {
        root->left = insert(root->left, key);
    } else if (key > root->key) {
        root->right = insert(root->right, key);
    }

    return root;
}

// Function to search for a key in the binary tree
TreeNode* search(TreeNode* root, int key) {
    if (root == NULL || root->key == key) {
        return root;  // If key is present at root or tree is empty
    }

    if (key < root->key) {
        return search(root->left, key);  // Key is smaller than root's key
    }

    return search(root->right, key);  // Key is greater than root's key
}

// Helper function to perform inorder traversal of the tree
void inorder_print(TreeNode* root) {
    if (root != NULL) {
        inorder_print(root->left);
        printf("%d ", root->key);
        inorder_print(root->right);
    }
}
