#include <stdio.h>
#include <stdlib.h>
#include "../include/btree.h"

// Create a new B-Tree node
BTreeNode* create_btree_node(int leaf) {
    BTreeNode *node = (BTreeNode *)malloc(sizeof(BTreeNode));
    node->leaf = leaf;
    node->keys = (int *)malloc((2 * MIN_DEGREE - 1) * sizeof(int));
    node->children = (BTreeNode **)malloc(2 * MIN_DEGREE * sizeof(BTreeNode *));
    node->num_keys = 0;
    return node;
}

// Create an empty B-Tree
BTree* create_btree() {
    BTree *tree = (BTree *)malloc(sizeof(BTree));
    tree->root = create_btree_node(1); // Root is a leaf initially
    return tree;
}

// Function to search for a key in the B-Tree
BTreeNode* btree_search(BTreeNode *x, int k) {
    int i = 0;

    // Find the first key greater than or equal to k
    while (i < x->num_keys && k > x->keys[i])
        i++;

    // If the found key is equal to k, return the current node
    if (i < x->num_keys && k == x->keys[i])
        return x;

    // If it's a leaf, key is not found
    if (x->leaf) {
        return NULL;
    } else {
        // Otherwise, search in the appropriate child node
        return btree_search(x->children[i], k);
    }
}


// Function to insert a key into the B-Tree
void btree_insert(BTree *T, int k) {
    BTreeNode *root = T->root;

    // If root is full, create a new root
    if (root->num_keys == 2 * MIN_DEGREE - 1) {
        BTreeNode *s = create_btree_node(0);
        s->children[0] = root;
        T->root = s;
        btree_split_child(s, 0);
        btree_insert_nonfull(s, k);
    } else {
        btree_insert_nonfull(root, k);
    }
}

// Insert a key into a non-full node
void btree_insert_nonfull(BTreeNode *x, int k) {
    int i = x->num_keys - 1;

    if (x->leaf) {
        // Shift keys to the right to make room for the new key
        while (i >= 0 && k < x->keys[i]) {
            x->keys[i + 1] = x->keys[i];
            i--;
        }
        x->keys[i + 1] = k;
        x->num_keys++;
    } else {
        // Find the child to recurse into
        while (i >= 0 && k < x->keys[i])
            i--;
        i++;

        // If the found child is full, split it
        if (x->children[i]->num_keys == 2 * MIN_DEGREE - 1) {
            btree_split_child(x, i);
            if (k > x->keys[i])
                i++;
        }
        btree_insert_nonfull(x->children[i], k);
    }
}

// Split a full child node
void btree_split_child(BTreeNode *x, int i) {
    BTreeNode *y = x->children[i];
    BTreeNode *z = create_btree_node(y->leaf);
    z->num_keys = MIN_DEGREE - 1;

    // Copy the last MIN_DEGREE - 1 keys from y to z
    for (int j = 0; j < MIN_DEGREE - 1; j++)
        z->keys[j] = y->keys[j + MIN_DEGREE];

    // If y is not a leaf, copy the last MIN_DEGREE children from y to z
    if (!y->leaf) {
        for (int j = 0; j < MIN_DEGREE; j++)
            z->children[j] = y->children[j + MIN_DEGREE];
    }

    y->num_keys = MIN_DEGREE - 1;

    // Shift children of x to make room for the new child
    for (int j = x->num_keys; j >= i + 1; j--)
        x->children[j + 1] = x->children[j];

    // Link z to x
    x->children[i + 1] = z;

    // Shift keys of x to make room for the new key
    for (int j = x->num_keys - 1; j >= i; j--)
        x->keys[j + 1] = x->keys[j];

    // Insert the middle key of y into x
    x->keys[i] = y->keys[MIN_DEGREE - 1];
    x->num_keys++;
}
