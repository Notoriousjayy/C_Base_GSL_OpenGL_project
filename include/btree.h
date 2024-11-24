#ifndef BTREE_H
#define BTREE_H

#define MIN_DEGREE  3  // Minimum degree for the B-Tree

// Structure for a B-Tree node
typedef struct BTreeNode {
    int *keys;                // Array of keys
    struct BTreeNode **children;  // Array of child pointers
    int num_keys;             // Number of keys currently in node
    int leaf;                 // Boolean, true if leaf node
} BTreeNode;

// Structure for the B-Tree
typedef struct BTree {
    BTreeNode *root;          // Pointer to the root node
} BTree;

// Function declarations
BTreeNode* create_btree_node(int leaf);
BTree* create_btree();
void btree_insert(BTree *T, int k);
void btree_insert_nonfull(BTreeNode *x, int k);
void btree_split_child(BTreeNode *x, int i);
BTreeNode* btree_search(BTreeNode *x, int k);

#endif
