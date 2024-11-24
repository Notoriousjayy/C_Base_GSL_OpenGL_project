#include <stdio.h>
#include "include/btree.h"
#include "include/strategy_btree_search.h"
#include "include/search_strategy.h"

int main() {
    BTree* tree = create_btree();

    btree_insert(tree, 10);
    btree_insert(tree, 20);
    btree_insert(tree, 5);
    btree_insert(tree, 6);
    btree_insert(tree, 12);
    btree_insert(tree, 30);
    btree_insert(tree, 7);
    btree_insert(tree, 17);

    SearchContext search_context;
    set_search_strategy(&search_context, strategy_btree_search);

    int key = 10;
    int result;
    execute_search_strategy(&search_context, tree->root, &key, &result);  // Pass address of key

    key = 15;
    execute_search_strategy(&search_context, tree->root, &key, &result);  // Pass address of key

    return 0;
}
