#ifndef ISOMORPH_H
#define ISOMORPH_H

#include <utils.h>
#include <automorph.h>
#include <time.h>

typedef struct Params{
    size_t n;
    int *degree;
    int *current_ord;
    int *fixed_ord;
    int *best_ord;
    bool best_exists;
    bool auto_found;
    Search_Node *top_search_node;
    int basis_ok;
    Search_Node *last_base_change;
} Params;

int compare_orders(Graph_m *g, int *best_ord, int *fixed_pts, unsigned int spt, unsigned int lpt);
void refine(Graph_l *g, Search_Node *s, Params *params);
void stabilise(Graph *g, Search_Node *s, Params *params);
bool is_isomorphic(Graph *g1, Graph *g2);
bool is_equal_ordering(Graph_m *graph_m1, Graph_m *graph_m2, int *ord1, int *ord2);
int* canonical_form(Graph *g);

Search_Node* create_search_node(int v);
void free_search_tree(Search_Node *s);
void print_search_node(Search_Node *s);
void print_perm(Permutation *p, int size);
#endif // ISOMORPH_H
