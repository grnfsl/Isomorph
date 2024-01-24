#ifndef GRAPH_H
#define GRAPH_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define TRIANGULAR_MATRIX   //non-triangular representation gives better performance but needs more space
#define ADJACENCY_LIST 1  // 0 for linked list representation and any other number for array representation

/* In is_isomorphic() function, each canonical_form() for each graph run in sequential
  if PARALLEL_CANON is false and parallel otherwise*/
#define PARALLEL_CANON 0

//----------------------------Graph adjacency list representation--------------------
#if ADJACENCY_LIST == 0
typedef struct Graph_node{
    unsigned int dest;
    struct Graph_node* next;
} Graph_node;

typedef struct Graph_list{
    Graph_node* head;
} Graph_list;

typedef struct Graph_l{
    unsigned long n;
    unsigned long e;
    Graph_list* array;
} Graph_l;

Graph_node* create_graph_node(unsigned int dest);

#else
typedef struct adjNode{
    unsigned long len;
    unsigned long d;
    unsigned int list[1];
}adjNode;

typedef struct Graph_l{
    unsigned long n;
    unsigned long e;
    adjNode **adjlist;
} Graph_l;

#endif

Graph_l* create_graph_l(const unsigned long n);
void add_edge_l(Graph_l *graph, unsigned int src, unsigned int dest);
void free_graph_l(Graph_l *graph);
void print_graph_l(const Graph_l *graph);


//----------------------------Graph matrix representation--------------------

typedef struct Graph_m{
    unsigned long n;
    unsigned long e;
    unsigned int size;
    int* matrix;
} Graph_m;

Graph_m* create_graph_m(unsigned long n);
bool add_edge_m(Graph_m *graph, unsigned int src, unsigned int dest);
void free_graph_m(Graph_m *graph);
void print_graph_m(Graph_m *graph);
Graph_l* mat_to_list(Graph_m *g);

void print_graph_by_ord(const Graph_m *g, int *a, bool full_m, bool triangular_m);

//----------------------------Graph struct hold both adjacency list and matrix representation--------------------

typedef struct Graph{
    unsigned long n;
    unsigned long e;
    Graph_l *graph_l;
    Graph_m *graph_m;
} Graph;

Graph* create_graph(int n);
void add_edge(Graph *g, unsigned int src, unsigned int dest);
void free_graph(Graph *g);

//---------------------------------Partition---------------------------------

typedef struct Cell{
    int *first_p, *last_p;
    bool counted, discrete;
    struct Cell* next;
} Cell;

typedef struct Partition{
    unsigned long n;
    Cell *head;
} Partition;

Partition* create_partition(unsigned long n, int *arr_v);
Partition* cpy_partition(Partition *partition);
Cell* create_cell(int *first_p, int *last_p, bool counted, bool discrete);
void split_cell(Cell *cell, int *degree); //The cell must already be sorted
void individualize_vertex(Cell *cell, int v); // 'v' is the vertex to be individualized from 'cell'
void sort_cell(Cell *cell);
void count_cells(Partition *partition, bool counted);
void free_partition(Partition *partition);
void print_cell(const Cell *cell);
void print_partition(const Partition *partition);

#endif // GRAPH_H

































