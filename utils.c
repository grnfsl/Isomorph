
/*  Date for last update: 09/01/2021
 *  Author: Goran Faisal; grnhry@gmail.com
 *
 *  g_V[] is storing the current ordering of the vertices of the graph.
 *  The way cells of the partition is represented is that each cell contains
 *  two pointer to array g_V[]. First pointer refer to the first element
 *  of the cell in g_V[] and second pointer refer to the last element of the
 *  cell in g_V[].
 */

#include <utils.h>
#include <time.h>
void qsortcell(int *arr, size_t n);

//----------------------------Graph adjacency list representation--------------------

#if ADJACENCY_LIST == 0     //Graph adjacency list representation using linked list
Graph_node* create_graph_node(unsigned int dest)
{
    Graph_node* new_node = (Graph_node*) malloc(sizeof (Graph_node));
    new_node->dest = dest;
    new_node->next = NULL;
    return new_node;
}

Graph_l* create_graph_l(size_t v)
{
    Graph_l* g = (Graph_l*) malloc(sizeof (Graph_l));

    g->n = v;
    g->e = 0;
    g->array = (Graph_list*) malloc(v * sizeof (Graph_list));

    size_t i;
    for(i = 0; i < v; ++i)
        g->array[i].head = NULL;

    return g;
}

void add_edge_l(Graph_l* graph, unsigned int src, unsigned int dest)
{
    Graph_node *new_node;

    new_node = create_graph_node(dest);
    new_node->next = graph->array[src].head;
    graph->array[src].head = new_node;

    //Since graph is undirected, an edge added from dest to src
    new_node = create_graph_node(src);
    new_node->next = graph->array[dest].head;
    graph->array[dest].head = new_node;
    ++graph->e;
}

void free_graph_l(Graph_l* g)
{
    size_t i;
    for(i = 0; i < g->n; ++i){
        Graph_node* head = g->array[i].head;
        Graph_node* tmp;
        while (head) {
            tmp = head;
            head = head->next;
            free(tmp);
        }
    }
    free(g->array);
    free(g);
}

void print_graph_l(const Graph_l* g)
{
    size_t v;
    for(v = 0; v < g->n; ++v){
        Graph_node* currenct_node = g->array[v].head;
        printf("%ld: ", v);
        while (currenct_node) {
            printf("%u ", currenct_node->dest);
            currenct_node = currenct_node->next;
        }
        printf("\n");
    }
}
#else       //Graph adjacency list representation using array

Graph_l* create_graph_l(size_t n)
{
    Graph_l* g = malloc(sizeof (Graph_l));
    g->adjlist = malloc(n * sizeof (adjNode*));
    g->n = n;
    g->e = 0;
    for(size_t i = 0; i < n; ++i){
        g->adjlist[i] = malloc(sizeof (adjNode));
        g->adjlist[i]->len = 1;
        g->adjlist[i]->d = 0;
    }
    return g;
}

bool is_edge_exists(adjNode **l, unsigned int src, unsigned int dest)
{
    size_t i;
    for(i = 0; i < l[src]->d; ++i)
        if(l[src]->list[i] == dest) return true;
    return false;
}

void add_edge_l(Graph_l* g, unsigned int src, unsigned int dest)
{
//    if(is_edge_exists(g->adjlist, src, dest)) return;

    if(g->adjlist[src]->d >= g->adjlist[src]->len){
        g->adjlist[src]->len *= 2;
        g->adjlist[src] = realloc(g->adjlist[src], sizeof(adjNode) + sizeof(int) * (g->adjlist[src]->len - 1));
    }
    g->adjlist[src]->list[g->adjlist[src]->d] = dest;
    ++g->adjlist[src]->d;

    if(g->adjlist[dest]->d >= g->adjlist[dest]->len){
        g->adjlist[dest]->len *= 2;
        g->adjlist[dest] = realloc(g->adjlist[dest], sizeof(adjNode) + sizeof(int) * (g->adjlist[dest]->len - 1));
    }
    g->adjlist[dest]->list[g->adjlist[dest]->d] = src;
    ++g->adjlist[dest]->d;

    ++g->e;
}

void free_graph_l(Graph_l *g)
{
    size_t i;
    for(i = 0; i < g->n; ++i)
        free(g->adjlist[i]);
    free(g->adjlist);
    free(g);
}

void print_graph_l(const Graph_l* g)
{
    size_t size;
    for(size_t i = 0; i < g->n; ++i){
        size = g->adjlist[i]->d;
        printf("%ld: ", i);
        for(size_t j = 0; j < size; ++j)
            printf("%u ", g->adjlist[i]->list[j]);
        printf("\n");
    }
}

#endif

//--------------------------------Graph matrix representation------------------------------

#ifdef TRIANGULAR_MATRIX

Graph_m* create_graph_m(size_t n)
//create triangular matrix because it is undirected graph
{
    size_t i;

    Graph_m* g = (Graph_m*)malloc(sizeof (Graph_m));
    g->n = n;
    g->e = 0;
    g->size = ((n-1) * n)/2;
    g->matrix = (int*)malloc(g->size * sizeof (int));

    for(i = 0; i < g->size; ++i)
        g->matrix[i] = 0;

    return g;
}

bool add_edge_m(Graph_m *g, unsigned int src, unsigned int dest)
//it consider that graph is undirected simple graph
{   //ignore when src == dest because diagonal should be 0, graph is simple undirected
    int v;
    if(src < dest)
        v = src*g->n + dest - ((src+1)*(src+2))/2;
    else if(src > dest)
        v = dest*g->n + src - ((dest+1)*(dest+2))/2;

    if(g->matrix[v] != 1){
        g->matrix[v] = 1;
        ++g->e;
        return true;
    }
    return false;
}

void print_graph_m(Graph_m* g)
{
    size_t i, j;
    for(i = 0; i < g->n; ++i){
        for(j = 0; j < g->n; ++j){
            if(i < j)
                printf("%d ", g->matrix[i*g->n + j - ((i+1)*(i+2))/2]);
            else if(i > j)
                printf("%d ", g->matrix[j*g->n + i - ((j+1)*(j+2))/2]);
            else
                printf("%d ", 0);
        }
        printf("\n");
    }
}

Graph_l* mat_to_list(Graph_m *g)
{
    size_t r, c, f;
    Graph_l *graph = create_graph_l(g->n);
    for(r = 0; r < g->n-1; ++r)
        for(c = r+1; c < g->n; ++c){
            f = r*g->n + c - ((r+1)*(r+2))/2;
            if(g->matrix[f] == 1)
                add_edge_l(graph, r, c);
        }
    return graph;
}

//print a 1d array representing an upper triagular adjacency matrix by ordering 'a'
void print_graph_by_ord(const Graph_m *g, int *a, bool full_m, bool triangular_m)
{
    size_t r, c, f;
    if(full_m){
        for(r = 0; r < g->n; ++r){
            for(c = 0; c < g->n; ++c){
                if(a[r] < a[c])
                    printf("%d ", g->matrix[a[r] * g->n + a[c] - ((a[r]+1)*(a[r]+2))/2]);
                else if(a[r] > a[c])
                    printf("%d ", g->matrix[a[c] * g->n + a[r] - ((a[c]+1)*(a[c]+2))/2]);
                else
                    printf("%d ", 0);
            }
            printf("\n");
        }
    }
    if(triangular_m){
        for(c = 1; c < g->n; ++c){
            for(r = 0; r < c; ++r){
                f = (a[r] < a[c]) ? a[r] * g->n + a[c] - ((a[r]+1)*(a[r]+2))/2
                                  : a[c] * g->n + a[r] - ((a[c]+1)*(a[c]+2))/2;
                printf("%d ", g->matrix[f]);
            }
        }
    }
    printf("\n");
}

#else

Graph_m* create_graph_m(size_t n)
{
    size_t i;

    Graph_m* g = (Graph_m*)malloc(sizeof (Graph_m));
    g->n = n;
    g->e = 0;
    g->matrix = (int*)malloc(n * n * sizeof (int));

    for(i = 0; i < n*n; ++i)
        g->matrix[i] = 0;

    return g;
}

bool add_edge_m(Graph_m *g, unsigned int src, unsigned int dest)
//it consider simple undirected graph
{
    if(g->matrix[src*g->n+dest] != 1){
        g->matrix[src*g->n+dest] = 1;
        g->matrix[dest*g->n+src] = 1;
        ++g->e;
        return true;
    }
    return false;
}

void print_graph_m(Graph_m* g)
{
    size_t i, j;
    for(i = 0; i < g->n; ++i){
        for(j = 0; j < g->n; ++j)
            printf("%d ", g->matrix[i*g->n+j]);
        printf("\n");
    }
}

Graph_l* mat_to_list(Graph_m *g)
{
    size_t r, c;
    Graph_l *graph = create_graph_l(g->n);
    for(r = 0; r < g->n-1; ++r)
        for(c = r+1; c < g->n; ++c){
            if(g->matrix[r*g->n+c] == 1)
                add_edge_l(graph, r, c);
        }
    return graph;
}

void print_graph_by_ord(const Graph_m *g, int *a, bool full_m, bool triangular_m)
{
    size_t r, c, f;
    if(full_m){
        for(r = 0; r < g->n; ++r){
            for(c = 0; c < g->n; ++c){
                printf("%d ", g->matrix[a[r] * g->n + a[c]]);
            }
            printf("\n");
        }
    }
    if(triangular_m){
        for(c = 1; c < g->n; ++c){
            for(r = 0; r < c; ++r){
                f = (a[r] < a[c]) ? a[r] * g->n + a[c] : a[c] * g->n + a[r];
                printf("%d ", g->matrix[f]);
            }
        }
    }
    printf("\n");
}

#endif

void free_graph_m(Graph_m* g)
{
    free(g->matrix);
    free(g);
}

//----------------A container Graph for both matrix representation and djacency list---------------------

Graph* create_graph(int n)
{
    Graph *g = malloc(sizeof (Graph));
    g->n = n;
    g->e = 0;
    g->graph_m = create_graph_m(n);
    g->graph_l = create_graph_l(n);

    return g;
}

void add_edge(Graph *g, unsigned int src, unsigned int dest)
{
    if(add_edge_m(g->graph_m, src, dest)){
        add_edge_l(g->graph_l, src, dest);
        g->e = g->graph_l->e;
    }
}

void free_graph(Graph *g)
{
    free_graph_m(g->graph_m);
    free_graph_l(g->graph_l);
    free(g);
}

//---------------------------------Partition---------------------------------

Partition* create_partition(size_t n, int *arr_v)
{
    Partition *partition = (Partition*) malloc(sizeof (Partition));

    partition->n = n;
    partition->head = create_cell(arr_v, arr_v + n - 1, false, false);

    return partition;
}

Partition* cpy_partition(Partition *partition)
{
    Partition *new_partition;
    Cell *cell, *new_cell;
    int *f_i, *l_i;

    new_partition = create_partition(partition->n, partition->head->first_p);

    cell = partition->head;
    new_cell = new_partition->head;

    new_cell->first_p = cell->first_p;
    new_cell->last_p = cell->last_p;

    cell = cell->next;
    while (cell) {
        f_i = cell->first_p;
        l_i = cell->last_p;
        new_cell->next = create_cell(f_i, l_i, cell->counted, cell->discrete);
        new_cell = new_cell->next;
        cell = cell->next;
    }
    return new_partition;
}

Cell* create_cell(int *first_p, int *last_p, bool counted, bool discrete)
{
    Cell *cell = (Cell*) malloc(sizeof (Cell));

    cell->counted = counted;
    cell->discrete = discrete;
    cell->first_p = first_p;
    cell->last_p = last_p;
    cell->next = NULL;

    return cell;
}

void split_cell(Cell *cell, int *degree)
//The cell must already be sorted
{
    Cell *new_cell, *prev_cell;
    int *i, *p_i, *f_i, *l_i;

    f_i = p_i = cell->first_p;
    l_i = cell->last_p;
    new_cell = cell;

    for(i = f_i+1; i <= l_i; ++i){
        if(degree[*i] != degree[*p_i]){ //if degree is not constant, split
            new_cell->last_p = p_i;
            new_cell->discrete = (new_cell->last_p == new_cell->first_p);
            prev_cell = new_cell;

            new_cell = create_cell(i, i, false, false);
            new_cell->next = prev_cell->next;
            prev_cell->next = new_cell;
        }
        p_i = i;
    }
    new_cell->last_p = p_i;
    new_cell->discrete = (new_cell->last_p == new_cell->first_p);
    //if original cell splited, then mark first cell uncounted
    //which is the original cell is used again as first cell
    if(cell != new_cell)    cell->counted = false;
}

void individualize_vertex(Cell *cell, int v) // 'v' is the vertex to be individualized from 'cell'
//used for stablising vertices
{
    Cell *new_cell;
    int *i, *p, tmp;

    for(i = cell->first_p; i <= cell->last_p; ++i){
        if(*i == v){
            //bring the individualize vertex to the beginning of the array
            tmp = *cell->first_p;
            *cell->first_p = *i;
            *i = tmp;

            p = cell->last_p;       //pointer for last vertex of the unsplitted cell
            cell->last_p = cell->first_p; //because the first cell contains one vertex
            cell->counted = false;
            cell->discrete = true;

            new_cell = create_cell(cell->last_p+1, p, false, false);
            new_cell->next = cell->next;
            cell->next = new_cell;
            break;
        }
    }
}

void count_cells(Partition *partition, bool counted)
//mark all the cells as 'counted'
{
    Cell *head_cell;

    head_cell = partition->head;
    while(head_cell){
        head_cell->counted = counted;
        head_cell = head_cell->next;
    }
}

void free_partition(Partition *partition)
{
    Cell *cell;
    while(partition->head){
        cell = partition->head;
        partition->head = partition->head->next;
        free(cell);
    }
    free(partition);
}

void print_cell(const Cell *cell)
{
    if(cell == NULL) printf("{NULL}");
    printf("{ ");
    int *i;
    for(i = cell->first_p; i <= cell->last_p; ++i)
        printf("%d ", *i);
    printf("}");
}

void print_partition(const Partition *partition)
{
    printf("( ");
    Cell *head_cell;

    head_cell = partition->head;
    while(head_cell){
        print_cell(head_cell);
        head_cell = head_cell->next;
    }
    printf(" )\n");
}















































