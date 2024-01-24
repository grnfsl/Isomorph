/* 
 * This file contains all neccesary functions for isomorphism program
 * is_isomorphic() function may work as interface for isomorphism program
 */

#include <isomorph.h>
#include <limits.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define BETTER 1
#define EQUAL 0
#define WORSE -1
#define NOT_EXIST -1

bool qsortcell(int *arr, int *degree, unsigned long n);

/******************************************************************************
 *  next_cell(partition, prev) returns next uncounted cell and assign prev
 *  to the cell before the returned uncounted cell.
 *  prev will be used in refine() function for removing discrte counted cell
 *  to avoid going through whole partition.
 *****************************************************************************/

Cell* next_cell(Partition *partition, Cell **prev)
//get the next uncounted cell
{
    Cell *cell;

    cell = partition->head;
    while (cell) {
        if(cell->counted == false) return cell;
        *prev = cell;
        cell = cell->next;
    }
    return NULL;
}

/******************************************************************************
 *  compare(g, best) compares best ordering found so far with the currrent
 *  obtained ordering.
 *  If equal, return 0. If better found, return 1. If worse, return -1.
 *  For accessing the 1D array that represents graph in upper triagular adjacency
 *  matrix, the formula {r * n + c - ((r + 1) * (r + 2)) / 2} is used, where
 *  r < c. r is the row and c is the column of the matrix of the graph, and n
 *  is the number of vertices of the graph.
 *****************************************************************************/
int compare_orderings(Graph_m *g, int *best_ord, int *fixed_ord, unsigned int spt, unsigned int lpt)
#ifdef TRIANGULAR_MATRIX
{
    unsigned int f, b, c, r;

    for(r = spt+1; r < lpt; ++r)
        for(c = 0; c < r; ++c){
            /*  if..else is used to switch raw and column because in some iterations
             *  column > raw which formula only accept when raw < column.
            */
            f = (fixed_ord[r] < fixed_ord[c]) ? fixed_ord[r] * g->n + fixed_ord[c] - ((fixed_ord[r]+1)*(fixed_ord[r]+2))/2
                    : fixed_ord[c] * g->n + fixed_ord[r] - ((fixed_ord[c]+1)*(fixed_ord[c]+2))/2;
            b = (best_ord[r] < best_ord[c]) ? best_ord[r] * g->n + best_ord[c] - ((best_ord[r]+1)*(best_ord[r]+2))/2
                    : best_ord[c] * g->n + best_ord[r] - ((best_ord[c]+1)*(best_ord[c]+2))/2;

            if(g->matrix[f] < g->matrix[b]) return BETTER;
            if(g->matrix[f] > g->matrix[b]) return WORSE;
        }
    return EQUAL;
}
#else
{
    unsigned int r, c;

    for(r = spt+1; r < lpt; ++r)
        for(c = 0; c < r; ++c){
            if(g->matrix[fixed_ord[r]*g->n+fixed_ord[c]] < g->matrix[best_ord[r]*g->n+best_ord[c]]) return BETTER;
            if(g->matrix[fixed_ord[r]*g->n+fixed_ord[c]] > g->matrix[best_ord[r]*g->n+best_ord[c]]) return WORSE;
        }
    return EQUAL;
}
#endif

/******************************************************************************
 *  refine(...) does the refinement of the partition, so when function returns,
 *  the partition is equitable.
 *****************************************************************************/

void refine(Graph_l *g, Search_Node *s, Params *params)
{
    count_cells(s->p, false);

    Cell *cell, *cell_p, *tmp, *prev;
    int *i;
    bool is_cell_discrete;
    unsigned long d;

    cell = s->p->head;
    while (cell) {
#if ADJACENCY_LIST == 0     //Graph adjacency list representation using linked list
        Graph_node *node;
        if(cell->discrete){

            node = g->array[*cell->first_p].head;
            while(node){
                ++params->degree[node->dest];
                node = node->next;
            }
            params->fixed_ord[s->nfixed] = *cell->first_p;
            ++s->nfixed;
            if(cell->first_p == cell->last_p){ //if discrete cell contains of one vertex
                if(cell == s->p->head)
                    s->p->head = cell->next;
                else
                    prev->next = cell->next;
                free(cell);
            }
            else ++cell->first_p; //else just shrink the size of cell by one vertex
        }
        else {
            //First count the adjacencies to cell 'cell'
            for(i = cell->first_p; i <= cell->last_p; ++i){
                node = g->array[*i].head;
                while(node){
                    ++params->degree[node->dest];
                    node = node->next;
                }
            }
            cell->counted = true;
        }
#else   //Graph adjacency list representation using array
        unsigned long j;
        if(cell->discrete){
            for(j = 0; j < g->adjlist[*cell->first_p]->d; ++j)
                ++params->degree[g->adjlist[*cell->first_p]->list[j]];

            params->fixed_ord[s->nfixed] = *cell->first_p;
            ++s->nfixed;
            if(cell->first_p == cell->last_p){ //if discrete cell contains of one vertex
                if(cell == s->p->head)
                    s->p->head = cell->next;
                else
                    prev->next = cell->next;
                free(cell);
            }
            else ++cell->first_p; //else just shrink the size of cell by one vertex
        }
        else {
            //First count the adjacencies to cell 'cell'
            for(i = cell->first_p; i <= cell->last_p; ++i)
                for(j = 0; j < g->adjlist[*i]->d; ++j)
                    ++params->degree[g->adjlist[*i]->list[j]];
            cell->counted = true;
        }
#endif
        cell_p = s->p->head;
        d = s->nfixed;
        unsigned long cell_size;
        while(cell_p){
            tmp = cell_p->next;
            cell_size = (int)((cell_p->last_p) - (cell_p->first_p))+1;
            if(cell_p->discrete)
                d += cell_size;
            else if(!cell_p->discrete){
                is_cell_discrete = qsortcell(cell_p->first_p, params->degree, cell_size);
                if(is_cell_discrete){
                    cell_p->discrete = true;
                    cell_p->counted = false;
                    d += cell_size;
                }
                else split_cell(cell_p, params->degree);
            }
            cell_p = tmp; //To avoid going through new generated cells
        }

        if(d == g->n){
            while(s->p->head){
                for(i = s->p->head->first_p; i <= s->p->head->last_p; ++i){
                    params->fixed_ord[s->nfixed] = *i;
                    ++s->nfixed;
                }
                tmp = s->p->head;
                s->p->head = s->p->head->next;
                free(tmp);
            }
            return;
        }
        cell = next_cell(s->p, &prev); //next uncounted cell
    }
}

/******************************************************************************
 Changes the basis of the automprhism group so as to become the same as the
 basis of the partitions.
 *****************************************************************************/

void change_base(int d, Search_Node *last_base_change)
{
    Group *grp;
    Search_Node *s;
    unsigned long n;

    n = last_base_change->p->n;
    s = last_base_change;
    grp = s->grp;
    if(grp == NULL || grp->st_grp_u == NULL) return;

    s = s->next;
    while (s) {
        if(s->depth <= d){
            Cell *c = s->p->head;
            for(int *i = c->first_p; i <= c->last_p; ++i){
                s->cell_orbits[*i] = -1;
            }
        }
        s = s->next;
    }

    free_cosetreps(grp->coset, n);
    free_cosetreps(grp->coset_inv, n);

    grp->coset = create_cosetreps(n);
    grp->coset_inv = create_cosetreps(n);
    grp->orbit_npts = 0;

    free_grps(grp->st_grp_u, n);
    grp->st_grp_u = NULL;

    Generator *tmp, *gen = grp->gens;
    grp->gens = NULL;
    while (gen) {
        add_gen(gen->perm, last_base_change);
        tmp = gen;
        gen = gen->next;
        free(tmp);
    }
}

/******************************************************************************
 Selects the next inequivalent vertex so as to be individualized
 *****************************************************************************/

int next_inequiv_vertex(const Cell *c, int *cell_orbits, int n)
{
    int *i;
    for(i = c->first_p; i <= c->last_p; ++i)
        if(cell_orbits[orbit_rep(*i, cell_orbits)] + n >= 0) return *i;
    return -1;
}

void set_on_best_path(Search_Node *s, bool val)
{
    Search_Node *tmp = s;
    while (tmp) {
        tmp->on_best_path = val;
        tmp = tmp->next;
    }
}

/******************************************************************************
 * stabilise(...) generates search tree and leafs of the tree are discrete
 * partitions which they correspond to the orderings of the vertices in the
 * adjacency matrix.
 *****************************************************************************/

void stabilise(Graph *g, Search_Node *s, Params *params)
{
    unsigned long j, m;
    int result;

    m = s->nfixed;
    s->on_best_path = false;
    refine(g->graph_l, s, params);

    if(params->best_exists)
        result = compare_orderings(g->graph_m, params->best_ord, params->fixed_ord, m, s->nfixed);
    else result = EQUAL;

    if(s->nfixed == params->n){ //if the partition is discrete
        if(params->best_exists){
            if(result == EQUAL){ 
                Permutation *y = malloc(params->n * sizeof (int));
                for(j = 0; j < params->n; ++j)
                    y[params->fixed_ord[j]] = params->best_ord[j];
                add_gen(y, params->top_search_node);
                params->auto_found = true;
            }
            else if(result == BETTER){ 
                memcpy(params->best_ord, params->fixed_ord, params->n * sizeof (int));
                set_on_best_path(params->top_search_node, true);
            } //if worse ignore
        }else{
            memcpy(params->best_ord, params->fixed_ord, params->n * sizeof (int));
            params->best_exists = true;
            set_on_best_path(params->top_search_node, true);
        }
    }
    else{
        if(result != WORSE){ 
            if(result == BETTER)
                params->best_exists = false; 

            Search_Node *s_u;
            Partition *upartition;
            Cell *cell;

            cell = s->p->head;

            if(s->next == NULL)
                s->next = create_search_node(params->n);
            s_u = s->next;

            int target_v = *s->p->head->first_p;
            while (target_v != NOT_EXIST) {
                upartition = cpy_partition(s->p);
                individualize_vertex(upartition->head, target_v);

                s->u = target_v;
                s_u->p = upartition;
                s_u->nfixed = s->nfixed;
                s_u->depth = s->depth+1;

                stabilise(g, s_u, params);

                if(params->auto_found){
                    if(!s->on_best_path) {
                        free_partition(upartition);
                        break;
                    }
                    params->auto_found = false;
                }

                if(s->depth > params->basis_ok)
                    change_base(s->depth, params->last_base_change);
                params->basis_ok = s->depth;
                params->last_base_change = s;

                s->cell_orbits[orbit_rep(target_v, s->cell_orbits)] -= params->n;

                target_v = next_inequiv_vertex(cell, s->cell_orbits, params->n);
                free_partition(upartition);
            }
        }
    }
    for(unsigned long t = m; t < s->nfixed; ++t)
        params->degree[params->fixed_ord[t]] = 0;
    s->u = INACTIVE; //mark this search node inactive
}

/******************************************************************************
 *  is_equal_ordering(...) compares two orderings.
 *****************************************************************************/

bool is_equal_orderings(Graph_m *graph_m1, Graph_m *graph_m2, int *ord1, int *ord2)
#ifdef TRIANGULAR_MATRIX
{
    unsigned long f, b, c, r;

    for(r = 1; r < graph_m1->n; ++r)
        for(c = 0; c < r; ++c){
            /*  if..else is used to switch raw and column because in some iterations
             *  column > raw which formula only accept when raw < column.
            */
            f = (ord1[r] < ord1[c]) ? ord1[r] * graph_m1->n + ord1[c] - ((ord1[r]+1)*(ord1[r]+2))/2
                    : ord1[c] * graph_m1->n + ord1[r] - ((ord1[c]+1)*(ord1[c]+2))/2;
            b = (ord2[r] < ord2[c]) ? ord2[r] * graph_m1->n + ord2[c] - ((ord2[r]+1)*(ord2[r]+2))/2
                    : ord2[c] * graph_m1->n + ord2[r] - ((ord2[c]+1)*(ord2[c]+2))/2;

            if(graph_m1->matrix[f] != graph_m2->matrix[b]) return false;
        }
    return true;
}
#else
{
    unsigned int r, c;

    for(r = 1; r < graph_m1->n; ++r)
        for(c = 0; c < r; ++c){
            if(graph_m1->matrix[ord1[r]*graph_m1->n+ord1[c]] != graph_m2->matrix[ord2[r]*graph_m2->n+ord2[c]])
                return false;
        }
    return true;
}
#endif

/******************************************************************************
 * canonical(...) returns the canonical ordering of a graph.
 *****************************************************************************/

int* canonical_form(Graph *g)
{
    int i, n;

    n = g->n;

    Params params;

    params.n = n;
    params.degree = (int*) calloc(n, sizeof (int)); //All vertices' degree must be initialized to 0
    params.current_ord = (int*) malloc(n * sizeof (int));
    params.best_ord = (int*) malloc(n * sizeof (int));
    params.fixed_ord = (int*) malloc(n * sizeof (int));
    params.basis_ok = INT_MAX;
    params.auto_found = false;
    params.best_exists = false;

    for(i = 0; i < n; ++i) params.current_ord[i] = i;

    Partition *p = create_partition(n, params.current_ord);
    Group *grp = identity_grp(n);
    params.top_search_node = create_search_node(n);
    params.top_search_node->p = p;
    params.top_search_node->grp = grp;

    stabilise(g, params.top_search_node, &params);

    free(params.degree);
    free(params.current_ord);
    free_search_tree(params.top_search_node);
    free(params.fixed_ord);

    return params.best_ord;
}

/******************************************************************************
 * is_isomorphic(...) checks if two given graphs are isomorphic.
 *****************************************************************************/

bool is_isomorphic(Graph *g1, Graph *g2)
{
    if(g1->n != g2->n || g1->e != g2->e) return false; 
    if(g1->n == 0 || g1->e == 0) return true;

    bool b;
    int *best1, *best2;   //best orderings for g1 and g2

#if PARALLEL_CANON == 1
    double start = omp_get_wtime();

    #pragma omp parallel sections
    {
       #pragma omp section
       {
         best1 = canonical_form(g1);
       }
       #pragma omp section
       {
          best2 = canonical_form(g2);
       }
    }
    double time_spent_parallel = omp_get_wtime() - start;
    printf("Time: %lf\n", time_spent_parallel);

#else

    clock_t begin = clock();
    best1 = canonical_form(g1);
    best2 = canonical_form(g2);
    clock_t end = clock();
    printf("time: %f\n", ((double)(end - begin) / CLOCKS_PER_SEC));

#endif

    b = is_equal_orderings(g1->graph_m, g2->graph_m, best1, best2);

    free(best1);
    free(best2);

    return b;
}

/******************************************************************************
 * Helper functions
 *****************************************************************************/

Search_Node* create_search_node(int v)
{
    Search_Node *s;

    s = malloc(sizeof (Search_Node));
    s->u = 0;
    s->p = NULL;
    s->grp = NULL;
    s->depth = 0;
    s->nfixed = 0;
    s->on_best_path = false;
    s->next = NULL;

    s->cell_orbits = malloc(v * sizeof (int));
    for(int i = 0; i < v; ++i) s->cell_orbits[i] = -1;

    return s;
}

void free_search_tree(Search_Node *node)
{
    Search_Node *tmp;
    free_grps(node->grp, node->p->n);
    free_partition(node->p);
    while(node){
        tmp = node;
        node = node->next;
        free(tmp->cell_orbits);
        free(tmp);
    }
}

void print_search_node(Search_Node *s)
{
    printf("u: %d\n", s->u);
    printf("nfixed: %ld\n", s->nfixed);
    print_partition(s->p);
}

void print_perm(Permutation *perm, int size)
{
    Permutation *p = malloc(size * sizeof (int));
    for(int i = 0; i < size; ++i)
        p[i] = perm[i];

    int tmp, i, j;
    i = j = 0;
    printf("(");
    while(i < size){
        if(p[i] != -1){
            printf(" %d ", i);
            tmp = p[i];
            p[i] = -1;
            i = tmp;
        }
        else {
            ++j;
            i = j;
            if(i == size || p[i] != -1) printf(")");
            if(i < size && p[i] != -1) printf("(");
        }
     }
    printf("\n");
    free(p);
}










































