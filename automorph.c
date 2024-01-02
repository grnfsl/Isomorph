/* Date for the last update: 14/01/2021
 * Author: Goran Faisal; grnhry@gmail.com
 *
 * This file contains an implementation of the Schreier–Sims algorithm, 
 * which is used in the stabilise(...) function in the isomorph.c file.
 */

#include <automorph.h>
#include <time.h>
void print_perm(Permutation *p, int size);

/******************************************************************************
 * Tests the membership of a generator in a group.
 *****************************************************************************/
bool is_grp_element(Permutation *y, Group *grp, int n)
{
    Permutation *p;
    int i, v;
    bool b;

    int tmp;
    tmp = y[n-1];
    y[n-1] = -1;

    for (i = 0; y[i] == i; ++i);
    y[n-1] = tmp;
    if(i == n-1 && y[i] == i) return true;

    v = y[grp->u];
    if(grp->coset[v] == NULL) return false;
    if(grp->coset_inv[v] == NULL) grp->coset_inv[v] = inverse_perm(grp->coset[v], n);
    p = multiply_perms(y, grp->coset_inv[v], n);

    b = is_grp_element(p, grp->st_grp_u, n);
    free(p);
    return b;
}

Permutation* multiply_perms(Permutation *a, Permutation *b, int n)
{
    Permutation *p;
    int i;

    p = malloc(n * sizeof (int));
    for(i = 0; i < n; ++i)
        p[i] = b[a[i]];
    return p;
}

Permutation* inverse_perm(Permutation *y, int n)
{
    Permutation *p;
    int i;

    p = malloc(n * sizeof (int));
    for(i = 0; i < n; ++i)
        p[y[i]] = i;
    return p;
}

Permutation* gen_multiply(Permutation *a, Permutation *b, Permutation *y, int n)
{
    Permutation *p;
    int i;

    p = malloc(n * sizeof (int));
    for(i = 0; i < n; ++i)
        p[i] = y[b[a[i]]];
    return p;
}

/******************************************************************************
 * Adds generator to the tower of stabilzer subgroups.
 *****************************************************************************/
void add_gen(Permutation *y, Search_Node *s)
{
    Permutation *p;
    Generator *gen;
    int i, m, w, v;
    size_t n;

    n = s->p->n;

    if(s->grp->st_grp_u == NULL){
        s->grp->st_grp_u = identity_grp(n);
        s->next->grp = s->grp->st_grp_u;

        if(s->u != INACTIVE)
            s->grp->u = s->u;
        else{
            s->grp->u = 0;
            while (y[s->grp->u] == s->grp->u) ++s->grp->u;
        }
        s->grp->coset[s->grp->u] = identity_perm(n);
        s->grp->orbit[s->grp->orbit_npts] = s->grp->u;
        ++s->grp->orbit_npts;
    }

    gen = create_gen(y);
    gen->next = s->grp->gens;
    s->grp->gens = gen;

    update_orbits(y, s->p->head, s->cell_orbits);

    m = s->grp->orbit_npts;
    i = 0;
    while (i < m) {
        v = s->grp->orbit[i];
        w = y[v];
        if(s->grp->coset[w] == NULL){
            s->grp->orbit[s->grp->orbit_npts] = w;
            ++s->grp->orbit_npts;
            s->grp->coset[w] = multiply_perms(s->grp->coset[v], y, n);
        }
        else{
            if(s->grp->coset_inv[w] == NULL)
                s->grp->coset_inv[w] = inverse_perm(s->grp->coset[w], n);
            p = gen_multiply(s->grp->coset[v], y, s->grp->coset_inv[w], n);
            if(!is_grp_element(p, s->grp->st_grp_u, n))
                add_gen(p, s->next);
            else free(p);
        }
        ++i;
    }

    while (i < s->grp->orbit_npts) {
        v = s->grp->orbit[i];
        gen = s->grp->gens;
        while (gen) {
            w = gen->perm[v];
            if(s->grp->coset[w] == NULL){
                s->grp->orbit[s->grp->orbit_npts] = w;
                ++s->grp->orbit_npts;
                s->grp->coset[w] = multiply_perms(s->grp->coset[v], gen->perm, n);
            }
            else{
                if(s->grp->coset_inv[w] == NULL)
                    s->grp->coset_inv[w] = inverse_perm(s->grp->coset[w], n);
                p = gen_multiply(s->grp->coset[v], gen->perm, s->grp->coset_inv[w], n);
                if(!is_grp_element(p, s->grp->st_grp_u, n))
                    add_gen(p, s->next);
                else free(p);
            }
            gen = gen->next;
        }
        ++i;
    }
}

void update_orbits(const Permutation *y, const Cell *c, int *cell_orbits)
{
    int *i, u, v, urep, vrep;
    for(i = c->first_p; i <= c->last_p; ++i){
        u = *i;
        v = y[u];
        urep = orbit_rep(u, cell_orbits);
        vrep = orbit_rep(v, cell_orbits);
        if(urep != vrep) merge_orbits(urep, vrep, cell_orbits);
    }
}

void merge_orbits(int urep, int vrep, int *cell_orbits)
{
    int usize, vsize, w;

    usize = -(cell_orbits[urep]);
    vsize = -(cell_orbits[vrep]);

    if(usize < vsize){
        cell_orbits[urep] = vrep;
        w = vrep;
    }
    else{
        cell_orbits[vrep] = urep;
        w = urep;
    }
    cell_orbits[w] = -(usize + vsize);
}

int orbit_rep(int v, int *cell_orbits)
{
    int w;

    if(cell_orbits[v] < 0) return v;
    w = orbit_rep(cell_orbits[v], cell_orbits);
    cell_orbits[v] = w;  //path compression
    return w;
}

//----------------------------------------Helper Functions-------------------------------

Generator* create_gen(Permutation *y)
{
    Generator *gen = malloc(sizeof (Generator));
    gen->perm = y;
    gen->next = NULL;
    return gen;
}

void free_gen(Generator *gen){
    Generator *tmp;

    while(gen){
        tmp = gen;
        gen = gen->next;
        free(tmp->perm);
        free(tmp);
    }
}

Coset** create_cosetreps(int n)
{
    Coset **c;
    int i;

    c = malloc(n * sizeof (Coset*));

    for(i = 0; i < n; ++i)
        c[i] = NULL;

    return c;
}

void free_cosetreps(Coset **c, int n)
{
    int i;
    for(i = 0; i < n; ++i)
        free(c[i]);
    free(c);
}

Group* identity_grp(int n)
{
    Group *grp;

    grp = malloc(sizeof (Group));

    grp->u = -1;
    grp->gens = NULL;

    grp->orbit = malloc(n * sizeof (int));
    grp->orbit_npts = 0;

    grp->coset = create_cosetreps(n);
    grp->coset_inv = create_cosetreps(n);

    grp->st_grp_u = NULL;

    return grp;
}

Permutation* identity_perm(int n)
{
    Permutation *perm = malloc(n * sizeof (int));

    for(int i = 0; i < n; ++i)
        perm[i] = i;
    return perm;
}

void free_grps(Group *grp, int n)
{
    Group *tmp;

    while(grp){
        free_gen(grp->gens);
        free_cosetreps(grp->coset, n);
        free_cosetreps(grp->coset_inv, n);
        free(grp->orbit);

        tmp = grp;
        grp = grp->st_grp_u;
        free(tmp);
    }
}

















































