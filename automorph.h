#ifndef AUTOMORPH_H
#define AUTOMORPH_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <utils.h>

#define INACTIVE -1 //to mark a node that it is inaactive when basis of group is longer than partition

typedef int Permutation, Coset;

#ifndef TYPES_GEN
#define TYPES_GEN

typedef struct Generator{
    Permutation *perm;
    struct Generator *next;
} Generator;

#endif // TYPES_GEN

#ifndef TYPES_GROUP
#define TYPES_GROUP

typedef struct Group{
    Generator *gens;
    int u;
    int orbit_npts;
    Permutation *orbit;
    Coset **coset;
    Coset **coset_inv;
    struct Group *st_grp_u;
} Group;

#endif // TYPES_GROUP

#ifndef TYPES_SN
#define TYPES_SN

typedef struct Search_Node{
    Partition *p;
    int u;
    Group *grp;
    int *cell_orbits;
    int depth;
    size_t nfixed;
    bool on_best_path;
    struct Search_Node *next;
} Search_Node;

#endif // TYPES_SN

bool is_grp_element(Permutation *y, Group *grp, int n);
Permutation* multiply_perms(Permutation *a, Permutation *b, int n);
Permutation* inverse_perm(Permutation *y, int n);
Permutation* gen_multiply(Permutation *a, Permutation *b, Permutation *y, int n);
void add_gen(Permutation *y, Search_Node *s);

void update_orbits(const Permutation *y, const Cell *c, int *cell_orbits);
void merge_orbits(int urep, int vrep, int *cell_orbits);
int orbit_rep(int v, int *cell_orbits);
//----------------------------------------Helper Functions-------------------------------

Generator* create_gen(Permutation *y);
void free_gen(Generator *gen);
Coset** create_cosetreps(int v);
void free_cosetreps(Coset **c, int n);
Group* identity_grp(int v); //create identity group
Permutation* identity_perm(int v);
void free_grps(Group *grp, int n);

#endif // AUTOMORPH_H
