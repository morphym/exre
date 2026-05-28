#include "exre.h"
#include <string.h>

/* Common Z pool, weighted toward organic chemistry. */
static const int Z_POOL[] = {
    1,1,1,1,1, 6,6,6,6, 7,7,7, 8,8,8, 9, 15, 16, 17, 11, 12, 19, 20, 35
};
static const int Z_POOL_N = (int)(sizeof(Z_POOL)/sizeof(Z_POOL[0]));

static int default_N(int Z){
    /* simple stable N guess */
    if (Z <= 1) return 0;
    if (Z <= 20) return Z;
    return (int)(Z * 1.3);
}

static void op_add_atom(Molecule* m, unsigned int* rng){
    if (m->n_atoms == 0 || m->n_atoms >= MAX_ATOMS) return;
    int Z = Z_POOL[xorshift32(rng) % Z_POOL_N];
    int idx = mol_add_atom(m, Z, default_N(Z), Z);
    if (idx < 0) return;
    int anchor = xorshift32(rng) % (m->n_atoms - 1);
    mol_add_bond(m, anchor, idx, 1);
}

static void op_remove_atom(Molecule* m, unsigned int* rng){
    if (m->n_atoms <= 1) return;
    int k = xorshift32(rng) % m->n_atoms;
    /* drop bonds touching k, then shift atoms */
    int w = 0;
    for (int i=0;i<m->n_bonds;i++){
        Bond b = m->bonds[i];
        if (b.i == k || b.j == k) continue;
        if (b.i > k) b.i--;
        if (b.j > k) b.j--;
        m->bonds[w++] = b;
    }
    m->n_bonds = w;
    for (int i=k;i<m->n_atoms-1;i++) m->atoms[i] = m->atoms[i+1];
    m->n_atoms--;
}

static void op_swap_atom(Molecule* m, unsigned int* rng){
    if (m->n_atoms == 0) return;
    int k = xorshift32(rng) % m->n_atoms;
    int Z = Z_POOL[xorshift32(rng) % Z_POOL_N];
    m->atoms[k].Z = Z;
    m->atoms[k].N = default_N(Z);
    m->atoms[k].e = Z;
}

static void op_add_bond(Molecule* m, unsigned int* rng){
    if (m->n_atoms < 2 || m->n_bonds >= MAX_BONDS) return;
    int i = xorshift32(rng) % m->n_atoms;
    int j = xorshift32(rng) % m->n_atoms;
    if (i == j) return;
    int order = 1 + (xorshift32(rng) % 3);
    mol_add_bond(m, i, j, order);
}

static void op_remove_bond(Molecule* m, unsigned int* rng){
    if (m->n_bonds == 0) return;
    int k = xorshift32(rng) % m->n_bonds;
    for (int i=k;i<m->n_bonds-1;i++) m->bonds[i] = m->bonds[i+1];
    m->n_bonds--;
}

static void op_bump_order(Molecule* m, unsigned int* rng){
    if (m->n_bonds == 0) return;
    int k = xorshift32(rng) % m->n_bonds;
    int delta = (xorshift32(rng) & 1) ? 1 : -1;
    int o = m->bonds[k].order + delta;
    if (o >= 1 && o <= 3) m->bonds[k].order = o;
}

static void op_mutate_nucleus(Molecule* m, unsigned int* rng){
    /* rare: tweak N or e to explore unknown isotopes / ions */
    if (m->n_atoms == 0) return;
    int k = xorshift32(rng) % m->n_atoms;
    int which = xorshift32(rng) % 2;
    int delta = (xorshift32(rng) & 1) ? 1 : -1;
    if (which == 0) m->atoms[k].N += delta;
    else            m->atoms[k].e += delta;
    if (m->atoms[k].N < 0) m->atoms[k].N = 0;
    if (m->atoms[k].e < 0) m->atoms[k].e = 0;
}

typedef void (*op_fn)(Molecule*, unsigned int*);
static op_fn OPS[] = {
    op_add_atom, op_remove_atom, op_swap_atom,
    op_add_bond, op_remove_bond, op_bump_order,
    op_mutate_nucleus
};
static const int N_OPS = (int)(sizeof(OPS)/sizeof(OPS[0]));

void mol_mutate(Molecule* m, unsigned int* rng){
    OPS[xorshift32(rng) % N_OPS](m, rng);
}
