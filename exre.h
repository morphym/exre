#ifndef EXRE_H
#define EXRE_H

#include <stdint.h>
#include <stddef.h>

#define MAX_ATOMS   64
#define MAX_BONDS  128
#define MAX_NAME    32

typedef struct {
    int Z;          /* protons    */
    int N;          /* neutrons   */
    int e;          /* electrons (total in atom; charge = Z - e) */
} Atom;

typedef struct {
    int i, j;       /* atom indices */
    int order;      /* 1,2,3 */
} Bond;

typedef struct {
    char  name[MAX_NAME];
    Atom  atoms[MAX_ATOMS];
    int   n_atoms;
    Bond  bonds[MAX_BONDS];
    int   n_bonds;
    double energy;  /* cached score (lower = more stable) */
} Molecule;

/* chem.c */
const char* element_symbol(int Z);
int   typical_valence(int Z);
int   max_valence(int Z);                        /* hypervalence cap */
int   is_noble_gas(int Z);
int   atom_period(int Z);
double electronegativity(int Z);
int   nucleus_stable(int Z, int N);
double bond_energy(int Z1, int Z2, int order);   /* kJ/mol, negative = stabilising */
double atom_self_energy(const Atom* a);

/* mol.c */
void   mol_clear(Molecule* m);
int    mol_add_atom(Molecule* m, int Z, int N, int e);
int    mol_add_bond(Molecule* m, int i, int j, int order);
int    mol_valence_of(const Molecule* m, int idx);
double mol_score(Molecule* m);                  /* recompute + cache */
int    mol_is_connected(const Molecule* m);
void   mol_print(const Molecule* m);
void   mol_formula(const Molecule* m, char* out, size_t cap);

/* mutate.c */
void   mol_mutate(Molecule* m, unsigned int* rng);

/* seeds.c */
int    load_seeds(Molecule* arr, int cap);

/* rng */
static inline unsigned int xorshift32(unsigned int* s){
    unsigned int x = *s; x ^= x<<13; x ^= x>>17; x ^= x<<5; *s = x; return x;
}
static inline double urand(unsigned int* s){
    return (xorshift32(s) & 0xFFFFFF) / (double)0x1000000;
}

#endif
