#include "exre.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

void mol_clear(Molecule* m){
    m->n_atoms = 0; m->n_bonds = 0; m->energy = 0.0; m->name[0] = 0;
}

int mol_add_atom(Molecule* m, int Z, int N, int e){
    if (m->n_atoms >= MAX_ATOMS) return -1;
    m->atoms[m->n_atoms] = (Atom){Z,N,e};
    return m->n_atoms++;
}

int mol_add_bond(Molecule* m, int i, int j, int order){
    if (m->n_bonds >= MAX_BONDS) return -1;
    if (i == j || i < 0 || j < 0 || i >= m->n_atoms || j >= m->n_atoms) return -1;
    for (int k=0;k<m->n_bonds;k++){
        Bond* b = &m->bonds[k];
        if ((b->i==i && b->j==j) || (b->i==j && b->j==i)) return -1;
    }
    m->bonds[m->n_bonds++] = (Bond){i,j,order};
    return 0;
}

int mol_valence_of(const Molecule* m, int idx){
    int v = 0;
    for (int k=0;k<m->n_bonds;k++){
        const Bond* b = &m->bonds[k];
        if (b->i==idx || b->j==idx) v += b->order;
    }
    return v;
}

/* connectivity via DFS */
static void dfs(const Molecule* m, int u, int* seen){
    seen[u] = 1;
    for (int k=0;k<m->n_bonds;k++){
        const Bond* b = &m->bonds[k];
        int v = -1;
        if (b->i==u) v = b->j;
        else if (b->j==u) v = b->i;
        if (v >= 0 && !seen[v]) dfs(m, v, seen);
    }
}
int mol_is_connected(const Molecule* m){
    if (m->n_atoms <= 1) return 1;
    int seen[MAX_ATOMS] = {0};
    dfs(m, 0, seen);
    for (int i=0;i<m->n_atoms;i++) if (!seen[i]) return 0;
    return 1;
}

/* shortest cycle through bond k (BFS i->j without using bond k). Returns 0 if no ring. */
static int shortest_ring_through(const Molecule* m, int kbond){
    int src = m->bonds[kbond].i, dst = m->bonds[kbond].j;
    int dist[MAX_ATOMS];
    for (int i=0;i<m->n_atoms;i++) dist[i] = -1;
    int queue[MAX_ATOMS], qh=0, qt=0;
    queue[qt++] = src; dist[src] = 0;
    while (qh < qt){
        int u = queue[qh++];
        for (int b=0;b<m->n_bonds;b++){
            if (b == kbond) continue;
            const Bond* bb = &m->bonds[b];
            int v = -1;
            if (bb->i == u) v = bb->j;
            else if (bb->j == u) v = bb->i;
            if (v < 0 || dist[v] >= 0) continue;
            dist[v] = dist[u] + 1;
            if (v == dst) return dist[v] + 1;
            queue[qt++] = v;
        }
    }
    return 0;
}

double mol_score(Molecule* m){
    double E = 0.0;

    /* (1) atom intrinsic energy (nuclear + ionisation) */
    for (int i=0;i<m->n_atoms;i++) E += atom_self_energy(&m->atoms[i]);

    /* (2) bond energies + ionic-lattice bonus when ΔEN > 1.7, saturated per atom */
    int ionic_used[MAX_ATOMS] = {0};
    for (int k=0;k<m->n_bonds;k++){
        const Bond* b = &m->bonds[k];
        int Z1 = m->atoms[b->i].Z, Z2 = m->atoms[b->j].Z;
        E += bond_energy(Z1, Z2, b->order);
        double den = fabs(electronegativity(Z1) - electronegativity(Z2));
        if (den > 1.7){
            int cap1 = typical_valence(Z1);
            int cap2 = typical_valence(Z2);
            if (ionic_used[b->i] < cap1 && ionic_used[b->j] < cap2){
                E += -400.0 * (den - 1.7);    /* one ionic credit per atom-valence-slot */
                ionic_used[b->i]++;
                ionic_used[b->j]++;
            }
        }
        /* (6) noble gas isolation */
        if (is_noble_gas(Z1) || is_noble_gas(Z2)) E += 1500.0;
    }

    /* (3) valence law: soft quadratic + hard hypervalence cap */
    int net_charge = 0;
    for (int i=0;i<m->n_atoms;i++){
        int Z = m->atoms[i].Z;
        int want = typical_valence(Z);
        int have = mol_valence_of(m, i);
        int hardcap = max_valence(Z);
        int d = have - want;
        E += 180.0 * d * d;
        if (have > hardcap){
            int over = have - hardcap;
            E += 700.0 * over * over;             /* expanded octet violation */
        }
        net_charge += (Z - m->atoms[i].e);
    }

    /* (4) charge balance — neutral species are favoured */
    E += 400.0 * net_charge * net_charge;

    /* (5) ring strain via cyclomatic search */
    int cyc = m->n_bonds - m->n_atoms + 1;
    if (cyc > 0){
        /* find smallest ring through each ring-bond and apply strain */
        int counted = 0;
        for (int k=0;k<m->n_bonds && counted < cyc;k++){
            int r = shortest_ring_through(m, k);
            if (r >= 3){
                if      (r == 3) E += 250.0;
                else if (r == 4) E += 120.0;
                else if (r == 5) E += 25.0;
                else if (r == 6) E -= 60.0;       /* small aromatic-style bonus */
                else if (r == 7) E += 20.0;
                counted++;
            }
        }
    }

    /* (7) size/entropy: linear + quadratic above ~8 atoms */
    E += 90.0 * m->n_atoms;
    if (m->n_atoms > 8){
        int over = m->n_atoms - 8;
        E += 25.0 * over * over;
    }

    /* connectivity */
    if (!mol_is_connected(m)) E += 1500.0;

    m->energy = E;
    return E;
}

void mol_formula(const Molecule* m, char* out, size_t cap){
    int counts[120] = {0};
    int has[120]    = {0};
    int unknown = 0;
    for (int i=0;i<m->n_atoms;i++){
        int Z = m->atoms[i].Z;
        if (Z >= 0 && Z < 120){ counts[Z]++; has[Z]=1; }
        else unknown++;
    }
    size_t n = 0;
    /* Hill order: C, H, then alphabetical */
    int order[120], nord = 0;
    if (has[6]) order[nord++] = 6;
    if (has[1]) order[nord++] = 1;
    for (int Z=2; Z<120; Z++){
        if (Z==6) continue;
        if (has[Z]) order[nord++] = Z;
    }
    for (int k=0;k<nord;k++){
        int Z = order[k];
        int c = snprintf(out+n, cap-n, "%s%s",
                         element_symbol(Z),
                         counts[Z]>1 ? "" : "");
        if (c < 0 || (size_t)c >= cap-n) return;
        n += c;
        if (counts[Z] > 1){
            c = snprintf(out+n, cap-n, "%d", counts[Z]);
            if (c < 0 || (size_t)c >= cap-n) return;
            n += c;
        }
    }
    if (unknown){
        snprintf(out+n, cap-n, "X%d", unknown);
    }
}

void mol_print(const Molecule* m){
    char f[128] = {0};
    mol_formula(m, f, sizeof f);
    printf("  %-16s  atoms=%d bonds=%d  E=%.1f kJ/mol  [%s]\n",
           f, m->n_atoms, m->n_bonds, m->energy,
           m->name[0] ? m->name : "?");
}
