#include "exre.h"
#include <math.h>
#include <string.h>

/* First 36 elements: H..Kr. Beyond that we fall back to derived rules. */
static const char* SYM[] = {
    "?",
    "H","He",
    "Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr"
};

const char* element_symbol(int Z){
    if (Z < 0) return "?";
    if (Z < (int)(sizeof(SYM)/sizeof(SYM[0]))) return SYM[Z];
    return "X";   /* unknown nucleus: treated by derived rules */
}

/* Period/column derivation for unknown atoms */
static void shell_position(int Z, int* period, int* group){
    /* shell capacities: 2,8,8,18,18,32,32 */
    static const int caps[] = {2,8,8,18,18,32,32};
    int z = Z, p = 0;
    while (p < 7 && z > caps[p]) { z -= caps[p]; p++; }
    *period = p + 1;
    *group  = z;          /* 1..cap-of-shell */
}

int atom_period(int Z){
    int p,g; shell_position(Z,&p,&g);
    return p;
}

int is_noble_gas(int Z){
    return Z==2 || Z==10 || Z==18 || Z==36 || Z==54 || Z==86;
}

int max_valence(int Z){
    /* hypervalent allowance: period >= 3 can expand octet */
    int v = typical_valence(Z);
    int p,g; shell_position(Z,&p,&g);
    if (p >= 3) {
        /* S up to 6, P up to 5, halogens up to 7 */
        if (Z == 16) return 6;
        if (Z == 15) return 5;
        if (Z == 17 || Z == 35) return 7;
        return v + 2;
    }
    return v;
}

int typical_valence(int Z){
    if (Z <= 0) return 0;
    int p,g; shell_position(Z, &p, &g);
    int cap = (p==1)?2:(p<=3?8:(p<=5?18:32));
    /* atoms favour 0 or full shell. Distance to nearest closure: */
    int give = g;
    int take = cap - g;
    int v = give < take ? give : take;
    /* d-block bookkeeping: clamp to 8 for chemistry purposes */
    if (v > 4 && cap > 8) v = (v % 8 == 0) ? 4 : (v > 4 ? 8 - (v % 8) : v);
    if (v > 8) v = 8;
    return v;
}

double electronegativity(int Z){
    /* Pauling-ish table for common elements; derived fallback otherwise. */
    static const double EN[] = {
        0.0,
        2.20, 0.0,           /* H, He */
        0.98,1.57,2.04,2.55,3.04,3.44,3.98,0.0,    /* Li..Ne */
        0.93,1.31,1.61,1.90,2.19,2.58,3.16,0.0,    /* Na..Ar */
        0.82,1.00,1.36,1.54,1.63,1.66,1.55,1.83,1.88,1.91,1.90,1.65,
        1.81,2.01,2.18,2.55,2.96,3.00
    };
    if (Z>0 && Z < (int)(sizeof(EN)/sizeof(EN[0])) && EN[Z] > 0.0) return EN[Z];
    /* derived: scales with Z_eff / period */
    int p,g; shell_position(Z,&p,&g);
    double v = 1.0 + (double)g / 8.0 * 3.0;
    v -= 0.15 * (p - 2);
    if (v < 0.5) v = 0.5;
    if (v > 4.0) v = 4.0;
    return v;
}

/* Nuclear stability: liquid-drop-ish band + magic numbers. */
int nucleus_stable(int Z, int N){
    if (Z <= 0 || N < 0) return 0;
    /* valley of stability: N/Z ~ 1 for light, ~1.5 for heavy */
    double ideal = (Z <= 20) ? Z : Z * (1.0 + (Z-20)/180.0);
    double dev = fabs((double)N - ideal) / (ideal + 1.0);
    if (dev > 0.25) return 0;
    if (Z > 92 && N > 146) return 0;          /* beyond known stability */
    return 1;
}

/* Bond energies (kJ/mol). Negative because bonding lowers energy. */
static double table_lookup(int Z1, int Z2, int order){
    if (Z1 > Z2){ int t=Z1; Z1=Z2; Z2=t; }
    struct E { int a,b,o; double kj; };
    static const struct E T[] = {
        {1,1,1,436},  {1,6,1,413},  {1,7,1,391},  {1,8,1,463},
        {1,9,1,565},  {1,15,1,322}, {1,16,1,347}, {1,17,1,431},
        {6,6,1,347},  {6,6,2,614},  {6,6,3,839},
        {6,7,1,305},  {6,7,2,615},  {6,7,3,891},
        {6,8,1,358},  {6,8,2,799},
        {6,9,1,485},  {6,17,1,328}, {6,16,1,272},
        {7,7,1,160},  {7,7,2,418},  {7,7,3,945},
        {7,8,1,201},  {7,8,2,607},
        {8,8,1,146},  {8,8,2,498},
        {8,1,1,463},
        {16,16,1,266},{16,8,2,523}, {16,1,1,347},
        {9,9,1,159},  {17,17,1,243},{35,35,1,193},
        {15,15,1,201},{15,8,2,544}, {15,1,1,322},
        {0,0,0,0}
    };
    for (int k=0; T[k].a; k++){
        if (T[k].a==Z1 && T[k].b==Z2 && T[k].o==order) return T[k].kj;
    }
    return 0.0;
}

double bond_energy(int Z1, int Z2, int order){
    double kj = table_lookup(Z1, Z2, order);
    if (kj == 0.0){
        /* fallback: ~250 kJ/mol single * order, modulated by EN compatibility */
        double en1 = electronegativity(Z1);
        double en2 = electronegativity(Z2);
        double base = 240.0 * order;
        double polar = 96.0 * fabs(en1 - en2);   /* polar covalent stabilises */
        kj = base + polar;
    }
    return -kj;  /* stabilising */
}

double atom_self_energy(const Atom* a){
    double e = 0.0;
    if (!nucleus_stable(a->Z, a->N)) e += 2000.0;          /* unstable nucleus */
    int q = a->Z - a->e;
    e += 300.0 * q * q;                                    /* charged atoms cost */
    return e;
}
