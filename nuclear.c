#include "nuclear.h"
#include <math.h>

/* Standard SEMF coefficients (Krane / Rohlf), MeV */
static const double aV = 15.8;
static const double aS = 17.8;
static const double aC = 0.711;
static const double aA = 23.7;
static const double aP = 11.18;

double semf_binding(int Z, int N){
    if (Z < 1 || N < 0) return 0.0;
    int A = Z + N;
    if (A < 1) return 0.0;
    double A23 = pow((double)A, 2.0/3.0);
    double A13 = pow((double)A, 1.0/3.0);
    double B = aV * A
             - aS * A23
             - aC * (double)Z * (Z - 1) / A13
             - aA * (double)(N - Z) * (N - Z) / A;
    double pairing = 0.0;
    int eZ = (Z % 2 == 0), eN = (N % 2 == 0);
    if (eZ && eN)         pairing =  aP / sqrt((double)A);
    else if (!eZ && !eN)  pairing = -aP / sqrt((double)A);
    return B + pairing;
}

static const int MAGIC[] = {2,8,20,28,50,82,126,184};
static const int NM = (int)(sizeof(MAGIC)/sizeof(MAGIC[0]));

double magic_bonus(int Z, int N){
    double b = 0;
    int magicZ = 0, magicN = 0;
    for (int k=0;k<NM;k++){
        if (Z == MAGIC[k]){ b += 1.5; magicZ = 1; }
        if (N == MAGIC[k]){ b += 1.5; magicN = 1; }
    }
    if (magicZ && magicN) b += 3.0;  /* doubly magic — extra stability */
    return b;
}

double binding_per_nucleon(int Z, int N){
    int A = Z + N;
    if (A <= 0) return 0;
    double B = semf_binding(Z, N) + magic_bonus(Z, N);
    return B / A;
}

double nucleus_energy(int Z, int N){
    return -(semf_binding(Z, N) + magic_bonus(Z, N));
}

int nucleus_bound(int Z, int N){
    if (Z < 1 || N < 0) return 0;
    int A = Z + N;
    if (A < 1) return 0;
    if (Z > 138) return 0;                       /* beyond predicted stability */
    double B = semf_binding(Z, N) + magic_bonus(Z, N);
    if (B <= 0) return 0;                        /* unbound */
    if (binding_per_nucleon(Z, N) < 1.0) return 0;
    /* approximate drip lines via N/Z ratio */
    double r = (double)N / (double)Z;
    if (Z <= 8){
        if (r < 0.4 || r > 2.5) return 0;
    } else if (Z <= 20){
        if (r < 0.8 || r > 1.8) return 0;
    } else if (Z <= 82){
        if (r < 1.0 || r > 1.8) return 0;
    } else {
        if (r < 1.3 || r > 1.7) return 0;
    }
    /* fission barrier: Z²/A > ~49 collapses */
    if ((double)Z * Z / A > 49.0) return 0;
    return 1;
}

int is_closed_shell_e(int e){
    static const int NE[] = {2,10,18,36,54,86,118};
    for (int k=0;k<7;k++) if (e == NE[k]) return 1;
    return 0;
}
