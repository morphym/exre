/* exre layer-1: stochastic element discovery via MCMC over (Z, N, e). */
#include "exre.h"
#include "nuclear.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAXZ 140
#define MAXN 260

static double best_BA[MAXZ][MAXN];
static int    seen[MAXZ][MAXN];

/* IUPAC systematic name for elements without a symbol in our table */
static const char* DIG[] = {
    "nil","un","bi","tri","quad","pent","hex","sept","oct","enn"
};
static void iupac_name(int Z, char* out, size_t cap){
    if (Z <= 0){ snprintf(out, cap, "neutronium"); return; }
    char buf[64] = {0};
    int digs[8], nd=0, z=Z;
    while (z > 0){ digs[nd++] = z % 10; z /= 10; }
    size_t n = 0;
    for (int i=nd-1; i>=0 && n<sizeof(buf)-1; i--){
        const char* d = DIG[digs[i]];
        while (*d && n<sizeof(buf)-1) buf[n++] = *d++;
    }
    buf[n] = 0;
    if (n > 0) buf[0] = (buf[0] >= 'a' && buf[0] <= 'z') ? buf[0] - 32 : buf[0];
    snprintf(out, cap, "%s-ium", buf);
}

static void element_name(int Z, char* out, size_t cap){
    const char* s = element_symbol(Z);
    if (!s || s[0] == 'X' || s[0] == '?') iupac_name(Z, out, cap);
    else snprintf(out, cap, "%s", s);
}

/* Total atom energy in MeV. */
static double total_energy(int Z, int N, int e){
    double E = nucleus_energy(Z, N);                  /* MeV */
    int q = Z - e;
    E += 0.005 * q * q;                                /* ~5 eV per charge^2 */
    if (is_closed_shell_e(e)) E -= 0.01;               /* 10 eV shell bonus */
    return E;
}

int main(int argc, char** argv){
    long iters     = (argc > 1) ? atol(argv[1]) : 500000;
    double T       = (argc > 2) ? atof(argv[2]) : 2.0;          /* MeV */
    unsigned int rng = (argc > 3) ? (unsigned)atol(argv[3]) : (unsigned)time(NULL);
    if (!rng) rng = 0xC0DE;

    FILE* LOG = fopen("exre_elements.log", "w");

    printf("exre layer-1: stochastic element discovery\n");
    printf("iters=%ld  T=%.2f MeV  rng=%u\n\n", iters, T, rng);
    if (LOG) fprintf(LOG, "# Z\tN\tA\te\tB/A_MeV\tname\n");

    /* multiple walkers spread across Z so we don't get stuck near iron */
    enum { NWALK = 16 };
    int  Z[NWALK], N[NWALK], e[NWALK];
    double Ec[NWALK];
    for (int w=0; w<NWALK; w++){
        Z[w] = 1 + w * 8;
        N[w] = Z[w];
        e[w] = Z[w];
        if (!nucleus_bound(Z[w], N[w])){ N[w] = Z[w] + 1; }
        Ec[w] = total_energy(Z[w], N[w], e[w]);
    }

    for (int z=0; z<MAXZ; z++) for (int n=0; n<MAXN; n++) best_BA[z][n] = -1;
    long discoveries = 0;

    for (long it = 0; it < iters; it++){
        int w = xorshift32(&rng) % NWALK;
        int Zc = Z[w], Nc = N[w], ec = e[w];
        int op = xorshift32(&rng) % 5;
        int d  = (xorshift32(&rng) & 1) ? +1 : -1;
        switch (op){
            case 0: Zc += d; ec += d; break;                 /* neutral Z step */
            case 1: Nc += d; break;
            case 2: ec += d; break;
            case 3: Zc += d; Nc += d; ec += d; break;        /* alpha-ish step */
            case 4:                                          /* random jump */
                Zc = 1 + (xorshift32(&rng) % 135);
                Nc = xorshift32(&rng) % 220;
                ec = Zc;
                break;
        }
        if (Zc < 1 || Zc >= MAXZ || Nc < 0 || Nc >= MAXN || ec < 0) continue;
        if (!nucleus_bound(Zc, Nc)) continue;

        double Enew = total_energy(Zc, Nc, ec);
        double dE = Enew - Ec[w];
        double p  = (dE <= 0) ? 1.0 : exp(-dE / T);
        if (urand(&rng) < p){
            Z[w] = Zc; N[w] = Nc; e[w] = ec; Ec[w] = Enew;
            double ba = binding_per_nucleon(Zc, Nc);
            if (!seen[Zc][Nc] || ba > best_BA[Zc][Nc]){
                seen[Zc][Nc] = 1;
                best_BA[Zc][Nc] = ba;
                discoveries++;
                if (LOG){
                    char nm[32]; element_name(Zc, nm, sizeof nm);
                    fprintf(LOG, "%d\t%d\t%d\t%d\t%.4f\t%s-%d\n",
                            Zc, Nc, Zc+Nc, ec, ba, nm, Zc+Nc);
                }
            }
        }
    }

    printf("discoveries (unique stable nuclides): %ld\n\n", discoveries);

    /* most stable nuclide per element Z */
    typedef struct { int Z, N; double ba; } Rec;
    Rec recs[MAXZ]; int nrec = 0;
    for (int z=1; z<MAXZ; z++){
        double best = -1; int bN = -1;
        for (int n=0; n<MAXN; n++){
            if (seen[z][n] && best_BA[z][n] > best){ best = best_BA[z][n]; bN = n; }
        }
        if (bN >= 0) recs[nrec++] = (Rec){z, bN, best};
    }
    /* sort by B/A descending */
    for (int i=0; i<nrec; i++)
      for (int j=i+1; j<nrec; j++)
        if (recs[j].ba > recs[i].ba){ Rec t=recs[i]; recs[i]=recs[j]; recs[j]=t; }

    printf("--- top-30 most stable nuclides (highest B/A) ---\n");
    int top = nrec < 30 ? nrec : 30;
    for (int i=0; i<top; i++){
        char nm[32]; element_name(recs[i].Z, nm, sizeof nm);
        int A = recs[i].Z + recs[i].N;
        int mz=0, mn=0;
        for (int k=0;k<8;k++){
            int M[] = {2,8,20,28,50,82,126,184};
            if (recs[i].Z == M[k]) mz=1;
            if (recs[i].N == M[k]) mn=1;
        }
        printf("%2d. Z=%-3d N=%-3d A=%-3d  B/A=%6.3f MeV  %-14s %s%s\n",
               i+1, recs[i].Z, recs[i].N, A, recs[i].ba, nm,
               mz?"[magic Z]":"", mn?"[magic N]":"");
    }

    printf("\n--- known elements rediscovered (Z = 1..36) ---\n");
    for (int z=1; z<=36; z++){
        for (int i=0;i<nrec;i++) if (recs[i].Z == z){
            char nm[16]; element_name(z, nm, sizeof nm);
            printf("  %-3s Z=%-3d  most-stable N=%-3d A=%-3d  B/A=%.3f MeV\n",
                   nm, z, recs[i].N, z + recs[i].N, recs[i].ba);
            break;
        }
    }

    printf("\n--- superheavy candidates (Z > 118): island-of-stability scan ---\n");
    int nsh = 0;
    for (int i=0; i<nrec; i++){
        if (recs[i].Z > 118){
            char nm[32]; element_name(recs[i].Z, nm, sizeof nm);
            int A = recs[i].Z + recs[i].N;
            printf("  Z=%-3d N=%-3d A=%-3d  B/A=%.3f MeV  %s-%d\n",
                   recs[i].Z, recs[i].N, A, recs[i].ba, nm, A);
            nsh++;
            if (nsh >= 15) break;
        }
    }
    if (!nsh) printf("  (none survived the bound + fission-barrier criteria)\n");

    if (LOG) fclose(LOG);
    return 0;
}
