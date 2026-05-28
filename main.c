#include "exre.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define POP_CAP   512
#define SEED_CAP   64

static Molecule pop[POP_CAP];
static int      pop_n;

static FILE* LOG;

static void emit(const Molecule* m, const char* tag){
    char f[128]={0};
    mol_formula(m, f, sizeof f);
    printf("[%s] %-16s  E=%.1f  atoms=%d bonds=%d  parent=%s\n",
           tag, f, m->energy, m->n_atoms, m->n_bonds,
           m->name[0]?m->name:"?");
    if (LOG){
        fprintf(LOG, "%s\t%s\tE=%.3f\tatoms=%d\tbonds=%d\tparent=%s\n",
                tag, f, m->energy, m->n_atoms, m->n_bonds,
                m->name[0]?m->name:"?");
        fflush(LOG);
    }
}

int main(int argc, char** argv){
    long iters     = (argc > 1) ? atol(argv[1]) : 200000;
    double T       = (argc > 2) ? atof(argv[2]) : 600.0;   /* "temperature" */
    unsigned int rng = (argc > 3) ? (unsigned)atol(argv[3]) : (unsigned)time(NULL);
    if (rng == 0) rng = 0xC0FFEE;

    LOG = fopen("exre.log", "w");
    if (!LOG){ perror("exre.log"); }

    pop_n = load_seeds(pop, SEED_CAP);
    printf("exre: stochastic chemical discovery\n");
    printf("seeds: %d   iters: %ld   T: %.1f   rng: %u\n\n", pop_n, iters, T, rng);
    if (LOG) fprintf(LOG, "# exre run T=%.3f iters=%ld rng=%u seeds=%d\n",
                     T, iters, rng, pop_n);

    /* baseline summary */
    double best_E = 1e18;
    int    best_i = 0;
    for (int i=0;i<pop_n;i++){
        emit(&pop[i], "seed");
        if (pop[i].energy < best_E){ best_E = pop[i].energy; best_i = i; }
    }
    printf("\nbest seed: %s (E=%.1f)\n\n", pop[best_i].name, best_E);

    long accepted = 0, rejected = 0, novel = 0;
    for (long it = 0; it < iters; it++){
        int src = xorshift32(&rng) % pop_n;
        Molecule cand = pop[src];
        cand.name[0] = 0;
        strncpy(cand.name, pop[src].name, MAX_NAME-1);

        /* apply 1-3 random mutations */
        int mut = 1 + (xorshift32(&rng) % 3);
        for (int k=0;k<mut;k++) mol_mutate(&cand, &rng);

        if (cand.n_atoms < 1) { rejected++; continue; }
        mol_score(&cand);

        double dE = cand.energy - pop[src].energy;
        double p  = (dE <= 0) ? 1.0 : exp(-dE / T);
        if (urand(&rng) < p){
            accepted++;
            /* novelty: cheaper than canonical form — check if E improves global best
               or if formula is new among current population. */
            char fnew[128]={0}; mol_formula(&cand, fnew, sizeof fnew);
            int seen = 0;
            for (int i=0;i<pop_n;i++){
                char fi[128]={0}; mol_formula(&pop[i], fi, sizeof fi);
                if (strcmp(fi, fnew)==0){ seen = 1; break; }
            }
            if (!seen && pop_n < POP_CAP){
                pop[pop_n] = cand;
                emit(&pop[pop_n], "novel");
                pop_n++;
                novel++;
            } else if (cand.energy < pop[src].energy - 1.0){
                pop[src] = cand;
            }
            if (cand.energy < best_E){
                best_E = cand.energy;
                emit(&cand, "best ");
            }
        } else {
            rejected++;
        }
    }

    printf("\n--- summary ---\n");
    printf("iters=%ld  accepted=%ld  rejected=%ld  novel=%ld  pop=%d\n",
           iters, accepted, rejected, novel, pop_n);
    printf("global best E=%.1f\n", best_E);
    if (LOG){
        fprintf(LOG, "# end iters=%ld accepted=%ld rejected=%ld novel=%ld pop=%d best=%.3f\n",
                iters, accepted, rejected, novel, pop_n, best_E);
        fclose(LOG);
    }
    return 0;
}
