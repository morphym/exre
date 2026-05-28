#ifndef NUCLEAR_H
#define NUCLEAR_H

/* Weizsäcker semi-empirical mass formula + shell corrections.
 * All energies in MeV. */
double semf_binding(int Z, int N);          /* total binding energy */
double magic_bonus(int Z, int N);            /* shell-closure bonus */
double binding_per_nucleon(int Z, int N);    /* (B + magic) / A */
double nucleus_energy(int Z, int N);         /* -B - magic, MeV (lower = more stable) */
int    nucleus_bound(int Z, int N);          /* drip-line + bound check */

/* Electronic shell closure (noble gas) */
int    is_closed_shell_e(int e);

#endif
