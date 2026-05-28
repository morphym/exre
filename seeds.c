#include "exre.h"
#include <string.h>

/* Helpers for compact seed construction. */
static int add(Molecule* m, int Z){
    int N = (Z <= 1) ? 0 : Z;
    return mol_add_atom(m, Z, N, Z);
}
static void bond(Molecule* m, int i, int j, int o){ mol_add_bond(m,i,j,o); }

/* Each builder returns a populated molecule. */
static void seed_H2(Molecule* m){
    mol_clear(m); strcpy(m->name,"H2");
    int a=add(m,1), b=add(m,1); bond(m,a,b,1);
}
static void seed_O2(Molecule* m){
    mol_clear(m); strcpy(m->name,"O2");
    int a=add(m,8), b=add(m,8); bond(m,a,b,2);
}
static void seed_N2(Molecule* m){
    mol_clear(m); strcpy(m->name,"N2");
    int a=add(m,7), b=add(m,7); bond(m,a,b,3);
}
static void seed_H2O(Molecule* m){
    mol_clear(m); strcpy(m->name,"water");
    int o=add(m,8), h1=add(m,1), h2=add(m,1);
    bond(m,o,h1,1); bond(m,o,h2,1);
}
static void seed_CO2(Molecule* m){
    mol_clear(m); strcpy(m->name,"CO2");
    int c=add(m,6), o1=add(m,8), o2=add(m,8);
    bond(m,c,o1,2); bond(m,c,o2,2);
}
static void seed_CH4(Molecule* m){
    mol_clear(m); strcpy(m->name,"methane");
    int c=add(m,6);
    for (int k=0;k<4;k++){ int h=add(m,1); bond(m,c,h,1); }
}
static void seed_NH3(Molecule* m){
    mol_clear(m); strcpy(m->name,"ammonia");
    int n=add(m,7);
    for (int k=0;k<3;k++){ int h=add(m,1); bond(m,n,h,1); }
}
static void seed_HCN(Molecule* m){
    mol_clear(m); strcpy(m->name,"HCN");
    int h=add(m,1), c=add(m,6), n=add(m,7);
    bond(m,h,c,1); bond(m,c,n,3);
}
static void seed_CO(Molecule* m){
    mol_clear(m); strcpy(m->name,"CO");
    int c=add(m,6), o=add(m,8); bond(m,c,o,3);
}
static void seed_NO(Molecule* m){
    mol_clear(m); strcpy(m->name,"NO");
    int n=add(m,7), o=add(m,8); bond(m,n,o,2);
}
static void seed_H2S(Molecule* m){
    mol_clear(m); strcpy(m->name,"H2S");
    int s=add(m,16), h1=add(m,1), h2=add(m,1);
    bond(m,s,h1,1); bond(m,s,h2,1);
}
static void seed_HCl(Molecule* m){
    mol_clear(m); strcpy(m->name,"HCl");
    int h=add(m,1), c=add(m,17); bond(m,h,c,1);
}
static void seed_HF(Molecule* m){
    mol_clear(m); strcpy(m->name,"HF");
    int h=add(m,1), f=add(m,9); bond(m,h,f,1);
}
static void seed_C2H6(Molecule* m){
    mol_clear(m); strcpy(m->name,"ethane");
    int c1=add(m,6), c2=add(m,6); bond(m,c1,c2,1);
    for (int k=0;k<3;k++){ int h=add(m,1); bond(m,c1,h,1); }
    for (int k=0;k<3;k++){ int h=add(m,1); bond(m,c2,h,1); }
}
static void seed_C2H4(Molecule* m){
    mol_clear(m); strcpy(m->name,"ethene");
    int c1=add(m,6), c2=add(m,6); bond(m,c1,c2,2);
    for (int k=0;k<2;k++){ int h=add(m,1); bond(m,c1,h,1); }
    for (int k=0;k<2;k++){ int h=add(m,1); bond(m,c2,h,1); }
}
static void seed_C2H2(Molecule* m){
    mol_clear(m); strcpy(m->name,"ethyne");
    int c1=add(m,6), c2=add(m,6); bond(m,c1,c2,3);
    int h1=add(m,1), h2=add(m,1);
    bond(m,c1,h1,1); bond(m,c2,h2,1);
}
static void seed_CH3OH(Molecule* m){
    mol_clear(m); strcpy(m->name,"methanol");
    int c=add(m,6), o=add(m,8), h1=add(m,1);
    bond(m,c,o,1); bond(m,o,h1,1);
    for (int k=0;k<3;k++){ int h=add(m,1); bond(m,c,h,1); }
}
static void seed_CH3NH2(Molecule* m){
    mol_clear(m); strcpy(m->name,"methylamine");
    int c=add(m,6), n=add(m,7);
    bond(m,c,n,1);
    for (int k=0;k<3;k++){ int h=add(m,1); bond(m,c,h,1); }
    for (int k=0;k<2;k++){ int h=add(m,1); bond(m,n,h,1); }
}
static void seed_HCHO(Molecule* m){
    mol_clear(m); strcpy(m->name,"formaldehyde");
    int c=add(m,6), o=add(m,8), h1=add(m,1), h2=add(m,1);
    bond(m,c,o,2); bond(m,c,h1,1); bond(m,c,h2,1);
}
static void seed_HCOOH(Molecule* m){
    mol_clear(m); strcpy(m->name,"formic_acid");
    int c=add(m,6), o1=add(m,8), o2=add(m,8), h1=add(m,1), h2=add(m,1);
    bond(m,c,o1,2); bond(m,c,o2,1); bond(m,o2,h1,1); bond(m,c,h2,1);
}
static void seed_C6H6(Molecule* m){
    mol_clear(m); strcpy(m->name,"benzene_kekule");
    int c[6];
    for (int k=0;k<6;k++) c[k]=add(m,6);
    for (int k=0;k<6;k++){
        int order = (k%2==0) ? 1 : 2;
        bond(m, c[k], c[(k+1)%6], order);
    }
    for (int k=0;k<6;k++){ int h=add(m,1); bond(m,c[k],h,1); }
}
static void seed_SO2(Molecule* m){
    mol_clear(m); strcpy(m->name,"SO2");
    int s=add(m,16), o1=add(m,8), o2=add(m,8);
    bond(m,s,o1,2); bond(m,s,o2,2);
}
static void seed_PH3(Molecule* m){
    mol_clear(m); strcpy(m->name,"phosphine");
    int p=add(m,15);
    for (int k=0;k<3;k++){ int h=add(m,1); bond(m,p,h,1); }
}
static void seed_SiH4(Molecule* m){
    mol_clear(m); strcpy(m->name,"silane");
    int s=add(m,14);
    for (int k=0;k<4;k++){ int h=add(m,1); bond(m,s,h,1); }
}
static void seed_HOOH(Molecule* m){
    mol_clear(m); strcpy(m->name,"hydrogen_peroxide");
    int o1=add(m,8), o2=add(m,8), h1=add(m,1), h2=add(m,1);
    bond(m,o1,o2,1); bond(m,o1,h1,1); bond(m,o2,h2,1);
}
static void seed_N2H4(Molecule* m){
    mol_clear(m); strcpy(m->name,"hydrazine");
    int n1=add(m,7), n2=add(m,7);
    bond(m,n1,n2,1);
    for (int k=0;k<2;k++){ int h=add(m,1); bond(m,n1,h,1); }
    for (int k=0;k<2;k++){ int h=add(m,1); bond(m,n2,h,1); }
}
static void seed_CH3Cl(Molecule* m){
    mol_clear(m); strcpy(m->name,"chloromethane");
    int c=add(m,6), cl=add(m,17);
    bond(m,c,cl,1);
    for (int k=0;k<3;k++){ int h=add(m,1); bond(m,c,h,1); }
}

typedef void (*seedfn)(Molecule*);
static seedfn TABLE[] = {
    seed_H2, seed_O2, seed_N2, seed_H2O, seed_CO2, seed_CH4, seed_NH3,
    seed_HCN, seed_CO, seed_NO, seed_H2S, seed_HCl, seed_HF,
    seed_C2H6, seed_C2H4, seed_C2H2, seed_CH3OH, seed_CH3NH2,
    seed_HCHO, seed_HCOOH, seed_C6H6, seed_SO2, seed_PH3, seed_SiH4,
    seed_HOOH, seed_N2H4, seed_CH3Cl
};

int load_seeds(Molecule* arr, int cap){
    int n = (int)(sizeof(TABLE)/sizeof(TABLE[0]));
    if (n > cap) n = cap;
    for (int i=0;i<n;i++){
        TABLE[i](&arr[i]);
        mol_score(&arr[i]);
    }
    return n;
}
