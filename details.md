# exre — law reference

The physics encoded in the program. Every constant in the source can be traced back to a line here.

---

## Layer 1 — nuclear laws (file: `nuclear.c`)

The total atom energy is

```
E_atom(Z, N, e) = − B(Z, N) − magic_bonus(Z, N) + electronic_correction(Z, e)
```

where lower energy = more stable. Energies in MeV.

### 1. Weizsäcker semi-empirical mass formula

Total nuclear binding energy:

```
B(Z, N) = a_V · A
        − a_S · A^(2/3)
        − a_C · Z(Z−1) / A^(1/3)
        − a_A · (N − Z)² / A
        + δ_pair(Z, N)
```

with `A = Z + N` and coefficients from Krane / Rohlf:

| coef | value (MeV) | physical meaning |
|------|-------------|------------------|
| a_V  | 15.8        | volume term — strong-force binding |
| a_S  | 17.8        | surface — nucleons at the boundary are under-bonded |
| a_C  | 0.711       | Coulomb — proton-proton electrostatic repulsion |
| a_A  | 23.7        | asymmetry — Pauli pressure penalising |N − Z| |
| a_P  | 11.18       | pairing — even-even bonus, odd-odd penalty |

### 2. Pairing term

```
δ_pair =  +a_P / √A    if Z even and N even
       =  −a_P / √A    if Z odd  and N odd
       =   0           otherwise
```

This is why even-even nuclides dominate the chart of stable isotopes.

### 3. Magic-number shell closures

Bonus added when Z or N is in `{2, 8, 20, 28, 50, 82, 126, 184}`:

```
magic_bonus = 1.5  per magic count
            + 3.0  extra if doubly magic
```

These correspond to closed nuclear shells in the shell model (Mayer / Jensen). The doubly-magic extra explains why ⁴He, ¹⁶O, ⁴⁰Ca, ⁴⁸Ca, ²⁰⁸Pb sit at peaks of the binding curve.

### 4. Drip-line constraints

A nucleus is rejected if it falls outside the rough valley of stability:

| Z range  | required N/Z |
|----------|--------------|
| 1–8      | 0.4 .. 2.5   |
| 9–20     | 0.8 .. 1.8   |
| 21–82    | 1.0 .. 1.8   |
| 83+      | 1.3 .. 1.7   |

### 5. Fission-barrier cutoff

```
if Z² / A > 49 : nucleus unbound
```

Above this threshold the Coulomb repulsion overwhelms the surface tension and the nucleus undergoes spontaneous fission (Bohr–Wheeler criterion).

### 6. Bound-test composite

A configuration is admitted to the population only if:

- `Z ≥ 1` and `Z ≤ 138`
- `B(Z, N) + magic_bonus > 0`
- `B/A ≥ 1.0 MeV` (anything weaker is essentially a resonance, not a nucleus)
- passes drip-line and fission tests

### 7. Electronic shell closures

A small bonus is given when the electron count is a noble-gas closure: `e ∈ {2, 10, 18, 36, 54, 86, 118}`. Non-neutral atoms pay a quadratic penalty in the charge `q = Z − e`.

### 8. Search dynamics

Sixteen parallel walkers (each over `(Z, N, e)`) with Metropolis acceptance:

```
p_accept = min(1, exp(−ΔE / T))    T defaults to 2.0 MeV
```

Mutation operators per step (uniform random choice):

- Z, e ← Z ± 1, e ± 1  (move along the periodic table, neutral)
- N ← N ± 1            (try a different isotope)
- e ← e ± 1            (try an ion)
- Z, N, e all ± 1      (an α-step in mass space)
- random jump          (escape local minima)

---

## Layer 2 — molecular laws (file: `mol.c`)

The total molecular energy:

```
E_mol = Σ atom_self_energy
      + Σ bond_energy(Z₁, Z₂, order)
      + Σ ionic_lattice_bonus           [saturated per atom]
      + Σ valence_penalty               [octet + hypervalence cap]
      + 400 · net_charge²
      + Σ ring_strain
      + 1500 · is_noble_gas             [per noble-gas-touching bond]
      + 90 · n_atoms + 25 · max(0, n−8)² [size/entropy]
      + 1500 · is_disconnected
```

Units kJ/mol. Lower = more stable. Each term in detail:

### 1. Atom self-energy (`atom_self_energy` in `chem.c`)

```
e_self  = 2000   if nucleus unstable     (from layer-1 test)
        + 300 · q²                       (per-atom charge cost)
```

### 2. Bond energy (`bond_energy` in `chem.c`)

A lookup table of measured bond dissociation enthalpies (Pauling / NIST) for the common pairs — H-H, C-C, C=C, C≡C, C-O, C=O, N-N, N=N, N≡N, O-H, C-H, etc.

For atom pairs missing from the table:

```
bond_kJ = 240 · order + 96 · |EN(Z₁) − EN(Z₂)|
```

The polar bonus on the second term comes from Pauling's electronegativity model: bonds between dissimilar atoms are more stabilising than between like atoms.

Returned negative (stabilising).

### 3. Ionic lattice bonus

When `|EN(Z₁) − EN(Z₂)| > 1.7`, the bond is essentially ionic and the dominant contribution is the Madelung lattice energy, not the covalent bond. We add:

```
ionic_extra = −400 · (ΔEN − 1.7)
```

**Saturation:** each atom has a fixed number of ionic credits equal to its typical valence — so a Na⁺ can claim one ionic bonus, Ca²⁺ can claim two, etc. This prevents the search from stacking infinite ionic bonds into a megasalt.

The 1.7 threshold comes from the common empirical rule that ΔEN > 1.7 marks the covalent → ionic transition (Pauling).

### 4. Octet / valence penalty

Each atom has a preferred valence `v_typ(Z)` (1 for H, 4 for C, 3 for N, 2 for O, etc.) derived from group position:

```
v_typ(Z) = min(group, shell_capacity − group)
```

Soft penalty:

```
+180 · (have − v_typ)²
```

Hard hypervalence cap `v_max(Z)`:

- period 1–2: `v_max = v_typ`
- period 3+: hypervalent expansion allowed (S up to 6, P up to 5, Cl/Br up to 7)

If `have > v_max`:

```
+700 · (have − v_max)²
```

These two terms together encode the octet rule and its principled exceptions.

### 5. Net charge balance

```
+400 · (Σ q_i)²
```

Free-floating ions cost energy. A neutral molecule is the ground state.

### 6. Ring strain

For each ring bond, the shortest cycle through that bond is found by BFS. Then:

| ring size | energy |
|-----------|--------|
| 3 | +250 |
| 4 | +120 |
| 5 |  +25 |
| 6 |  −60 |
| 7 |  +20 |

The 6-ring bonus is a small aromatic-like reward (Hückel's rule, simplified). Three- and four-rings are penalised for angle strain (cyclopropane / cyclobutane).

### 7. Noble-gas isolation

Any bond involving He, Ne, Ar, Kr, Xe, Rn adds +1500. These atoms exist as monoatomic gases; the search shouldn't bond them.

### 8. Size / entropy

```
+90 · n_atoms
+25 · max(0, n_atoms − 8)²
```

Real chemistry penalises giant clusters via translational entropy and aggregation kinetics. Without this term the MCMC builds 64-atom megaclusters that exploit unbounded bond-energy accumulation. The 8-atom soft floor is the typical small-molecule regime.

### 9. Connectivity

Disconnected components cost +1500. The walker can still split a molecule, but only as an intermediate.

### 10. Search dynamics

Single walker, Metropolis on the population:

```
parent  ← random pick from current population
cand    ← copy of parent
        ← apply 1–3 random mutation operators
ΔE      ← E(cand) − E(parent)
accept  ← (ΔE ≤ 0) ∨ (urand() < exp(−ΔE / T))
```

Mutation operators (`mutate.c`):

- add atom (with one anchor bond)
- remove atom
- swap atom (change Z of an existing position)
- add bond
- remove bond
- bump bond order ±1
- mutate nucleus (tweak N or e)

Default `T = 600 kJ/mol`. Lower T = sharper convergence to known chemistry. Higher T = more exotic exploration.

---

## Constants summary

| Where | Symbol | Value | Source |
|-------|--------|-------|--------|
| nuclear.c | a_V | 15.8 MeV | Krane, *Introductory Nuclear Physics* |
| nuclear.c | a_S | 17.8 MeV | Krane |
| nuclear.c | a_C | 0.711 MeV | Krane |
| nuclear.c | a_A | 23.7 MeV | Krane |
| nuclear.c | a_P | 11.18 MeV | Rohlf |
| nuclear.c | magic single | 1.5 MeV | tuned shell-correction estimate |
| nuclear.c | magic double | +3.0 MeV | tuned |
| nuclear.c | fission cutoff | Z²/A < 49 | Bohr–Wheeler |
| chem.c    | bond table | kJ/mol | Pauling / NIST tabulated BDE |
| chem.c    | EN-fallback base | 240 · order | empirical fit |
| chem.c    | EN-fallback polar | 96 · ΔEN | Pauling polar model |
| mol.c     | ionic threshold | ΔEN > 1.7 | Pauling rule |
| mol.c     | ionic bonus | 400 · (ΔEN − 1.7) | tuned to put NaCl, CaO in correct wells |
| mol.c     | octet penalty | 180 · d² | tuned |
| mol.c     | hypervalence | 700 · over² | hard cap, tuned |
| mol.c     | charge balance | 400 · q² | tuned |
| mol.c     | ring strain 3,4 | 250, 120 | matches cyclopropane / cyclobutane strain energies |
| mol.c     | aromatic-6 bonus | −60 | small Hückel hint |
| mol.c     | size linear | 90 · n | translational entropy proxy |
| mol.c     | size quadratic | 25 · (n−8)² above 8 | anti-runaway |
| mol.c     | noble gas | +1500 | hard exclusion |

---

## References

- Krane, K. *Introductory Nuclear Physics*, Wiley, 1988 — SEMF
- Rohlf, J. *Modern Physics from α to Z⁰*, Wiley, 1994 — pairing
- Mayer & Jensen, *Elementary Theory of Nuclear Shell Structure*, 1955 — magic numbers
- Bohr & Wheeler, *Phys. Rev.* 56, 426 (1939) — fission barrier
- Pauling, L. *The Nature of the Chemical Bond*, Cornell UP, 1960 — electronegativity, ionic threshold
- Hückel, E. *Z. Physik* 70, 204 (1931) — aromatic stability
- NIST Chemistry WebBook — tabulated bond dissociation energies
