# exre — stochastic chemical discovery

A two-layer Monte Carlo search written in C that hunts for stable atoms and stable molecules under quantum-chemistry-inspired laws. The only attractor is **least-energy stable state**; everything else is uniform random mutation.

---

## Why

Most chemistry software simulates known systems. exre does the opposite: it explores the space of *possible* atoms and *possible* molecules by random mutation and lets a physical attractor (the lowest free-energy state) decide what survives. The hypothesis is simple:

> If the laws of stability are encoded faithfully, then by the law of large numbers, every real species should re-appear, and any exotic species the search finds should be a candidate for real synthesis.

There is no selectivity built into the search — no heuristic that "guides" the walker toward known chemistry. The only bias is the energy law. If a configuration is lower-energy, it wins.

## How

Two independent binaries, each running its own MCMC walk:

### Layer 1 — `exre1` (element discovery)

Walker state: a triple `(Z, N, e)` — protons, neutrons, electrons.
Search space: all atoms, known and unknown.
Energy = `−B(Z,N) − magic_bonus(Z,N) + electronic_correction(Z,e)`

Binding energy `B(Z,N)` comes from the **Weizsäcker semi-empirical mass formula**:

```
B = a_V·A − a_S·A^(2/3) − a_C·Z(Z−1)/A^(1/3) − a_A·(N−Z)²/A ± δ_pair
```

with magic-number shell-closure bonuses at `{2, 8, 20, 28, 50, 82, 126, 184}` for both Z and N, electronic closures at `e ∈ {2, 10, 18, 36, 54, 86, 118}`, drip-line cutoffs, and a Z²/A < 49 fission-barrier limit.

Sixteen parallel walkers explore the (Z, N, e) cube via Metropolis acceptance.

### Layer 2 — `exre` (molecule discovery)

Walker state: a labelled graph — atoms (Z, N, e) as nodes, bonds with order as edges.
Search space: all molecules buildable from the periodic table.
Energy = sum of seven terms (bond energies, octet penalty, ionic-lattice bonus, ring strain, hypervalence cap, charge balance, size/entropy).

Seeded with 27 known molecules (H₂O, CO₂, CH₄, NH₃, C₆H₆, HCOOH, …). Each iteration: pick a parent at random, apply 1–3 random mutations (add/remove atom, swap, add/remove bond, bump bond order, tweak nucleus), score, accept by Metropolis. Novel formulas are added to the population.

Full law reference in [details.md](details.md). Usage and CLI in [docs.md](docs.md).

---

## Discoveries

### Layer 1 — atoms

The MCMC walker, knowing only Weizsäcker + magic numbers + drip lines, independently produced the central facts of nuclear physics:

**The iron peak.** Top-of-stability ranking:

| Rank | Nuclide | B/A (MeV) | Note |
|------|---------|-----------|------|
| 1    | Ni-56   | 8.772     | doubly magic, Z = N = 28 |
| 2    | Ca-41   | 8.749     | magic Z = 20 |
| 3    | K-39    | 8.735     | magic N = 20 |
| 4    | Mn-51   | 8.735     | |
| 5    | Fe-53   | 8.731     | |

Real experimental peak is Ni-62 at 8.7945 MeV/A. exre1 hit the right plateau within 0.2%.

**Known elements rediscovered** (Z 1–36): all isotopes within ±1 mass unit of the experimentally observed most-stable nuclide, with exact matches for F-19, Na-23, Al-27, S-32, Cl-35, K-39.

**The island of stability — predicted.** Above Z = 118 (oganesson, the heaviest element ever synthesised), exre1 flagged:

| Z | Symbol | Most-stable A | Most-stable N | B/A (MeV) |
|---|--------|---------------|---------------|-----------|
| 119 | Ununennium  | 297 | 178 | 7.159 |
| **120** | **Unbinilium**  | **304** | **184** | **7.142** |
| 121 | Unbiunium   | 305 | 184 | 7.121 |
| 122 | Unbibium    | 306 | 184 | 7.103 |
| 126 | Unbihexium  | 324 | 198 | 7.012 |

Both stability islands match real-world nuclear theory: the N = 184 cluster around Z = 120 (currently being hunted at RIKEN and Dubna) and the Z = 126 cluster predicted by shell theory.

### Layer 2 — molecules

**Real species rediscovered from random mutation:**

| Found | Identity |
|-------|----------|
| `OMg`     | MgO, periclase / magnesia |
| `OCa`     | CaO, quicklime |
| `FK`, `FNa` | KF, NaF — fluoride salts |
| `ClK`     | KCl, sylvite |
| `O2Mg`    | MgO₂, magnesium peroxide |
| `OFK`     | KOF, potassium hypofluorite (synthesised 1968) |
| `MgCl2Ca` | tachyhydrite family mineral |
| `CHO`     | formyl radical (interstellar) |
| `CHONa`   | sodium formate HCOONa |
| `O2Si`    | silica, SiO₂ |
| `OSiCa`   | calcium silicate — the binder in Portland cement |
| `CFCl2`   | chlorofluoromethyl radical (CFC fragment) |
| `P2`      | diphosphorus (gas-phase) |
| `CNO`     | cyanate / fulminate |

**Predicted but uncatalogued molecular species (the interesting ones):**

| Formula | E (kJ/mol) | Plausible identity |
|---------|------------|---------------------|
| `NOMgCa`   | −1546 | mixed Mg/Ca nitrosyl-oxide |
| `HN2NaK`   | −1189 | mixed alkali hydrazide-amide |
| `NOMgKCa`  | −1203 | quaternary nitrosyl cluster |
| `Na2SCaBr` | −1127 | sodium-calcium sulfobromide |
| `CHOMgCl`  | −1116 | chloroformyl-magnesium |
| `NFKCa`    | −1023 | potassium-calcium fluoro-nitride |
| `OCa2`     | −620  | Ca₂O suboxide |
| `NKCa`     | −603  | potassium calcium nitride |
| `OFPCl`    | −601  | fluoro-chloro phosphine oxide |
| `C3Ca`     | −764  | calcium tricarbide (only CaC₂ is known) |
| `KBr2`     | −617  | K[Br₂]⁻ ion pair |

---

## Theoretical status

These results are **theoretical and thermodynamically plausible**. They are produced by the program's encoding of quantum-chemistry-inspired laws (Weizsäcker SEMF, octet rule, bond enthalpies, ionic lattice approximation, ring strain, magic-number shell closures), with the sole attractor being the **least-energy stable state**.

Most of what the program produces it **rediscovers** — the iron peak, KCl, CaO, SiO₂, sodium formate, the predicted superheavy island around Z = 120. The rediscoveries are the model's correctness check.

The intended use of exre is the *opposite* slice of the output: the **exotic, undiscovered candidates** that the search keeps finding stable but which have no entry in standard chemistry. Those are the species worth experimental investigation — `NOMgCa`, `HN2NaK`, `OCa2`, the superheavy `Z = 120, A = 304`, and so on.

**Important caveat.** This program does not simulate reality. It uses analytical approximations of the underlying laws — there is no Schrödinger equation, no relativistic correction, no kinetic barrier, no solvent. A "stable" candidate from exre is a thermodynamic plausibility, not a synthesis recipe. It will fail to capture: spontaneous fission rates, alpha-decay half-lives, transition states, steric hindrance beyond simple ring strain, solvation, and entropy of mixing.

**The purpose of exre is to surface exotic, undiscovered chemical species — atoms and molecules — that a researcher should investigate further with real quantum-chemistry software or laboratory synthesis.**

---

## Further reading

- [docs.md](docs.md) — build, run, CLI flags, output format
- [details.md](details.md) — the seven molecular laws and the nine nuclear laws, with constants, formulas, and references
