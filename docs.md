# exre — usage

## Build

```sh
make
```

Produces two binaries:

- `exre`  — layer 2, molecule discovery
- `exre1` — layer 1, element discovery

Clean with `make clean`.

## Layer 2 — molecule discovery (`exre`)

```sh
./exre [iters] [T] [seed]
```

| arg   | default | meaning |
|-------|---------|---------|
| iters | 200000  | number of MCMC steps |
| T     | 600.0   | "temperature" in kJ/mol — controls acceptance of energy-uphill mutations |
| seed  | time()  | RNG seed; pass an integer for reproducible runs |

### Examples

```sh
./exre                       # default run
./exre 500000                # half a million iterations
./exre 1000000 400 42        # cooler, reproducible
./exre 2000000 200 31337     # long, cold — sharp convergence
```

### Output

stdout shows every seed loaded, every novel formula discovered, and every new global-best molecule. Each line:

```
[novel] OMg              E=-676.5  atoms=2 bonds=1  parent=O2
[best ] CNO3F2SiCa4      E=-6479.8 atoms=12 bonds=11 parent=silane
```

A tab-separated log is written to `exre.log` with one row per discovery:

```
novel<TAB>OMg<TAB>E=-676.480<TAB>atoms=2<TAB>bonds=1<TAB>parent=O2
```

Final summary:

```
iters=...  accepted=...  rejected=...  novel=...  pop=...
global best E=...
```

### Tuning

- **High T (800+)** — wider exploration, more diverse novelties, weaker convergence to global minimum.
- **Low T (200–300)** — tighter convergence on real stable species; fewer exotic suggestions.
- **More iters** — discovery curve plateaus around 500k–1M; longer runs mostly polish the population.

### Inspecting discoveries

```sh
# top stable molecules of any size
grep -E "^(novel|best)" exre.log | sort -t= -k2 -g | head -30

# stable small molecules only (≤ 5 atoms)
grep -E "^novel" exre.log | \
  awk -F'\t' '{split($3,a,"="); gsub("atoms=","",$4); print a[2]"\t"$2"\t"$4}' | \
  sort -g | awk -F'\t' '$3+0 <= 5' | head -30
```

## Layer 1 — element discovery (`exre1`)

```sh
./exre1 [iters] [T] [seed]
```

| arg   | default  | meaning |
|-------|----------|---------|
| iters | 500000   | number of MCMC steps |
| T     | 2.0      | temperature in MeV |
| seed  | time()   | RNG seed |

### Examples

```sh
./exre1                       # default run
./exre1 1000000 1.5 42        # one million steps, cooler
./exre1 2000000 1.0 31337     # long, sharp scan
```

### Output

stdout reports:

1. **Top 30 most stable nuclides** sorted by B/A, with magic-number annotations
2. **Known elements rediscovered** (Z = 1..36) with most-stable A
3. **Superheavy candidates** (Z > 118) — the island of stability

Log: `exre_elements.log`, one row per nuclide:

```
Z<TAB>N<TAB>A<TAB>e<TAB>B/A_MeV<TAB>name-A
28<TAB>28<TAB>56<TAB>28<TAB>8.7716<TAB>Ni-56
120<TAB>184<TAB>304<TAB>120<TAB>7.1423<TAB>Unbinil-ium-304
```

Names beyond Z = 36 use IUPAC systematic naming (Unbinilium, Unbihexium, …).

### Tuning

- **T ≈ 2.0 MeV** — balanced; finds both iron peak and superheavy island
- **T ≈ 0.5 MeV** — sharper attraction to global minimum (clusters near Ni-56)
- **T ≈ 5.0 MeV** — wide exploration; surfaces more drip-line candidates

### Filtering for undiscovered candidates

```sh
# most stable nuclide per superheavy element (Z > 118)
awk -F'\t' '$1+0 > 118 {
  if ($5+0 > best[$1]) { best[$1]=$5; bestN[$1]=$2; bestA[$1]=$3 }
} END {
  for (z in best) printf "Z=%-3d  N=%-3d  A=%-3d  B/A=%.4f MeV\n",
    z, bestN[z], bestA[z], best[z]
}' exre_elements.log | sort -k2 -n
```

## File layout

```
exre/
├── exre.h        ── shared atom/bond/molecule structs, RNG
├── chem.c        ── periodic table, valence, EN, bond energies
├── mol.c         ── molecule graph ops, scoring, formula print
├── mutate.c      ── seven mutation operators
├── seeds.c       ── 27 built-in known molecules
├── main.c        ── layer-2 MCMC driver
├── nuclear.h     ── nuclear physics API
├── nuclear.c     ── Weizsäcker SEMF + magic numbers + drip lines
├── disc_elem.c   ── layer-1 MCMC driver
├── Makefile
├── README.md
├── docs.md       ── (this file)
└── details.md    ── full law reference
```

For the physical laws behind each term, see [details.md](details.md).
