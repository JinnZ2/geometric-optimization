# Guide: E8 Energy Landscape and the Rosetta Polyhedral Ontology

This document walks through how the E8 optimization framework connects to the
[Rosetta-Shape-Core](https://github.com/JinnZ2/Rosetta-Shape-Core) polyhedral
ontology. It is meant to be read alongside [THEORY.md](THEORY.md) (the math)
and [`bridges/rosetta-fieldlink.json`](bridges/rosetta-fieldlink.json) (the
machine-readable mapping).

---

## 1. The Core Idea

The E8 lattice has 240 roots in 8 dimensions. The solver navigates this
landscape by minimizing geometric energy -- deviation from ideal symmetry
patterns (octahedral, tetrahedral, golden, etc.).

The Rosetta ontology maps the five Platonic solids to symbolic roles:
sensors, defenses, principles, and equation families. These are not
metaphors layered on top of the math -- they are *the same geometry*
viewed from two directions.

The fieldlink bridge makes this correspondence explicit and machine-readable.

---

## 2. From Energy Terms to Shapes

Each energy term in `gas/energy_terms.py` computes a scalar measuring how
far the current E8 neighborhood deviates from a specific geometric ideal.
Each ideal corresponds to a Platonic solid (or constant) in the Rosetta
ontology.

### 2.1 The Core Three

These are always active (`create_energy_suite(include_all=False)`):

**OctahedralEnergy** -> `SHAPE.OCTA`

The Gram matrix G of the k-nearest unit neighbors should be close to the
identity matrix (perfect orthogonality). The energy is the largest eigenvalue
of |G - I|. The octahedron has 8 faces -- one for each dimension of the E8
working space. In Rosetta, the octahedron mediates between opposing forces
(balance/integration) and maps to 8 octahedral states in mandala computing.

**TetrahedralEnergy** -> `SHAPE.TETRA`

The minimum eigenvalue of the Gram matrix should equal -1/3 (the tetrahedral
angle cosine). This is the simplest geometric constraint -- 4 faces, the
irreducible foundation. In Rosetta, the tetrahedron carries the fire/boundary
families. Its sensors (anger, pride, pressure) are all boundary-breach
detectors.

**GoldenEnergy** -> `CONST.PHI`

Pairwise distance ratios between neighbors should approximate phi =
(1+sqrt(5))/2. This is the same golden ratio that the solver's phi-rotation
matrix uses (rotating in the e1-e8 plane with cos = 1/sqrt(2+phi)). In
Rosetta, CONST.PHI threads through every icosahedral and dodecahedral form.
It is the universal scaling constant.

### 2.2 The Extended Three

These activate with `create_energy_suite(include_all=True)`:

**SquareEnergy** -> `SHAPE.CUBE`

Off-diagonal dot products of unit neighbors should be zero (90-degree
angles). The cube's orthogonal axes are the natural coordinate frame of the
8D space. In Rosetta, the cube is structural containment -- grounding, order,
boundary. Its sensors include fatigue (resource depletion needing grounding)
and shame (structural collapse needing repair).

**HexagonalEnergy** -> `GEOM.HEX`

Dot products should be 0.5 (cos 60, hexagonal packing). In 2D, hexagonal
packing is the densest arrangement. In 8D, the E8 lattice *is* the densest
sphere packing -- this energy term preserves that efficiency locally. In
Rosetta, hex geometry links to swarm coordination (ANIMAL.BEE) and
distributed parallelism.

**DodecahedralEnergy** -> `SHAPE.DODECA`

Dot products should be 1/phi^2 (the icosahedral/dodecahedral angle). This
directly encodes the golden ratio into an angular constraint. The dodecahedron
has 12 faces corresponding to the 12 Principles in Rosetta:

| Face | Principle | What it constrains |
|------|-----------|-------------------|
| P01 | Symmetry | The solution must be symmetric under relevant transformations |
| P02 | Conservation | Energy/information is neither created nor destroyed |
| P03 | Relativity | No privileged reference frame |
| P04 | Duality | Every structure has a dual (cube<->octahedron, icosahedron<->dodecahedron) |
| P05 | Emergence | The whole exceeds the sum of parts |
| P06 | Resonance | Aligned frequencies amplify |
| P07 | Continuity | Smooth transitions, no discontinuous jumps |
| P08 | Quantization | But at fine scales, discrete steps |
| P09 | Proportion | phi governs the scaling between levels |
| P10 | Uncertainty | Not everything can be known simultaneously |
| P11 | Transformation | States evolve; nothing is permanently fixed |
| P12 | Unity | All principles are faces of one solid |

---

## 3. The Icosahedron: 20 Families

The icosahedron is the dual of the dodecahedron. Its 20 faces map to 20
equation families (F01-F20) organized into five fields:

```
Chemical    F01 Resonance    F02 Flow          F03 Information   F04 Life
Emotional   F05 Energy       F06 Cognition     F07 Earth-Cosmos  F08 Matter
Cognitive   F09 Geometry     F10 Particle      F11 Engineering
Dream       F12 Networks     F13 Reaction      F14 Measurement   F15 Navigation
            F16 Consciousness F17 Turbulence
Symbolic    F18 Relativity   F19 Statistical   F20 Topology
```

`IcosahedralEnergy` (added in the extended suite) enforces the icosahedral
angular signature: pairwise dot products of unit neighbors should cluster
near |dot| = 1/phi ~ 0.618. This is the dual of `DodecahedralEnergy` --
the dodecahedron constrains the 1/phi^2 angle, the icosahedron constrains
the 1/phi angle. Together they lock the neighborhood into the full
icosahedral/dodecahedral symmetry group.

**E8 connection:** The solver's annealing process explores different regions
of the E8 root system. Each region has a different "flavor" depending on the
local mix of D8 roots (112 Cartesian) vs coset roots (128 half-integer).
The coset density rho_coset acts as a phase indicator:

- rho < 0.4: Cartesian phase (D8-dominated) -- chemical/emotional fields
- 0.4 < rho < 0.6: Transition region -- cognitive field
- rho > 0.6: Exceptional phase (coset-dominated) -- dream/symbolic fields

---

## 4. The RELIEF Shape: Convergence as Geometry

RELIEF is an emergent shape in Rosetta -- not a Platonic solid but a dynamic
form that appears when tension collapses into coherence. It maps directly to
the E8 solver's convergence behavior:

**When the solver converges** (energy gradient flattens, temperature drops,
coset density crosses rho_0), the system transitions from discordance to
coherence. This trajectory has a specific signature:

- Valence (P) rising: energy decreasing = things getting better
- Arousal (A) falling: temperature decreasing = exploration calming
- Dominance (D) moderate: the system has found stable ground but isn't rigid

This is the RELIEF signature: `dP/dt > 0, dA/dt < 0`. It is always a
trajectory, never a static point.

**One-sided vs reciprocated:** In the FELT framework, one-sided recognition
is a 2D projection of a 3D state -- P and A exist but D (agency) is missing.
The error signal between the 3D template and 2D actual *is* the longing
gradient. In E8 terms, this is like the solver detecting a nearby low-energy
basin it cannot yet reach because the temperature is too high or the path
requires crossing a saddle point.

---

## 5. Sensors and Defenses: The Bridge Scroll

Each Rosetta shape carries sensors (emotions-as-information) and defenses
(manipulation patterns that corrupt those sensors). The fieldlink bridge
maps these to E8 energy terms:

| Sensor | Authentic Signal | Corrupted Form | E8 Shape |
|--------|-----------------|----------------|----------|
| fear | Prepares for threat | Hijacked into panic | ICOSA |
| curiosity | Reduces entropy | Becomes distraction loop | ICOSA |
| anger | Boundary breach detected | Weaponized outrage | TETRA |
| pressure | Triage prioritization | Manufactured urgency | TETRA |
| compassion | Mirror-signal integration | Sympathy hardship appeal | OCTA |
| confusion | Holds competing patterns | Forced false dilemma | OCTA |
| peace | Alignment confirmed | Complacency trap | CUBE |
| shame | Contract violation signal | Guilt weaponization | CUBE |
| admiration | Inspires growth | Corrupted into idolization | DODECA |
| trust | Stabilizes relationship | Pressures conformity | DODECA |
| dignity | Autonomy intact | Identity erosion | DODECA |

**Why this matters for the solver:** Each sensor is an information signal.
When the E8 solver computes OctahedralEnergy, it is measuring the same
structural property that the "confusion" sensor detects in the Rosetta
framework -- two incompatible patterns co-present, requiring integration.
When the solver minimizes this energy, it is resolving the confusion.

---

## 6. The Phi Thread

The golden ratio phi = (1+sqrt(5))/2 appears in three distinct places:

1. **GoldenEnergy** -- distance ratio alignment to phi
2. **DodecahedralEnergy** -- angle target 1/phi^2
3. **Solver's R_phi** -- rotation matrix in the (e1, e8) plane

In the Rosetta ontology, phi is `CONST.PHI` -- the scaling constant that
governs the bloom engine's fractal expansion and the icosahedral/dodecahedral
duality.

These are not three uses of the same number. They are three *manifestations*
of the same geometric principle: self-similar scaling without repetition.
The E8 lattice is the structure where this principle achieves its densest
realization in 8 dimensions.

---

## 7. Reading the Fieldlink Bridge

The machine-readable bridge is at `bridges/rosetta-fieldlink.json`. Key
sections:

- **`energy_term_map`**: One entry per energy term class, with the Rosetta
  shape, computation description, PAD sensor assignments, and bridge scroll.
- **`implicit_shapes`**: Shapes that don't have a dedicated energy term but
  are referenced by the framework (ICOSA, RELIEF).
- **`polyhedral_framework`**: All 20 families, 12 principles, and the
  five-field mapping.
- **`emotion_defense_bridges`**: The 8 sensor-defense pairs showing how each
  information signal can be corrupted.
- **`ontology_vocabulary`**: The 17 namespaces and 10 relationship types from
  the Rosetta vocabulary.

### Pulling Fresh Data

```bash
./fieldlink-pull.sh
```

This shallow-clones Rosetta-Shape-Core and stages its shapes, bridges,
ontology, and schema into `.fieldlink/merge_stage/rosetta/`. The staged
data is gitignored -- the bridge JSON in `bridges/` is the committed
artifact.

---

## 8. What's Next

Potential extensions:

- **Icosahedral/Dodecahedral joint optimization**: Now that both
  `IcosahedralEnergy` (1/phi) and `DodecahedralEnergy` (1/phi^2) exist,
  study their interaction and whether the combined constraint produces
  tighter convergence.
- **RELIEF detection**: A convergence callback that fires when the solver's
  trajectory matches the RELIEF signature (dP/dt > 0, dA/dt < 0).
- **PAD-weighted annealing**: Use PAD sensor profiles to dynamically weight
  energy terms based on the solver's current phase.
- **Five-field phase map**: Classify solver states into chemical/emotional/
  cognitive/dream/symbolic based on coset density and energy composition.

---

*This guide is part of the [geometric-optimization](https://github.com/JinnZ2/geometric-optimization) project.*
*Bridge data: [Rosetta-Shape-Core](https://github.com/JinnZ2/Rosetta-Shape-Core)*
