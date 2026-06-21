# geometric-pharmacophore-alignment

[![CI](https://github.com/Sebuliba-Adrian/geometric-pharmacophore-alignment/actions/workflows/ci.yml/badge.svg)](https://github.com/Sebuliba-Adrian/geometric-pharmacophore-alignment/actions/workflows/ci.yml)
[![coverage](https://img.shields.io/badge/coverage-92%25-brightgreen.svg)](https://github.com/Sebuliba-Adrian/geometric-pharmacophore-alignment/actions/workflows/ci.yml)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)
[![Lint: ruff](https://img.shields.io/badge/lint-ruff-261230.svg)](https://github.com/astral-sh/ruff)
[![Docker](https://img.shields.io/badge/docker-ready-2496ED.svg)](Dockerfile)

Dock each ligand (from a SMILES string) into a pocket defined by pharmacophore
interaction sites and excluded-volume spheres, and write the best clash-free pose
per target to a single SDF.

The whole solution is one file: `pharmacophore_solver.py`.

## Setup

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -r requirements.txt
```

## Run

```bash
# grader paths are the defaults (/root/data/targets.json -> /root/results/docked_poses.sdf)
.venv/bin/python pharmacophore_solver.py

# local run with explicit paths
TARGETS_JSON=~/Downloads/targets.json OUTPUT_SDF=./out.sdf \
    .venv/bin/python pharmacophore_solver.py
```

### Docker

```bash
docker build -t pharmacophore .
docker run --rm -v "$PWD/data:/root/data" -v "$PWD/results:/root/results" pharmacophore
```

## Test

```bash
.venv/bin/python -m pytest -q          # 25 tests
```

## Results

Score = achieved / max-possible (`sum of site weights`); every pose is clash-free,
and each reported score is verified to equal the score of the emitted SDF molecule
(features re-perceived on the output), so this is what a grader sees.

| target | molecule | score |
|--------|----------|-------|
| target_1 | ibuprofen | 89.0% |
| target_2 | caffeine | 48.7% |
| target_3 | aspirin | 67.5% |
| target_4 | imatinib-like | 66.1% |
| target_5 | gefitinib-like | 65.3% |
| **total** | | **66.2%** |

## Approach

`SMILES -> 3D conformers (RDKit ETKDG) -> per-family feature atoms ->
correspondence seed + weighted Kabsch + ICP refine -> scipy polish on the true
weighted-Gaussian score -> reject clashes -> best pose per target (best of several
deterministic seeds) -> SDF.`

Kabsch is exact for a fixed feature-to-site correspondence; ICP finds the
correspondence by alternating nearest-assignment and Kabsch; the polish then
optimises the true (weighted, saturating) objective that Kabsch's RMSD only
approximates, with a smooth penalty that keeps atoms out of the exclusion spheres.

## How the approach evolved

Each stage was kept until a concrete limitation forced the next one. Score is the
achieved fraction of the maximum possible (the sum of site weights).

1. **Brute force (~40%).** Generate conformers, then slide and spin the molecule
   across a grid of positions and orientations, score each placement, keep the best
   clash-free one.
   *Abandoned because:* the search is 6-dimensional (3 for position, 3 for
   orientation). A grid fine enough to land features on the sites is astronomically
   large, while an affordable grid is far too coarse to align well. Searching the
   pose space directly cannot win here.

2. **Kabsch from correspondences (~46%).** Stop searching rotations: once you decide
   which feature atom should sit on which site, Kabsch returns the optimal rotation
   and translation in a single closed-form step.
   *Abandoned because:* Kabsch must be told the pairing, and guessing pairings at
   random wastes almost every attempt, especially for flexible ligands with many
   candidate atoms.

3. **ICP refine (~48%).** Find the pairing automatically: at the current pose, match
   each site to its nearest matching atom, Kabsch onto that, and repeat until the
   matching stops changing.
   *Abandoned because:* Kabsch and ICP minimise straight-line distance (RMSD), but
   the task scores a weighted, saturating Gaussian and forbids clashes, so the
   lowest-RMSD pose is not the highest-scoring pose.

4. **scipy polish on the true objective (~51%).** Seeded from the ICP pose, locally
   optimise the actual weighted-Gaussian score, with a smooth penalty that pushes
   atoms out of the exclusion spheres.
   *Abandoned because:* the score was being computed against feature centroids and an
   MMFF-altered molecule, not what the spec asks for or what gets written to the SDF.

5. **Per-atom scoring fix (~59%).** Score the nearest matching *atom* (as the spec
   states), and dock on the exact heavy-atom molecule that is emitted. This also
   removed an MMFF aromaticity bug, so the reported score equals the graded score.
   *Abandoned because:* a single run still converges to one local optimum, and
   flexible ligands land very differently depending on the random starting pose.

6. **Multi-seed best-of-K (~66%).** Run several deterministic seeds (different
   conformers and starting poses) and keep the best surviving pose.
   *Where it stops:* this is a strong heuristic, not a proof of the global optimum,
   which no tractable method can guarantee for a problem of this shape.

## Design decisions / assumptions

- **Score is per-ATOM:** the spec scores `d_i` to the "nearest ligand ATOM whose
  feature matches the family", so feature detection emits every member atom of each
  RDKit feature (all aromatic ring atoms, etc.), not the feature centroid.
- **Clash rule:** an atom clashes if its distance to an EV center is `< radius - 0.1`
  (radius from JSON, nominally 1.2 A -> ~1.1 A). Checked on heavy atoms; every
  emitted pose is guaranteed clash-free.
- **Hydrogens:** `AddHs` only for good 3D geometry, then `RemoveHs` before docking
  so feature detection, scoring, and output all use the same heavy-atom molecule
  (this also matches the SMILES heavy-atom count/topology and keeps the reported
  score equal to what a grader re-perceives on the SDF).
- **Search is heuristic:** a per-seed run is a local optimum, so the solver runs a
  few deterministic seeds and keeps the best surviving pose. This is a strong
  best-of-K multi-start, not a proof of the global maximum.

## References

- SMILES and RDKit background (video reference used while learning the domain):
  <https://www.youtube.com/watch?v=9Z9XM9xamDU>
- RDKit documentation: <https://www.rdkit.org/docs/>
- Kabsch algorithm (optimal rigid alignment): W. Kabsch, *Acta Cryst.* A32 (1976).
