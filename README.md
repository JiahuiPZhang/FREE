# Forcefield Reinforced Evolvable Encoder (FREE)

A physics-informed **8-dimensional physicochemical encoding pipeline** for amino acids based directly on Amber topology (`prmtop`) files.

The pipeline extracts molecular-mechanics information from Amber force-field parameters and converts each amino-acid residue into an 8-dimensional feature vector describing:

1. Electrostatics
2. Sterics
3. Dispersion interactions
4. Hydrogen-bond adaptability
5. Torsional rigidity
6. Conformational entropy

A **Terminal Effect Balancing** mechanism is included to reduce systematic differences caused by whether the amino acid is represented as an internal residue, an N-terminal residue, or a C-terminal residue.

The implementation is based on `ParmEd` and therefore obtains atomic charges, Lennard-Jones parameters, bonding topology, and dihedral parameters directly from Amber topology files.

---

## Overview

The final encoding contains eight physicochemical features:

| ID | Feature | Category | Physical meaning |
|---|---|---|---|
| F1 | Net Charge | Electrostatics | Intrinsic residue net charge |
| F2 | Charge Variance | Electrostatics | Distribution of atomic partial charges |
| F3 | van der Waals Volume Proxy | Sterics | Approximate steric volume |
| F4 | Dispersion Capacity | Dispersion | Lennard-Jones dispersion interaction capacity |
| F5 | Donor Strength | H-bond Adaptability | Electrostatic strength of hydrogen-bond donors |
| F6 | Acceptor Strength | H-bond Adaptability | Electrostatic strength of hydrogen-bond acceptors |
| F7 | Torsional Rigidity | Flexibility & Entropy | Effective torsional energy barriers |
| F8 | Entropy Proxy | Flexibility & Entropy | Number of accessible torsional minima |

The resulting representation can be written as

$$
\mathbf{x}
=
(F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8).
$$

---

# 1. Terminal Effect Balancing

A free amino acid placed at the N- or C-terminus contains additional terminal atoms and terminal charges that are absent when the same amino acid is located inside a peptide chain.

To reduce this topological bias, the pipeline automatically determines the residue topology from the `prmtop` filename and applies explicit terminal corrections.

Three supported states are defined.

| Filename pattern | Topological state | Charge bias | Boundary atoms masked |
|---|---|---:|---|
| `_caps` | Internal | 0 | None |
| `_Ncap` | C-terminal residue | -1 | `OXT` |
| `_Ccap` | N-terminal residue | +1 | `H1`, `H2`, `H3` |

### Internal residue

A topology containing

```text
_caps
