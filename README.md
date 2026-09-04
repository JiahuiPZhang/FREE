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
| --- | --- | --- | --- |
| F1 | Net Charge | Electrostatics | Intrinsic residue net charge |
| F2 | Charge Variance | Electrostatics | Distribution of atomic partial charges |
| F3 | van der Waals Volume Proxy | Sterics | Approximate steric volume |
| F4 | Dispersion Capacity | Dispersion | Lennard-Jones dispersion interaction capacity |
| F5 | Donor Strength | H-bond Adaptability | Electrostatic strength of hydrogen-bond donors |
| F6 | Acceptor Strength | H-bond Adaptability | Electrostatic strength of hydrogen-bond acceptors |
| F7 | Torsional Rigidity | Flexibility and Entropy | Effective torsional energy barriers |
| F8 | Entropy Proxy | Flexibility and Entropy | Number of accessible torsional minima |

The resulting representation can be written as:

```math
\mathbf{x}=(F_1,F_2,F_3,F_4,F_5,F_6,F_7,F_8)
```

---

# 1. Terminal Effect Balancing

A free amino acid placed at the N- or C-terminus contains additional terminal atoms and terminal charges that are absent when the same amino acid is located inside a peptide chain.

To reduce this topological bias, the pipeline automatically determines the residue topology from the `prmtop` filename and applies explicit terminal corrections.

Three supported states are defined:

| Filename pattern | Topological state | Charge bias | Boundary atoms masked |
| --- | --- | ---: | --- |
| `_caps` | Internal | 0 | None |
| `_Ncap` | C-terminal residue | -1 | `OXT` |
| `_Ccap` | N-terminal residue | +1 | `H1`, `H2`, `H3` |

### Internal residue

A topology containing:

```text
_caps
```

is interpreted as an internally capped amino acid.

No terminal charge correction or boundary-atom masking is applied.

### C-terminal residue

A topology containing:

```text
_Ncap
```

is interpreted as an amino acid whose N-terminus is capped by ACE while its C-terminus remains free.

The free C-terminal carboxyl group contributes an approximately `-1` systematic charge.

Therefore:

```text
charge_bias = -1.0
```

The terminal atom:

```text
OXT
```

is excluded from boundary-sensitive features.

### N-terminal residue

A topology containing:

```text
_Ccap
```

is interpreted as an amino acid whose C-terminus is capped by NME while its N-terminus remains free.

The free N-terminal ammonium group contributes an approximately `+1` systematic charge.

Therefore:

```text
charge_bias = +1.0
```

The terminal protons:

```text
H1
H2
H3
```

are excluded from boundary-sensitive features.

---

# 2. Eight-Dimensional Encoding

## F1 — Intrinsic Net Charge

### Category

**Electrostatics**

The raw residue charge is calculated from the sum of all Amber atomic partial charges:

```math
Q_{\mathrm{raw}}=\sum_i q_i
```

The terminal contribution is then explicitly removed:

```math
F_1=Q_{\mathrm{balanced}}
=Q_{\mathrm{raw}}-Q_{\mathrm{terminal\ bias}}
```

In the implementation:

```python
raw_charge = sum(atom.charge for atom in residue.atoms)
balanced_charge = raw_charge - charge_bias
```

The terminal charge biases are:

```text
Internal residue:      0
N-terminal residue:   +1
C-terminal residue:   -1
```

The final `F1` value is the balanced charge.

### Physical interpretation

`F1` describes the intrinsic electrostatic charge state of the amino-acid residue after explicitly removing the systematic charge introduced by free peptide termini.

---

## F2 — Charge Variance

### Category

**Electrostatics**

The second electrostatic feature is the variance of atomic partial charges:

```math
F_2=\operatorname{Var}(q_1,q_2,\ldots,q_N)
```

It is calculated as:

```python
charges = [
    atom.charge
    for atom in residue.atoms
    if atom.name not in boundary_atoms
]

F2 = np.var(charges)
```

Boundary atoms associated with terminal effects are removed before calculating the variance.

For an N-terminal residue:

```text
H1, H2, H3
```

are removed.

For a C-terminal residue:

```text
OXT
```

is removed.

### Physical interpretation

`F2` measures the heterogeneity of the residue's atomic charge distribution and acts as a descriptor of intrinsic residue polarity.

A residue containing strongly separated positive and negative partial charges produces a larger charge variance than a residue with a more uniform charge distribution.

---

# 3. Steric and Dispersion Features

## F3 — van der Waals Volume Proxy

### Category

**Sterics and Dispersion**

The steric size of the residue is approximated using the Lennard-Jones sigma parameters of all atoms:

```math
F_3=\sum_i \sigma_i^3
```

Implementation:

```python
F3 = sum(atom.sigma**3 for atom in residue.atoms)
```

### Physical interpretation

The Lennard-Jones parameter sigma represents a characteristic atomic length scale.

Its cube provides a volume-like quantity:

```math
V_i^{\mathrm{proxy}}=\sigma_i^3
```

The sum over all atoms therefore serves as a simple proxy for the effective van der Waals volume of the residue.

This is a force-field-derived descriptor rather than a geometrically calculated molecular volume.

---

## F4 — Dispersion Capacity

### Category

**Sterics and Dispersion**

The dispersion feature is obtained by summing the Lennard-Jones epsilon parameters over heavy atoms:

```math
F_4=\sum_{i\in\mathrm{heavy\ atoms}}\epsilon_i
```

Implementation:

```python
F4 = sum(
    atom.epsilon
    for atom in residue.atoms
    if atom.atomic_number > 1
)
```

Hydrogen atoms are excluded.

### Physical interpretation

The Lennard-Jones potential is:

```math
V_{\mathrm{LJ}}(r)
=
4\epsilon
\left[
\left(\frac{\sigma}{r}\right)^{12}
-
\left(\frac{\sigma}{r}\right)^6
\right]
```

The parameter epsilon determines the depth of the Lennard-Jones potential well.

The sum of heavy-atom epsilon values is therefore used as a residue-level descriptor of dispersion interaction capacity.

---

# 4. Hydrogen-Bond Adaptability

## F5 — Donor Strength

### Category

**Hydrogen-bond adaptability**

Hydrogen-bond donor strength is calculated from polar hydrogen atoms.

A hydrogen is treated as a donor hydrogen when it is covalently bonded to:

```text
N
O
S
```

corresponding to atomic numbers:

```python
[7, 8, 16]
```

For each such hydrogen, its Amber partial charge is added:

```math
F_5=\sum_{H_{\mathrm{polar}}}q_H
```

Implementation concept:

```python
donor_strength = 0.0

for atom in residue.atoms:

    if atom.name in boundary_atoms:
        continue

    if atom.atomic_number == 1:

        for bond in atom.bonds:

            partner = (
                bond.atom1
                if bond.atom2 == atom
                else bond.atom2
            )

            if partner.atomic_number in [7, 8, 16]:
                donor_strength += atom.charge
                break
```

For N-terminal residues, the additional terminal protons:

```text
H1
H2
H3
```

are masked.

### Physical interpretation

Polar hydrogen atoms generally carry positive partial charges.

`F5` therefore measures the total electrostatic contribution from hydrogen atoms capable of participating in hydrogen-bond donation.

---

## F6 — Acceptor Strength

### Category

**Hydrogen-bond adaptability**

Potential acceptor atoms are identified as nitrogen or oxygen atoms:

```python
atomic_number in [7, 8]
```

The absolute values of their atomic partial charges are summed:

```math
F_6=\sum_{i\in\{N,O\}}\lvert q_i\rvert
```

Implementation:

```python
acceptor_strength = 0.0

for atom in residue.atoms:

    if atom.name in boundary_atoms:
        continue

    if atom.atomic_number in [7, 8]:
        acceptor_strength += abs(atom.charge)
```

For C-terminal residues, the terminal atom:

```text
OXT
```

is masked.

### Physical interpretation

`F6` represents the total electrostatic magnitude associated with nitrogen and oxygen atoms that can act as hydrogen-bond acceptors.

The use of absolute partial charge prevents negatively charged acceptor atoms from cancelling one another in the residue-level sum.

---

# 5. Flexibility and Conformational Entropy

The last two features are calculated from Amber bonding and torsional parameters.

The procedure first identifies the residue's `CA` atom and calculates the topological distance of every atom from `CA`.

For an atom `i`:

```math
d_i=
\text{minimum number of covalent bonds between CA and atom }i
```

For a bond connecting atoms `i` and `j`, the bond depth is:

```math
d_{ij}=\max(d_i,d_j)
```

Only bonds satisfying the following conditions are analyzed:

- Both atoms belong to the target residue.
- Both atoms are heavy atoms.
- The bond is not a pure backbone bond.
- The bond has associated Amber dihedral terms.

The following atoms are treated as backbone atoms:

```text
N
CA
C
O
H
HA
```

Bonds connecting two backbone atoms are skipped.

---

## F7 — Torsional Rigidity

### Category

**Flexibility and Entropy**

For every selected bond, all Amber dihedral terms whose central bond corresponds to that bond are collected.

For a dihedral term, the torsional energy is evaluated as:

```math
E(\phi)=k_\phi\left[1+\cos(n\phi-\delta)\right]
```

If multiple Amber dihedral terms contribute to the same bond, they are summed:

```math
E_{\mathrm{bond}}(\phi)
=
\sum_j
k_{\phi,j}
\left[
1+\cos(n_j\phi-\delta_j)
\right]
```

The torsional profile is evaluated at 360 uniformly spaced angles between:

```math
0\leq\phi<2\pi
```

For each bond, the effective torsional barrier is:

```math
\Delta E=E_{\max}-E_{\min}
```

The contribution of the bond to `F7` is weighted by its topological distance from the alpha carbon:

```math
R_k=\frac{\Delta E_k}{d_k}
```

The complete feature is:

```math
F_7=\sum_k\frac{\Delta E_k}{d_k}
```

Implementation:

```python
delta_E = np.max(energy_profile) - np.min(energy_profile)

rigidity_F7 += (
    1.0 / bond_depth
) * delta_E
```

### Physical interpretation

`F7` represents an effective residue-level torsional rigidity.

A bond with a large torsional energy barrier contributes strongly to the rigidity score.

The contribution is additionally weighted by:

```math
\frac{1}{d_k}
```

Therefore, torsional barriers closer to the alpha carbon contribute more strongly than those located farther away in the residue topology.

---

## F8 — Entropy Proxy

### Category

**Flexibility and Entropy**

The same torsional energy profile used for `F7` is also used to estimate the number of accessible conformational minima.

For each rotatable bond, the number of local minima is counted:

```math
m_k=
\text{number of local minima in }E_k(\phi)
```

A point `i` is considered a local minimum when:

```math
E_i<E_{i-1}
```

and:

```math
E_i<E_{i+1}
```

Periodic boundary conditions are used for the angular energy profile.

At least one minimum is enforced:

```math
m_k\geq 1
```

The entropy contribution of each bond is:

```math
S_k=\ln(m_k)
```

The final entropy proxy is:

```math
F_8=\sum_k\ln(m_k)
```

Implementation:

```python
m_k = 0

for i in range(n_points):

    if (
        energy_profile[i] < energy_profile[i - 1]
        and
        energy_profile[i]
        < energy_profile[(i + 1) % n_points]
    ):
        m_k += 1

m_k = max(1, m_k)

entropy_F8 += math.log(m_k)
```

### Physical interpretation

A torsional degree of freedom containing multiple energy minima can access more distinct conformational states than one containing only a single minimum.

The logarithm of the number of minima is therefore used as a simple conformational entropy proxy.

---

# 6. Complete 8D Representation

The final amino-acid representation is:

```math
\boxed{
\mathbf{x}
=
\left[
Q_{\mathrm{balanced}},
\operatorname{Var}(q),
\sum_i\sigma_i^3,
\sum_{i\in\mathrm{heavy}}\epsilon_i,
\sum_{H_{\mathrm{polar}}}q_H,
\sum_{i\in\{N,O\}}\lvert q_i\rvert,
\sum_k\frac{\Delta E_k}{d_k},
\sum_k\ln(m_k)
\right]
}
```

Equivalently:

```text
F1  Intrinsic net charge
F2  Charge variance
F3  van der Waals volume proxy
F4  Dispersion capacity
F5  Hydrogen-bond donor strength
F6  Hydrogen-bond acceptor strength
F7  Torsional rigidity
F8  Conformational entropy proxy
```

---

# 7. Requirements

The implementation requires Python and the following packages:

```text
ParmEd
NumPy
```

Install them with:

```bash
pip install parmed numpy
```

The following standard Python modules are also used:

```text
os
math
json
```

---

# 8. Input Files

The main input is an Amber topology file:

```text
*.prmtop
```

An optional coordinate file can also be supplied:

```text
*.pdb
```

The topology is loaded using:

```python
struct = pmd.load_file(
    prmtop_path,
    xyz=pdb_path
)
```

The residue used for feature extraction is the first residue whose residue name is not one of:

```python
['ACE', 'NME', 'NMA', 'CH3']
```

---

# 9. Filename Convention

The topology filename determines how terminal balancing is performed.

Recommended naming patterns are:

```text
XXX_caps.prmtop
XXX_Ncap.prmtop
XXX_Ccap.prmtop
```

Here, `XXX` identifies the amino acid or modified amino acid.

The interpretation used by the code is:

```text
XXX_caps.prmtop
    -> Internal residue

XXX_Ncap.prmtop
    -> ACE-capped N terminus
    -> Free C-terminal residue

XXX_Ccap.prmtop
    -> NME-capped C terminus
    -> Free N-terminal residue
```

If none of these patterns is detected, the topology is assigned:

```text
Unknown
```

with:

```text
charge_bias = 0
boundary_atoms = []
```

---

# 10. Basic Usage

Import the feature extraction function:

```python
from prmtop_functional import extract_amber_features
```

Extract features from a topology:

```python
features = extract_amber_features(
    "ALA_caps.prmtop"
)
```

If a coordinate file is available:

```python
features = extract_amber_features(
    "ALA_caps.prmtop",
    "ALA_caps.pdb"
)
```

The returned object is a Python dictionary containing metadata and all eight features.

---

# 11. Saving the Results

The module provides a JSON output function:

```python
from prmtop_functional import save_features_to_json
```

Example:

```python
features = extract_amber_features(
    "ALA_caps.prmtop",
    "ALA_caps.pdb"
)

save_features_to_json(
    features,
    "ALA_features.json"
)
```

The generated JSON file contains both feature values and human-readable descriptions.

---

# 12. Output Structure

The output dictionary follows four major physicochemical categories:

```json
{
    "Metadata": {
        "residue_name": "...",
        "topological_state": "...",
        "masked_boundary_atoms": []
    },

    "Electrostatics": {
        "net_charge": {
            "raw_value": 0.0,
            "terminal_bias": 0.0,
            "value": 0.0,
            "desc": "Intrinsic net charge. Explicitly balanced by subtracting terminal bias."
        },

        "charge_variance": {
            "value": 0.0,
            "desc": "Variance of partial charges, with boundary atoms masked to reflect intrinsic polarity."
        }
    },

    "Sterics_And_Dispersion": {
        "vdw_volume_proxy": {
            "value": 0.0,
            "desc": "Sum of LJ sigma cubed."
        },

        "dispersion_capacity": {
            "value": 0.0,
            "desc": "Sum of heavy-atom LJ epsilons."
        }
    },

    "H_Bond_Adaptability": {
        "donor_strength": {
            "value": 0.0,
            "desc": "Sum of partial charges of polar H. N-terminal protons are explicitly masked."
        },

        "acceptor_strength": {
            "value": 0.0,
            "desc": "Sum of absolute charges of atoms with lone pairs. C-terminal OXT is explicitly masked."
        }
    },

    "Flexibility_And_Entropy": {
        "torsional_rigidity": {
            "value": 0.0,
            "desc": "Sum of real dihedral barriers weighted by 1/topological_depth."
        },

        "entropy_proxy": {
            "value": 0.0,
            "desc": "Sum of ln(m_k) over rotatable bonds."
        }
    }
}
```

---

# 13. Feature Summary

To maximize compatibility with GitHub and VS Code, mathematical definitions in this table are written as plain code text rather than LaTeX.

| Feature | Mathematical definition | Amber information used |
| --- | --- | --- |
| **F1 Net Charge** | `Q_raw - Q_bias` | Atomic partial charges |
| **F2 Charge Variance** | `Var(q_i)` | Atomic partial charges |
| **F3 vdW Volume Proxy** | `sum_i sigma_i^3` | LJ sigma |
| **F4 Dispersion Capacity** | `sum_heavy epsilon_i` | LJ epsilon |
| **F5 Donor Strength** | `sum q_(H,polar)` | Atomic charges and bonding topology |
| **F6 Acceptor Strength** | `sum_(N,O) abs(q_i)` | Atomic charges and atom identity |
| **F7 Torsional Rigidity** | `sum_k DeltaE_k / d_k` | Amber dihedral parameters and topology |
| **F8 Entropy Proxy** | `sum_k ln(m_k)` | Torsional energy landscapes |

---

# 14. Physical Organization of the Encoding

The eight dimensions can be conceptually grouped as:

```text
                    8D Physicochemical Encoding
                              |
        +---------------------+---------------------+
        |                     |                     |
  Electrostatics       Sterics/Dispersion      H-bonding
        |                     |                     |
   +----+----+           +----+----+           +----+----+
   F1        F2          F3        F4          F5        F6
 Charge    Charge       vdW     Dispersion    Donor    Acceptor
          Variance     Volume    Capacity     Strength  Strength

                              |
                    Flexibility / Entropy
                              |
                         +----+----+
                         F7        F8
                      Torsional   Entropy
                       Rigidity    Proxy
```

Thus, the representation combines information from four different aspects of molecular mechanics:

```text
Electrostatic interactions
        +
Steric and dispersion interactions
        +
Hydrogen-bond adaptability
        +
Conformational flexibility
```

into a compact eight-dimensional numerical representation.

---

# 15. Source Code

The core API is:

```python
extract_amber_features(
    prmtop_path,
    pdb_path=None
)
```

It returns the full feature dictionary.

JSON output is provided by:

```python
save_features_to_json(
    feature_dict,
    output_path
)
```

---

## License

Add the appropriate license for your project here.

For example:

```text
MIT License
```

---

## Citation

If this encoding scheme is used in published research, please cite the corresponding repository, software release, or publication associated with this implementation.

## FREE PTM AA Database
https://drive.google.com/drive/folders/1rFo6Xxff8AzS9buVBT7cWi0rctDpzXAL?usp=sharing
