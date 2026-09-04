"""
prmtop_functional.py (v1.1)

基于 Amber 拓扑文件 (prmtop) 的氨基酸物理化学特征 8D 提取管线。
包含端基效应平衡机制 (Terminal Effect Balancing)，自动解析拓扑位置，
消除暴漏在 N 端或 C 端的系统性偏差，还原氨基酸本征理化性质。
"""

import os
import math
import parmed as pmd
import numpy as np
import json

# =====================================================================
# 辅助模块：边界状态解析与原子掩码
# =====================================================================

def _parse_terminal_state(filename):
    """
    解析文件名，判断氨基酸在肽链中的拓扑位置，并返回补偿参数。
    返回: (拓扑状态, 电荷偏差基线, 需要屏蔽的边界原子名列表)
    """
    if '_caps' in filename:
        return "Internal", 0.0, []
    elif '_Ncap' in filename:
        # N端被 ACE 封端，C端自由 (拓扑学上的 C-terminal residue)
        # 显式考虑 C 端羧基 (C-terminal carboxyls, COO-) 带来的 -1 系统电荷
        return "C_terminal", -1.0, ['OXT']
    elif '_Ccap' in filename:
        # C端被 NME 封端，N端自由 (拓扑学上的 N-terminal residue)
        # 显式考虑 N 端质子 (N-terminal protons, NH3+) 带来的 +1 系统电荷
        return "N_terminal", 1.0, ['H1', 'H2', 'H3']
    return "Unknown", 0.0, []

# =====================================================================
# 维度一：静电特征 (Electrostatics)
# =====================================================================

def _get_net_charge_balanced(residue, charge_bias):
    """F1: 提取原始净电荷，并显式减去边界偏差以平衡为本征电荷"""
    raw_charge = sum(atom.charge for atom in residue.atoms)
    balanced_charge = raw_charge - charge_bias
    return float(raw_charge), float(balanced_charge)

def _get_charge_variance_balanced(residue, boundary_atoms):
    """F2: 提取电荷方差。跳过极端极化的边界原子，反映残基本征极性"""
    charges = [atom.charge for atom in residue.atoms if atom.name not in boundary_atoms]
    return float(np.var(charges))

# =====================================================================
# 维度二：空间与疏水特征 (Sterics & Dispersion)
# =====================================================================

def _get_vdw_volume_proxy(residue):
    """F3: 近似范德华体积代理"""
    return sum(atom.sigma**3 for atom in residue.atoms)

def _get_dispersion_capacity(residue):
    """F4: 疏水色散力极化率 (仅计重原子)"""
    return sum(atom.epsilon for atom in residue.atoms if atom.atomic_number > 1)

# =====================================================================
# 维度三：氢键适配度 (H-bond Adaptability)
# =====================================================================

def _get_donor_strength_balanced(residue, boundary_atoms):
    """F5: 供体静电强度。屏蔽 N 端的额外质子，平衡供体信号"""
    donor_strength = 0.0
    for atom in residue.atoms:
        if atom.name in boundary_atoms:
            continue
        if atom.atomic_number == 1:
            for bond in atom.bonds:
                partner = bond.atom1 if bond.atom2 == atom else bond.atom2
                if partner.atomic_number in [7, 8, 16]:
                    donor_strength += atom.charge
                    break
    return donor_strength

def _get_acceptor_strength_balanced(residue, boundary_atoms):
    """F6: 受体静电强度。屏蔽 C 端的 OXT 原子，平衡受体信号"""
    acceptor_strength = 0.0
    for atom in residue.atoms:
        if atom.name in boundary_atoms:
            continue
        if atom.atomic_number in [7, 8]:
            acceptor_strength += abs(atom.charge)
    return acceptor_strength

# =====================================================================
# 维度四：柔性与构象熵 (Flexibility & Entropy) - 保持与 v1.0 一致
# =====================================================================

def _get_flexibility_and_entropy(residue, struct):
    ca_atom = next((a for a in residue.atoms if a.name == 'CA'), None)
    if not ca_atom:
        return 0.0, 0.0

    depths = {}
    queue = [(ca_atom, 0)]
    visited = {ca_atom}
    while queue:
        curr, d = queue.pop(0)
        depths[curr] = d
        for partner in curr.bonded_atoms:
            if partner in residue.atoms and partner not in visited:
                visited.add(partner)
                queue.append((partner, d + 1))

    backbone_names = {'N', 'CA', 'C', 'O', 'H', 'HA'}
    rigidity_F7, entropy_F8 = 0.0, 0.0
    processed_bonds = set()
    
    for bond in struct.bonds:
        a1, a2 = bond.atom1, bond.atom2
        if (a1 in residue.atoms) and (a2 in residue.atoms) and \
           (a1.atomic_number > 1) and (a2.atomic_number > 1):
            
            if a1.name in backbone_names and a2.name in backbone_names:
                continue
                
            bond_pair = frozenset([a1.idx, a2.idx])
            if bond_pair in processed_bonds:
                continue
            processed_bonds.add(bond_pair)

            bond_depth = max(depths.get(a1, 1), depths.get(a2, 1))
            angles = np.linspace(0, 2 * np.pi, 360, endpoint=False)
            energy_profile = np.zeros_like(angles)
            has_dihedral = False

            for dih in struct.dihedrals:
                if dih.type is None:
                    continue
                if (dih.atom2 == a1 and dih.atom3 == a2) or (dih.atom2 == a2 and dih.atom3 == a1):
                    has_dihedral = True
                    energy_profile += dih.type.phi_k * (1 + np.cos(dih.type.per * angles - dih.type.phase))

            if has_dihedral:
                delta_E = np.max(energy_profile) - np.min(energy_profile)
                rigidity_F7 += (1.0 / bond_depth) * delta_E

                m_k = 0
                n_points = len(energy_profile)
                for i in range(n_points):
                    if energy_profile[i] < energy_profile[i - 1] and energy_profile[i] < energy_profile[(i + 1) % n_points]:
                        m_k += 1
                
                m_k = max(1, m_k)
                entropy_F8 += math.log(m_k)

    return float(rigidity_F7), float(entropy_F8)


# =====================================================================
# 核心 API 函数
# =====================================================================

def extract_amber_features(prmtop_path, pdb_path=None):
    """统一特征提取入口，附带自动边界平衡功能"""
    struct = pmd.load_file(prmtop_path, xyz=pdb_path)
    filename = os.path.basename(prmtop_path)
    
    # 1. 自动解析拓扑学位置与补偿常数
    term_state, charge_bias, boundary_atoms = _parse_terminal_state(filename)
    
    target_res = None
    for res in struct.residues:
        if res.name not in ['ACE', 'NME', 'NMA', 'CH3']:
            target_res = res
            break
            
    if target_res is None:
        raise ValueError("Cannot find a valid amino acid residue in the topology.")

    # 2. 提取特征，并传入补偿机制
    raw_q, bal_q = _get_net_charge_balanced(target_res, charge_bias)
    f2 = _get_charge_variance_balanced(target_res, boundary_atoms)
    f3 = _get_vdw_volume_proxy(target_res)
    f4 = _get_dispersion_capacity(target_res)
    f5 = _get_donor_strength_balanced(target_res, boundary_atoms)
    f6 = _get_acceptor_strength_balanced(target_res, boundary_atoms)
    f7, f8 = _get_flexibility_and_entropy(target_res, struct)

    # 3. 封装带注释的数据结构 (显式展示平衡过程)
    feature_dict = {
        "Metadata": {
            "residue_name": target_res.name,
            "topological_state": term_state,
            "masked_boundary_atoms": boundary_atoms
        },
        "Electrostatics": {
            "net_charge": {
                "raw_value": raw_q,
                "terminal_bias": charge_bias,
                "value": bal_q,
                "desc": "Intrinsic net charge. Explicitly balanced by subtracting terminal bias."
            },
            "charge_variance": {
                "value": f2,
                "desc": "Variance of partial charges, with boundary atoms masked to reflect intrinsic polarity."
            }
        },
        "Sterics_And_Dispersion": {
            "vdw_volume_proxy": {
                "value": f3,
                "desc": "Sum of LJ sigma cubed."
            },
            "dispersion_capacity": {
                "value": f4,
                "desc": "Sum of heavy-atom LJ epsilons."
            }
        },
        "H_Bond_Adaptability": {
            "donor_strength": {
                "value": f5,
                "desc": "Sum of partial charges of polar H. N-terminal protons are explicitly masked."
            },
            "acceptor_strength": {
                "value": f6,
                "desc": "Sum of absolute charges of atoms with lone pairs. C-terminal OXT is explicitly masked."
            }
        },
        "Flexibility_And_Entropy": {
            "torsional_rigidity": {
                "value": f7,
                "desc": "Sum of real dihedral barriers weighted by 1/topological_depth."
            },
            "entropy_proxy": {
                "value": f8,
                "desc": "Sum of ln(m_k) over rotatable bonds."
            }
        }
    }
    
    return feature_dict


def save_features_to_json(feature_dict, output_path):
    """保存为 JSON"""
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(feature_dict, f, indent=4, ensure_ascii=False)
    print(f"[Success] Saved balanced features to: {output_path}")