import os
import subprocess
import json
import math
import numpy as np
import parmed as pmd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA

# ==========================================================
# 1. 氨基酸分类定义 (基于生理 pH 7.4)
# ==========================================================
AA_CLASSES = {
    'Basic': {'aas': ['ARG', 'LYS', 'HIE'], 'color': 'red', 'label': 'Basic (+)'},
    'Acidic': {'aas': ['ASP', 'GLU'], 'color': 'blue', 'label': 'Acidic (-)'},
    'Polar': {'aas': ['ASN', 'GLN', 'SER', 'THR', 'CYS', 'TYR'], 'color': 'green', 'label': 'Polar Uncharged'},
    'Hydrophobic': {'aas': ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP', 'GLY'], 'color': 'gray', 'label': 'Hydrophobic'}
}

PHOS_DEPROT = ['SEP', 'TPO', 'PTR']  # 非质子化 (-2)
PHOS_PROT = ['S1P', 'T1P', 'Y1P']    # 质子化 (-1)

STD_AAS = [aa for cat in AA_CLASSES.values() for aa in cat['aas']]
ALL_AAS = STD_AAS + PHOS_DEPROT + PHOS_PROT

PDB_DIR = "ff19SB_pdbs"


# ==========================================================
# 2. 自动化 tleap 参数生成
# ==========================================================
def generate_amber_params(res_name):
    pdb_path = os.path.join(PDB_DIR, f"{res_name}_caps.pdb")
    prmtop_path = os.path.join(PDB_DIR, f"{res_name}_caps.prmtop")
    inpcrd_path = os.path.join(PDB_DIR, f"{res_name}_caps.inpcrd")
    
    if not os.path.exists(pdb_path): return None
        
    tleap_script = f"""
source leaprc.protein.ff19SB
source leaprc.phosaa19SB
mol = loadpdb {pdb_path}
saveamberparm mol {prmtop_path} {inpcrd_path}
quit
"""
    process = subprocess.Popen(
        ['tleap', '-f', '-'], 
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    process.communicate(tleap_script)
    return prmtop_path if os.path.exists(prmtop_path) else None


# ==========================================================
# 3. 特征提取引擎
# ==========================================================
def extract_8d_features(prmtop_path, pdb_path):
    struct = pmd.load_file(prmtop_path, xyz=pdb_path)
    target_res = next((res for res in struct.residues if res.name not in ['ACE', 'NME', 'NMA', 'CH3']), None)
    if not target_res: return None
    
    net_charge = sum(a.charge for a in target_res.atoms)
    charge_var = np.var([a.charge for a in target_res.atoms])
    vdw_vol = sum(a.sigma**3 for a in target_res.atoms)
    dispersion = sum(a.epsilon for a in target_res.atoms if a.atomic_number > 1)
    
    donor_str = 0.0
    for a in target_res.atoms:
        if a.atomic_number == 1:
            for bond in a.bonds:
                partner = bond.atom1 if bond.atom2 == a else bond.atom2
                if partner.atomic_number in [7, 8, 16]:
                    donor_str += a.charge
                    break
                    
    acceptor_str = sum(abs(a.charge) for a in target_res.atoms if a.atomic_number in [7, 8])
    
    ca_atom = next((a for a in target_res.atoms if a.name == 'CA'), None)
    if not ca_atom: return None
    
    depths = {ca_atom: 0}
    queue = [(ca_atom, 0)]
    visited = {ca_atom}
    while queue:
        curr, d = queue.pop(0)
        for bond in curr.bonds:
            partner = bond.atom1 if bond.atom2 == curr else bond.atom2
            if partner in target_res.atoms and partner not in visited:
                visited.add(partner)
                depths[partner] = d + 1
                queue.append((partner, d + 1))

    rigidity_f7, entropy_f8 = 0.0, 0.0
    processed_bonds = set()
    for bond in struct.bonds:
        a1, a2 = bond.atom1, bond.atom2
        if a1 in target_res.atoms and a2 in target_res.atoms and a1.atomic_number > 1 and a2.atomic_number > 1:
            if a1.name in {'N', 'CA', 'C', 'O', 'H'} and a2.name in {'N', 'CA', 'C', 'O', 'H'}: continue
                
            bond_pair = frozenset([a1.idx, a2.idx])
            if bond_pair in processed_bonds: continue
            processed_bonds.add(bond_pair)

            bond_depth = max(depths.get(a1, 1), depths.get(a2, 1))
            angles = np.linspace(0, 2 * np.pi, 360, endpoint=False)
            energy_prof = np.zeros_like(angles)
            has_dih = False

            for dih in struct.dihedrals:
                if dih.type and ((dih.atom2 == a1 and dih.atom3 == a2) or (dih.atom2 == a2 and dih.atom3 == a1)):
                    has_dih = True
                    energy_prof += dih.type.phi_k * (1 + np.cos(dih.type.per * angles - dih.type.phase))

            if has_dih:
                rigidity_f7 += (1.0 / bond_depth) * (np.max(energy_prof) - np.min(energy_prof))
                m_k = max(1, sum(1 for i in range(len(energy_prof)) if energy_prof[i] < energy_prof[i - 1] and energy_prof[i] < energy_prof[(i + 1) % len(energy_prof)]))
                entropy_f8 += math.log(m_k)

    return [net_charge, charge_var, vdw_vol, dispersion, donor_str, acceptor_str, rigidity_f7, entropy_f8]


# ==========================================================
# 4. 置信椭圆绘制辅助函数
# ==========================================================
def add_confidence_ellipse(ax, x, y, color, chi2_val=4.605):
    if len(x) < 2: return
        
    cov = np.cov(x, y)
    if np.isnan(cov).any() or np.isinf(cov).any(): return
        
    eigvals, eigvecs = np.linalg.eigh(cov)
    eigvals = np.maximum(eigvals, 0) 
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    
    angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
    width, height = 2 * np.sqrt(chi2_val * eigvals)
    
    ellip = Ellipse(xy=(np.mean(x), np.mean(y)), width=width, height=height, 
                    angle=angle, facecolor='none', edgecolor=color, 
                    linestyle=(0, (4, 3)), linewidth=1.5, zorder=1)
    
    ax.add_patch(ellip)


# ==========================================================
# 主流程
# ==========================================================
def main():
    print("1. 开始提取 8D 向量...")
    raw_data = {}
    for aa in ALL_AAS:
        prmtop = generate_amber_params(aa)
        if prmtop:
            vec = extract_8d_features(prmtop, prmtop.replace('.prmtop', '.pdb'))
            if vec: raw_data[aa] = vec

    print("2. 基于【全体氨基酸 (含磷酸化)】执行 Z-score 归一化...")
    # 核心修改点：这里的基准面扩大到了 ALL_AAS
    all_matrix = np.array([raw_data[aa] for aa in ALL_AAS if aa in raw_data])
    mean = np.mean(all_matrix, axis=0)
    std = np.std(all_matrix, axis=0)
    std[std == 0] = 1e-8 

    norm_data = {aa: (np.array(vec) - mean) / std for aa, vec in raw_data.items()}
    valid_aas = list(norm_data.keys())
    X = np.array([norm_data[aa] for aa in valid_aas])
    
    print("3. 执行 PCA...")
    pca = PCA(n_components=3)
    X_pca = pca.fit_transform(X)
    ev_ratio = pca.explained_variance_ratio_ * 100

    print("4. 生成可视化图像 (附带 90% 置信椭圆)...")
    fig, axes = plt.subplots(3, 1, figsize=(9, 21))
    pairs = [(0, 1), (0, 2), (1, 2)]

    for ax, (idx_x, idx_y) in zip(axes, pairs):
        # 4a. 为标准氨基酸分类绘制置信椭圆
        for class_name, info in AA_CLASSES.items():
            class_aas = [aa for aa in info['aas'] if aa in valid_aas]
            if len(class_aas) >= 2:
                indices = [valid_aas.index(aa) for aa in class_aas]
                add_confidence_ellipse(
                    ax, 
                    X_pca[indices, idx_x], 
                    X_pca[indices, idx_y], 
                    color=info['color']
                )

        # 4b. 绘制所有散点
        for i, aa in enumerate(valid_aas):
            x, y = X_pca[i, idx_x], X_pca[i, idx_y]
            
            if aa in PHOS_DEPROT:
                color, marker = 'black', '^'
            elif aa in PHOS_PROT:
                color, marker = 'black', 'v'
            else:
                marker = 'o'
                color = next(info['color'] for info in AA_CLASSES.values() if aa in info['aas'])

            ax.scatter(x, y, color=color, marker=marker, s=120, edgecolors='black', linewidth=1, zorder=2)
            ax.annotate(aa, (x, y), xytext=(5, 5), textcoords='offset points', fontsize=10, zorder=3)
            
        ax.set_xlabel(f'PC{idx_x + 1} ({ev_ratio[idx_x]:.1f}%)', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'PC{idx_y + 1} ({ev_ratio[idx_y]:.1f}%)', fontsize=12, fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.6, zorder=0)

    # 4c. 构建图例
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=c['color'], markersize=11, markeredgecolor='k', label=c['label'])
        for c in AA_CLASSES.values()
    ] + [
        Line2D([0], [0], marker='^', color='w', markerfacecolor='black', markersize=11, label='Phosphorylated (-2)'),
        Line2D([0], [0], marker='v', color='w', markerfacecolor='black', markersize=11, label='Phosphorylated (-1)'),
        Line2D([0], [0], color='black', linestyle=(0, (4, 3)), linewidth=1.5, label='90% Confidence Ellipse')
    ]
    
    fig.legend(handles=legend_elements, loc='lower center', ncol=3, bbox_to_anchor=(0.5, 0.02), fontsize=12)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.08, hspace=0.2) 
    
    # 核心修改点：调整了输出文件的名称
    plt.savefig('pca_conf_all_norm.png', dpi=300, bbox_inches='tight')
    print("5. 图像已保存至 pca_conf_all_norm.png")

if __name__ == "__main__":
    main()