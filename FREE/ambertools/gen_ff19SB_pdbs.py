import os
import subprocess

# 1. 配置参数与力场定义
FORCEFIELD_MAIN = "leaprc.protein.ff19SB"
FORCEFIELD_PHOS = "leaprc.phosaa19SB"
OUTPUT_DIR = "amino_acid_pdbs"

# 2. 标准氨基酸及其质子化变体 (Amber残基名, 1字母, 3字母CCD, 描述)
STANDARD_AAS = [
    ("ALA", "A", "ALA", "Neutral"),
    ("ARG", "R", "ARG", "Protonated (+1)"),
    ("ASN", "N", "ASN", "Neutral"),
    ("ASP", "D", "ASP", "Deprotonated (-1)"),
    ("ASH", "D", "ASP", "Protonated (Neutral)"),
    ("CYS", "C", "CYS", "Protonated (Neutral)"),
    ("CYX", "C", "CYS", "Deprotonated/SS-bond capable"),
    ("CYM", "C", "CYS", "Deprotonated (-1)"),
    ("GLN", "Q", "GLN", "Neutral"),
    ("GLU", "E", "GLU", "Deprotonated (-1)"),
    ("GLH", "E", "GLU", "Protonated (Neutral)"),
    ("GLY", "G", "GLY", "Neutral"),
    ("HID", "H", "HIS", "Delta-protonated (Neutral)"),
    ("HIE", "H", "HIS", "Epsilon-protonated (Neutral)"),
    ("HIP", "H", "HIS", "Doubly-protonated (+1)"),
    ("ILE", "I", "ILE", "Neutral"),
    ("LEU", "L", "LEU", "Neutral"),
    ("LYS", "K", "LYS", "Protonated (+1)"),
    ("LYN", "K", "LYS", "Neutral (Deprotonated)"),
    ("MET", "M", "MET", "Neutral"),
    ("PHE", "F", "PHE", "Neutral"),
    ("PRO", "P", "PRO", "Neutral"),
    ("SER", "S", "SER", "Neutral"),
    ("THR", "T", "THR", "Neutral"),
    ("TRP", "W", "TRP", "Neutral"),
    ("TYR", "Y", "TYR", "Neutral"),
    ("VAL", "V", "VAL", "Neutral"),
]

# 3. phosaa19SB 官方支持的 6 种磷酸化残基完整列表
PHOS_AAS = [
    ("SEP", "S", "SEP", "Phosphoserine (-2)"),
    ("S1P", "S", "SEP", "Phosphoserine (-1, protonated)"),
    ("TPO", "T", "TPO", "Phosphothreonine (-2)"),
    ("T1P", "T", "TPO", "Phosphothreonine (-1, protonated)"),
    ("PTR", "Y", "PTR", "Phosphotyrosine (-2)"),
    ("Y1P", "Y", "PTR", "Phosphotyrosine (-1, protonated)"),
]

CAP_TYPES = [
    ("caps", "N-ACE / C-NME", "Double capped"),
    ("Ncap", "N-ACE / C-Free", "N-terminal capped only"),
    ("Ccap", "N-Free / C-NME", "C-terminal capped only"),
]

def build_sequence(res_name, cap_type):
    """根据盖帽类型构建 sequence 字符串"""
    if cap_type == "caps":
        return f"{{ ACE {res_name} NME }}"
    elif cap_type == "Ncap":
        return f"{{ ACE {res_name} }}"
    elif cap_type == "Ccap":
        return f"{{ {res_name} NME }}"
    return ""

def test_residue_in_leap(res_name):
    """测试当前残基是否被当前力场支持"""
    test_script = f"""
source {FORCEFIELD_MAIN}
source {FORCEFIELD_PHOS}
mol = sequence {{ {res_name} }}
quit
"""
    process = subprocess.Popen(
        ["tleap", "-f", "-"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    stdout, _ = process.communicate(input=test_script)
    return "Unknown residue" not in stdout and "Could not find residue" not in stdout

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_records = []

    all_amino_acids = [(aa, "Standard") for aa in STANDARD_AAS] + [(aa, "Phosphorylated") for aa in PHOS_AAS]

    print("=== 开始进行 tleap 自动生成 ===")

    for (res_name, one_letter, ccd_code, protonation), category in all_amino_acids:
        # 验证残基合法性
        if not test_residue_in_leap(res_name):
            print(f"[Skipped] 残基 {res_name} 不存在于当前的 {FORCEFIELD_PHOS} 力场中，已跳过。")
            continue

        for cap_suffix, cap_desc, cap_detail in CAP_TYPES:
            filename = f"{res_name}_{cap_suffix}.pdb"
            filepath = os.path.join(OUTPUT_DIR, filename)
            seq = build_sequence(res_name, cap_suffix)

            tleap_script = f"""
source {FORCEFIELD_MAIN}
source {FORCEFIELD_PHOS}
mol = sequence {seq}
savepdb mol {filepath}
quit
"""
            process = subprocess.Popen(
                ["tleap", "-f", "-"],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            stdout, stderr = process.communicate(input=tleap_script)

            if process.returncode != 0 or not os.path.exists(filepath):
                print(f"[Error] 生成失败: {filename}")
                continue

            print(f"[Success] 成功生成: {filename}")

            all_records.append({
                "filename": filename,
                "res_name": res_name,
                "category": category,
                "one_letter": one_letter,
                "ccd_code": ccd_code,
                "cap_suffix": cap_suffix,
                "cap_desc": cap_desc,
                "protonation": protonation
            })

    # 生成 README.md 说明文档
    readme_path = os.path.join(OUTPUT_DIR, "README.md")
    with open(readme_path, "w", encoding="utf-8") as f:
        f.write("# 氨基酸 PDB 结构库说明文档 (README.md)\n\n")
        f.write(f"本数据库由 Amber `tleap` 自动生成，所用力场：`{FORCEFIELD_MAIN}` 与 `{FORCEFIELD_PHOS}`。\n\n")
        f.write("| 文件名 (File Name) | 类别 (Category) | Amber残基名 | 单字母 (1-Letter) | CCD三字母 (3-Letter) | 封端状态 (Cap Status) | 质子化/磷酸化状态描述 |\n")
        f.write("|---|---|---|---|---|---|---|\n")

        for rec in all_records:
            f.write(f"| `{rec['filename']}` | {rec['category']} | `{rec['res_name']}` | {rec['one_letter']} | {rec['ccd_code']} | {rec['cap_desc']} | {rec['protonation']} |\n")

    print(f"\n全部生成完毕！说明文档存放在: {readme_path}")

if __name__ == "__main__":
    main()