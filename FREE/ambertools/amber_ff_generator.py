import os
import glob
import subprocess
import shutil

def main():
    input_dir = 'optimized'
    output_dir = 'ff_output'
    
    # 检查输入文件夹是否存在
    if not os.path.exists(input_dir):
        print(f"未找到输入文件夹: {input_dir}。请确保已经在当前目录下。")
        return

    # 创建主输出文件夹
    os.makedirs(output_dir, exist_ok=True)
    
    # 获取所有的 PDB 文件
    pdb_files = glob.glob(os.path.join(input_dir, '*.pdb'))
    if not pdb_files:
        print(f"在 {input_dir} 中未找到 .pdb 文件！")
        return
        
    print(f"找到 {len(pdb_files)} 个 PDB 文件，准备生成 AMBER 力场文件...")
    
    # 严格顺序处理每个分子和对应的子文件夹，避免目录跳转和中间文件读写发生冲突
    for pdb_path in pdb_files:
        filename = os.path.basename(pdb_path)
        base_name, _ = os.path.splitext(filename)
        
        # 1. 判断电荷
        if base_name.endswith('_caps'):
            charge = 0
        elif base_name.endswith('_N'):
            charge = 1
        elif base_name.endswith('C'): # 匹配 C 结尾
            charge = -1
        else:
            print(f"[{filename}] 跳过：文件名不符合设定规则。")
            continue
            
        print(f"\n=========================================")
        print(f"开始处理: {filename} (电荷: {charge})")
        
        # 2. 为当前分子创建专用的子文件夹
        sub_dir = os.path.join(output_dir, base_name)
        os.makedirs(sub_dir, exist_ok=True)
        
        # 获取绝对路径，方便在子目录中调用
        abs_pdb_path = os.path.abspath(pdb_path)
        abs_sub_dir = os.path.abspath(sub_dir)
        
        # 记录当前工作目录，并切换到子文件夹
        # 这样 antechamber 产生的 sqm.out, ANTECHAMBER_* 等大量临时文件都会被限制在这个子文件夹中
        original_cwd = os.getcwd()
        os.chdir(abs_sub_dir)
        
        try:
            # -----------------------------------------------------
            # 步骤 A: 运行 antechamber 生成带 AM1-BCC 电荷的 mol2 文件
            # -----------------------------------------------------
            mol2_out = f"{base_name}.mol2"
            antechamber_cmd = [
                "antechamber",
                "-i", abs_pdb_path,
                "-fi", "pdb",
                "-o", mol2_out,
                "-fo", "mol2",
                "-c", "bcc",      # 使用 AM1-BCC 电荷
                "-nc", str(charge),
                "-at", "gaff2"    # 显式指定使用 GAFF2 原子类型
            ]
            
            print(f"  -> 正在执行: {' '.join(antechamber_cmd)}")
            subprocess.run(antechamber_cmd, check=True, stdout=subprocess.DEVNULL) # 隐藏标准输出，保持终端整洁
            
            # -----------------------------------------------------
            # 步骤 B: 运行 parmchk2 检查缺失参数并生成 frcmod 文件
            # -----------------------------------------------------
            frcmod_out = f"{base_name}.frcmod"
            parmchk_cmd = [
                "parmchk2",
                "-i", mol2_out,
                "-f", "mol2",
                "-o", frcmod_out,
                "-s", "2"         # -s 2 代表 GAFF2，这是与前面 antechamber 匹配的关键
            ]
            
            print(f"  -> 正在执行: {' '.join(parmchk_cmd)}")
            subprocess.run(parmchk_cmd, check=True, stdout=subprocess.DEVNULL)
            
            print(f"  ---> 成功！所需力场文件已生成于: {sub_dir}")
            
        except subprocess.CalledProcessError as e:
            print(f"  ---> [错误] AMBER 工具执行失败，详情: {e}")
        except FileNotFoundError:
            print(f"  ---> [严重错误] 未找到 antechamber 或 parmchk2！请检查 AMBER 环境变量 (AMBERHOME) 是否已正确配置。")
            os.chdir(original_cwd)
            break
        finally:
            # 无论成功与否，必须切回主目录以便处理下一个分子
            os.chdir(original_cwd)

if __name__ == '__main__':
    main()
