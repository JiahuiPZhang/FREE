import os
import glob
import sys
import traceback
from ase.io import read, write
from pyscf import gto, dft
from pyscf.geomopt.geometric_solver import optimize

def main():
    # 1. 创建输出文件夹
    out_dir = 'optimized'
    os.makedirs(out_dir, exist_ok=True)
    
    # 2. 设置日志重定向
    log_file = open('batch_opt.log', 'w')
    err_file = open('batch_opt.err', 'w')
    
    # 备份原始的 stdout 和 stderr
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    
    # 重定向标准输出和错误到文件
    sys.stdout = log_file
    sys.stderr = err_file
    
    try:
        # 3. 获取当前目录下所有的 pdb 文件
        pdb_files = glob.glob('*.pdb')
        
        if not pdb_files:
            print("当前目录下未找到 .pdb 文件！", file=original_stdout)
            return
            
        for pdb_file in pdb_files:
            filename = os.path.basename(pdb_file)
            base_name, _ = os.path.splitext(filename)
            
            # 4. 根据文件名规则判定电荷
            if base_name.endswith('_caps'):
                charge = 0
            elif base_name.endswith('_N'):
                charge = 1
            elif base_name.endswith('C'): # 匹配 C 结尾（包括 _C）
                charge = -1
            else:
                print(f"[{filename}] 跳过：文件名不符合设定规则。")
                continue
                
            print(f"\n=========================================")
            print(f"开始处理: {filename} (设定电荷: {charge})")
            print(f"=========================================")
            
            try:
                # 5. 读取并初始化分子
                mol_ase = read(pdb_file)
                pyscf_atoms = [[atom.symbol, atom.position] for atom in mol_ase]
                
                mol = gto.M(
                    atom=pyscf_atoms,
                    basis='def2-SVP',
                    charge=charge
                )
                
                # 【极其重要】：自动推断自旋多重度
                # 对于氨基酸残基的质子化/去质子化态，得失质子(H+)不改变电子数
                # 这保证了无论是哪种终端状态，自旋都能被正确设定
                mol.spin = mol.nelectron % 2
                mol.build()
                
                print(f"体系电子数: {mol.nelectron}, 自旋(Spin): {mol.spin}")
                
                # 6. DFT 设置与优化
                if mol.spin == 0:
                    mf = dft.RKS(mol)
                else:
                    mf = dft.UKS(mol)
                
                mf.xc = 'B3LYP'
                
                # 调用 geomeTRIC 优化
                mol_eq = optimize(mf)
                
                # 7. 导出结构
                # 提取坐标并导回 ASE (单位转换：Bohr -> Angstrom)
                opt_coords = mol_eq.atom_coords() * 0.529177
                mol_ase.set_positions(opt_coords)
                
                # 导出 pdb。保留为 pdb 格式对 antechamber 识别残基和原子名最友好
                out_path = os.path.join(out_dir, f"{base_name}.pdb")
                write(out_path, mol_ase)
                
                print(f"\n---> 成功: 结构已保存至 {out_path}")
                
            except Exception as e:
                # 捕获单个分子的报错，防止中断后续分子的计算
                print(f"\n---> 失败: {filename} 优化出错！详情见 .err 文件", file=original_stdout)
                print(f"Error in {filename}:", file=sys.stderr)
                traceback.print_exc(file=sys.stderr)
                
    finally:
        # 8. 恢复控制台输出并关闭文件
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_file.close()
        err_file.close()
        print("批量计算脚本运行结束。日志已分别写入 batch_opt.log 和 batch_opt.err")

if __name__ == '__main__':
    main()
