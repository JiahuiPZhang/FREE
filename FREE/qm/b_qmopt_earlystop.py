import os
import glob
import sys
import traceback
import types
from ase.io import read, write
from pyscf import gto, dft
from pyscf.geomopt.geometric_solver import optimize

# 自定义异常，用于跳出优化循环
class EnergyPlateauException(Exception):
    pass

def run_optimization_with_plateau_check(mf, maxsteps=100, plateau_steps=10, plateau_thresh=1.5e-5):
    """
    带能量平缓检测的几何优化器。
    通过拦截 PySCF 的 Scanner 实现每次迭代的能量记录。
    """
    # 将 DFT 对象转换为梯度计算扫描器
    scanner = mf.nuc_grad_method().as_scanner()
    energy_history = []
    
    # 备份原始的 __call__ 方法
    original_call = scanner.__call__
    
    def custom_call(self, mol):
        # 调用原始计算器获取能量 (e) 和梯度 (g)
        e, g = original_call(mol)
        energy_history.append(e)
        
        # 将当前分子结构备份到 scanner 实例中
        self._last_mol = mol.copy()
        
        # 检查是否满足平缓条件
        if len(energy_history) >= plateau_steps:
            recent_e = energy_history[-plateau_steps:]
            e_diff = max(recent_e) - min(recent_e)
            
            if e_diff < plateau_thresh:
                # 抛出异常以打断 geomeTRIC 的内部循环
                raise EnergyPlateauException(
                    f"连续 {plateau_steps} 步的能量极差 ({e_diff:.2e} a.u.) 小于阈值 ({plateau_thresh} a.u.)"
                )
        return e, g

    # 动态绑定自定义的计算方法
    scanner.__call__ = types.MethodType(custom_call, scanner)
    
    try:
        # 尝试进行正常的几何优化
        mol_eq = optimize(scanner, maxsteps=maxsteps)
        return mol_eq, "NORMAL_CONVERGED"
        
    except EnergyPlateauException as e:
        # 如果触发了提前结束条件，捕获异常并提取最后一步的合理结构
        print(f"\n      [警告/提示] 触发自定义结束标准: {str(e)}")
        return scanner._last_mol, "PLATEAU_REACHED"


def main():
    out_dir = 'optimized'
    os.makedirs(out_dir, exist_ok=True)
    
    log_file = open('batch_opt_fallback.log', 'w')
    err_file = open('batch_opt_fallback.err', 'w')
    
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    
    sys.stdout = log_file
    sys.stderr = err_file
    
    try:
        pdb_files = glob.glob('*.pdb')
        if not pdb_files:
            print("当前目录下未找到 .pdb 文件！", file=original_stdout)
            return
            
        for pdb_file in pdb_files:
            filename = os.path.basename(pdb_file)
            base_name, _ = os.path.splitext(filename)
            
            if base_name.endswith('_caps'):
                charge = 0
            elif base_name.endswith('_N'):
                charge = 1
            elif base_name.endswith('C'):
                charge = -1
            else:
                continue
                
            print(f"\n=========================================")
            print(f"开始处理: {filename} (设定电荷: {charge})")
            
            try:
                mol_ase = read(pdb_file)
                pyscf_atoms = [[atom.symbol, atom.position] for atom in mol_ase]
                
                mol = gto.M(atom=pyscf_atoms, basis='def2-SVP', charge=charge)
                mol.spin = mol.nelectron % 2
                mol.build()
                
                if mol.spin == 0:
                    mf = dft.RKS(mol)
                else:
                    mf = dft.UKS(mol)
                mf.xc = 'B3LYP'
                
                # --- 使用自定义的优化器 ---
                mol_eq, status = run_optimization_with_plateau_check(mf, maxsteps=100)
                
                # 同步坐标并输出
                opt_coords = mol_eq.atom_coords() * 0.529177249
                mol_ase.set_positions(opt_coords)
                
                out_path = os.path.join(out_dir, f"{base_name}.pdb")
                write(out_path, mol_ase)
                
                print(f"---> 优化状态: {status}")
                print(f"---> 成功保存至: {out_path}")
                
                # 实时刷新日志，防止阻塞
                log_file.flush()
                
            except Exception as e:
                print(f"\n---> 失败: {filename} 优化出错！", file=original_stdout)
                print(f"Error in {filename}:", file=sys.stderr)
                traceback.print_exc(file=sys.stderr)
                err_file.flush()
                
    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_file.close()
        err_file.close()
        print("备用批量计算完成。")

if __name__ == '__main__':
    main()