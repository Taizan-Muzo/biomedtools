"""
使用示例
"""

from biomedtools import Molecule, Protein, Pocket
from biomedtools import molecule_fingerprint_similarity
from biomedtools import viz, io

# ========== 分子分析示例 ==========

def molecule_example():
    print("=" * 50)
    print("分子分析示例")
    print("=" * 50)
    
    # 创建分子（阿司匹林）
    aspirin = Molecule.from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
    print(f"分子: {aspirin}")
    print(f"原子数: {aspirin.get_num_atoms()}")
    
    # 计算性质
    print(f"\n性质计算:")
    print(f"  QED (药物相似性): {aspirin.calc_qed():.3f}")
    print(f"  SA (合成可及性): {aspirin.calc_sa():.3f}")
    print(f"  LogP (脂溶性): {aspirin.calc_logp():.3f}")
    print(f"  Lipinski (规则满足): {aspirin.calc_lipinski()}/5")
    
    # 保存
    # aspirin.save_sdf("./aspirin.sdf")
    # aspirin.save_binary("./aspirin.pkl")


# ========== 分子相似度示例 ==========

def similarity_example():
    print("\n" + "=" * 50)
    print("分子相似度示例")
    print("=" * 50)
    
    ethanol = Molecule.from_smiles("CCO")
    propanol = Molecule.from_smiles("CCCO")
    benzene = Molecule.from_smiles("c1ccccc1")
    
    sim1 = molecule_fingerprint_similarity(ethanol, propanol)
    sim2 = molecule_fingerprint_similarity(ethanol, benzene)
    
    print(f"乙醇 vs 丙醇 相似度: {sim1:.3f}")
    print(f"乙醇 vs 苯 相似度: {sim2:.3f}")
    print(f"结论: 乙醇和丙醇更相似（都是醇类）")


# ========== 蛋白质分析示例 ==========

def protein_example():
    print("\n" + "=" * 50)
    print("蛋白质分析示例")
    print("=" * 50)
    
    # KRAS蛋白片段
    seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQRVEDAFYTLVREIRQYRLKKISKEEKTPGCVKIKKCIIM"
    
    kras = Protein.from_fasta(seq, name="KRAS")
    print(f"蛋白质: {kras.name}")
    print(f"序列长度: {len(kras)}")
    print(f"序列前50个氨基酸: {kras.sequence[:50]}...")
    
    # 保存
    # kras.save_fasta("./kras.fasta")
    # kras.save_pdb("./kras.pdb")  # 需要3D结构


# ========== 批量处理示例 ==========

def batch_example():
    print("\n" + "=" * 50)
    print("批量处理示例")
    print("=" * 50)
    
    # 创建分子列表
    smiles_list = [
        "CCO",      # 乙醇
        "CCCO",     # 丙醇
        "CCCCO",    # 丁醇
        "c1ccccc1", # 苯
        "CC(=O)O",  # 乙酸
    ]
    
    molecules = [Molecule.from_smiles(s) for s in smiles_list]
    
    print(f"处理了 {len(molecules)} 个分子:\n")
    print(f"{'SMILES':<15} {'QED':<8} {'LogP':<8} {'Lipinski':<10}")
    print("-" * 45)
    
    for mol in molecules:
        print(f"{str(mol):<15} {mol.calc_qed():<8.3f} {mol.calc_logp():<8.3f} {mol.calc_lipinski():<10}")


if __name__ == "__main__":
    molecule_example()
    similarity_example()
    protein_example()
    batch_example()
    
    print("\n" + "=" * 50)
    print("示例完成!")
    print("=" * 50)
