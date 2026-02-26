# BioMedTools

軽量な生物医药計算ツールライブラリ - ディープラーニング不要

## 特徴

- ✅ **Pure Calculation, No Models** - PyTorch/TensorFlowなどのディープラーニングフレームワークに依存しない
- ✅ **Lightning Fast Startup** - ミリ秒級のインポート、大規模モデルの読み込み不要
- ✅ **Fully Offline** - コア機能はネットワーク不要（オプションのデータベースクエリを除く）
- ✅ **RDKit Powered** - 業界標準のRDKit化学情報学ライブラリをベースに構築

## 機能

| モジュール | 機能 |
|-----------|------|
| `molecule` | 分子処理、プロパティ計算(QED/SA/LogP/Lipinski)、類似度計算 |
| `protein` | タンパク質シーケンス処理、PDB解析、シーケンスアラインメント |
| `pocket` | タンパク質ポケット定義、結合サイト分析 |
| `viz` | 2D/3D可視化 |
| `io` | フォーマット変換(SDF/PDB/PKL/FASTA) |

## インストール

```bash
pip install rdkit numpy scipy biopython
```

## クイックスタート

```python
from biomedtools import Molecule, Protein

# Molecular analysis
mol = Molecule.from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
print(f"QED: {mol.calc_qed():.3f}")      # Drug-likeness
print(f"SA: {mol.calc_sa():.3f}")        # Synthetic accessibility
print(f"LogP: {mol.calc_logp():.3f}")    # Lipophilicity

# Protein analysis
protein = Protein.from_fasta("MTEYKLVVV...")
print(f"Length: {len(protein)}")
protein.save_pdb("./protein.pdb")
```

## API Reference

### Molecule Class

```python
# Create from different formats
mol = Molecule.from_smiles("CCO")                    # From SMILES
mol = Molecule.from_selfies("[C][C][O]")             # From SELFIES
mol = Molecule.from_sdf_file("./molecule.sdf")       # From SDF file

# Property calculation (all offline, RDKit-based)
mol.calc_qed()        # QED drug-likeness score (0-1)
mol.calc_sa()         # SA synthetic accessibility (0-10)
mol.calc_logp()       # LogP partition coefficient
mol.calc_lipinski()   # Lipinski's Rule of Five (0-5)
mol.get_num_atoms()   # Number of atoms

# Save
mol.save_sdf("./output.sdf")      # SDF format
mol.save_binary("./output.pkl")   # Binary format
```

### Protein Class

```python
# Create from different formats
protein = Protein.from_fasta("MTEYK...")             # From FASTA sequence
protein = Protein.from_pdb_file("./protein.pdb")     # From PDB file

# Properties
len(protein)              # Sequence length
protein.sequence          # Amino acid sequence
protein.residues          # Residue list

# Save
protein.save_pdb("./output.pdb")      # PDB format
protein.save_fasta("./output.fasta")  # FASTA format
protein.save_binary("./output.pkl")   # Binary format
```

### Molecular Similarity

```python
from biomedtools import molecule_fingerprint_similarity

mol1 = Molecule.from_smiles("CCO")   # Ethanol
mol2 = Molecule.from_smiles("CCCO")  # Propanol

# Calculate Tanimoto similarity
similarity = molecule_fingerprint_similarity(mol1, mol2, fingerprint_type="morgan")
# Returns: 0.0 to 1.0
```

### Visualization

```python
from biomedtools import viz

# 2D visualization (RDKit only, fully offline)
viz.draw_molecule_2d(mol, "./molecule.png", size=(400, 400))

# Grid visualization for multiple molecules
viz.draw_molecules_grid([mol1, mol2, mol3], "./grid.png", mols_per_row=4)

# 3D visualization (requires PyMOL)
viz.visualize_protein_3d("./protein.pdb", "./protein.png")
```

## Comparison with OpenBioMed

| Feature | OpenBioMed | BioMedTools |
|---------|-----------|-------------|
| Deep Learning Models | Required | Not Required |
| PyTorch | Required | Not Required |
| RDKit | Required | Required |
| Startup Time | Slow (model loading) | Fast |
| Offline Capability | Partial | Full (core features) |
| Use Case | AI-driven drug discovery | Basic cheminformatics analysis |

## Dependencies

### Required
- **RDKit** (>=2022.03.1) - Core cheminformatics functionality
- **NumPy** (>=1.20.0) - Numerical computations
- **SciPy** (>=1.7.0) - Scientific computing

### Optional
- **Biopython** (>=1.79) - Protein sequence alignment
- **PyMOL** (open-source) - 3D visualization
- **imageio** (>=2.9.0) - Image processing

## License

MIT
