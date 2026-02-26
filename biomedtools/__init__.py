"""
BioMedTools - 軽量な生物医药計算ツールライブラリ / Lightweight Biomedical Calculation Tools

機能 / Features:
- 分子処理・プロパティ計算 / Molecule processing & property calculation
- タンパク質シーケンス処理 / Protein sequence processing  
- ポケット分析 / Pocket analysis
- 可視化 / Visualization
- IOツール / IO tools
"""

from .molecule import (
    Molecule,
    molecule_fingerprint_similarity,
    check_identical_molecules,
    calc_mol_diversity,
    calc_mol_rmsd,
)

from .protein import (
    Protein,
    protein_sequence_similarity,
    Residue,
    AA_NAME_SYM,
    AA_NAME_NUMBER,
    BACKBONE_NAMES,
)

from .pocket import (
    Pocket,
    estimate_ligand_atom_num,
)

from . import viz
from . import io

__version__ = "0.1.0"
__all__ = [
    # Molecule / 分子
    "Molecule",
    "molecule_fingerprint_similarity",
    "check_identical_molecules",
    "calc_mol_diversity",
    "calc_mol_rmsd",
    # Protein / タンパク質
    "Protein",
    "protein_sequence_similarity",
    "Residue",
    "AA_NAME_SYM",
    "AA_NAME_NUMBER",
    "BACKBONE_NAMES",
    # Pocket / ポケット
    "Pocket",
    "estimate_ligand_atom_num",
    # Submodules / サブモジュール
    "viz",
    "io",
]
