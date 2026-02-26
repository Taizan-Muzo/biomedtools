"""
ポケット（結合サイト）分析モジュール / Pocket (Binding Site) Analysis Module

タンパク質ポケット定義・結合サイト分析ライブラリ
Protein pocket definition and binding site analysis library
"""

from typing import List, Optional
from typing_extensions import Self

from datetime import datetime
import numpy as np
import os
import pickle
import re

try:
    from scipy import spatial
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    spatial = None


class Pocket:
    """
    Protein pocket (binding site) / タンパク質ポケット（結合サイト）
    
    Represents a binding pocket within a protein structure.
    タンパク質構造内の結合ポケットを表現。
    """
    
    def __init__(self, name: Optional[str] = None) -> None:
        self.atoms = []
        self.conformer = None
        self.orig_indices = []
        self.orig_protein = None
        self.name = name

    @classmethod
    def from_protein_subseq(cls, protein, indices: List[int]) -> Self:
        """
        Create pocket from protein residue indices.
        タンパク質の残基インデックスからポケットを作成。
        """
        from .protein import Protein
        if not isinstance(protein, Protein):
            raise TypeError("protein must be a Protein instance")
        
        pocket = cls()
        pocket.orig_indices = indices
        pocket.orig_protein = protein
        
        for i in indices:
            if i < 0 or i >= len(protein.residues):
                continue
            residue = protein.residues[i]
            for atom_i in residue.atoms:
                if atom_i < len(protein.all_atom):
                    atom = protein.all_atom[atom_i]
                    pocket.atoms.append(atom)
                    if "pos" in atom:
                        if pocket.conformer is None:
                            pocket.conformer = []
                        pocket.conformer.append(atom["pos"])
        
        if pocket.conformer:
            pocket.conformer = np.array(pocket.conformer)
        
        return pocket

    @classmethod
    def from_protein_ref_ligand(cls, protein, ligand, radius: float = 10.0) -> Self:
        """
        Define pocket based on reference ligand.
        参照リガンドに基づいてポケットを定義。
        
        Selects residues within specified radius of ligand atoms.
        リガンド原子の指定半径内の残基を選択。
        """
        from .molecule import Molecule
        from .protein import Protein
        
        if not isinstance(protein, Protein):
            raise TypeError("protein must be a Protein instance")
        if not isinstance(ligand, Molecule):
            raise TypeError("ligand must be a Molecule instance")
        
        if ligand.conformer is None:
            raise ValueError("Ligand must have 3D conformer")
        
        sel_idx = []
        for pos in ligand.conformer:
            for i, residue in enumerate(protein.residues):
                if hasattr(residue, 'center_of_mass'):
                    dist = np.linalg.norm(residue.center_of_mass - pos)
                    if dist < radius and i not in sel_idx:
                        sel_idx.append(i)
        
        return cls.from_protein_subseq(protein, sel_idx)

    @classmethod
    def from_binary_file(cls, file: str) -> Self:
        """
        Load pocket from binary file.
        バイナリファイルからポケットを読み込む。
        """
        return pickle.load(open(file, "rb"))

    def _add_name(self) -> None:
        """
        Auto-generate pocket name.
        ポケット名を自動生成。
        """
        if self.name is None:
            self.name = "pocket_" + re.sub(r"[-:.]", "_",
                                           datetime.now().isoformat(sep="_", timespec="milliseconds"))

    def get_num_atoms(self) -> int:
        """
        Get number of atoms in pocket.
        ポケット内の原子数を取得。
        """
        return len(self.atoms)

    def get_center(self) -> Optional[np.ndarray]:
        """
        Get pocket center coordinates.
        ポケット中心座標を取得。
        """
        if self.conformer is not None and len(self.conformer) > 0:
            return np.mean(self.conformer, axis=0)
        return None

    def get_radius(self) -> float:
        """
        Get pocket radius (max distance from center).
        ポケット半径（中心からの最大距離）を取得。
        """
        center = self.get_center()
        if center is None or self.conformer is None:
            return 0.0
        distances = np.linalg.norm(self.conformer - center, axis=1)
        return float(np.max(distances))

    def save_binary(self, file: Optional[str] = None, overwrite: bool = False) -> str:
        """
        Save as binary file.
        バイナリファイルとして保存。
        """
        self._add_name()
        if file is None:
            file = f"./{self.name}.pkl"
        if not os.path.exists(file) or overwrite:
            with open(file, "wb") as f:
                pickle.dump(self, f)
        return file

    def save_pdb(self, file: Optional[str] = None, overwrite: bool = False) -> str:
        """
        Save as PDB file (via original protein).
        PDBファイルとして保存（元のタンパク質経由）。
        """
        if self.orig_protein is None:
            raise ValueError("Cannot save PDB without original protein")
        self._add_name()
        if file is None:
            file = f"./{self.name}.pdb"
        return self.orig_protein.save_pdb(file, overwrite, self.orig_indices)

    def __str__(self) -> str:
        if self.orig_protein is not None:
            indices_str = ",".join([str(idx+1) for idx in self.orig_indices[:5]])
            if len(self.orig_indices) > 5:
                indices_str += "..."
            return f"Pocket({self.get_num_atoms()} atoms, residues {indices_str})"
        return f"Pocket({self.get_num_atoms()} atoms)"

    def __repr__(self) -> str:
        return self.__str__()


def estimate_ligand_atom_num(pocket: Pocket) -> int:
    """
    Estimate suitable ligand atom count for this pocket.
    このポケットに適したリガンド原子数を推定。
    
    Based on pocket spatial size using simple heuristic.
    ポケット空間サイズに基づく簡易ヒューリスティック。
    """
    if not SCIPY_AVAILABLE:
        # Simple estimate based on atom count / 原子数に基づく簡易推定
        return min(pocket.get_num_atoms() // 3, 50)
    
    if pocket.conformer is None or len(pocket.conformer) < 2:
        return 10
    
    # Calculate distances between pocket atoms / ポケット原子間の距離を計算
    try:
        dist = spatial.distance.pdist(pocket.conformer, metric='euclidean')
        if len(dist) == 0:
            return 10
        
        # Use max distance as spatial size estimate / 最大距離を空間サイズの推定値として使用
        space_size = np.max(dist)
        
        # Heuristic: approximately 15-20 atoms per 10 Angstroms
        # ヒューリスティック：約10Åあたり15-20原子
        estimated = int(space_size * 1.5)
        return max(5, min(estimated, 100))  # Limit to 5-100 / 5-100に制限
    except Exception:
        return 10
