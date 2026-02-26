"""
タンパク質処理モジュール / Protein Processing Module

タンパク質シーケンス処理・PDB解析ライブラリ
Protein sequence processing and PDB parsing library
"""

from typing import Any, Dict, Iterator, List, Optional, Tuple
from typing_extensions import Self

from datetime import datetime
import numpy as np
import os
import pickle
import re

try:
    from rdkit import Chem
    from rdkit.Chem import PeriodicTable
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    # Simple periodic table fallback / 簡易周期表フォールバック
    class SimplePeriodicTable:
        def GetAtomicNumber(self, symbol):
            table = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'S': 16, 'P': 15}
            return table.get(symbol, 6)
        def GetAtomicWeight(self, symbol):
            table = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.065, 'P': 30.974}
            return table.get(symbol, 12.0)
    PeriodicTable = SimplePeriodicTable

# Amino acid name mapping / アミノ酸名マッピング
AA_NAME_SYM = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
    'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}

AA_NAME_NUMBER = {k: i for i, (k, _) in enumerate(AA_NAME_SYM.items())}
BACKBONE_NAMES = ["CA", "C", "N", "O"]

ptable = PeriodicTable()


class Residue(dict):
    """
    Amino acid residue / アミノ酸残基
    """
    def __init__(self, name: str = None, atoms: list = None, chain: str = "A",
                 segment: str = "", res_id: int = 1, res_insert_id: str = "",
                 chain_res_id: str = None) -> None:
        self.name = name
        self.atoms = atoms or []
        self.chain = chain
        self.segment = segment
        self.res_id = res_id
        self.res_insert_id = res_insert_id
        self.chain_res_id = chain_res_id


def _enumerate_pdb_lines(lines: List[str]) -> Iterator[Dict[str, Any]]:
    """
    Parse PDB file lines.
    PDBファイル行を解析。
    """
    for line in lines:
        if len(line) < 6:
            continue
        record_type = line[0:6].strip()
        
        if record_type == 'ATOM':
            element_symb = line[76:78].strip().capitalize() if len(line) >= 78 else ""
            if not element_symb:
                element_symb = line[13:14] if len(line) > 13 else "C"
            if element_symb == 'D':
                element_symb = 'H'
            
            yield {
                'type': 'ATOM',
                'atom_id': int(line[6:11]) if len(line) >= 11 else 0,
                'atom_name': line[12:16].strip() if len(line) >= 16 else "",
                'res_name': line[17:20].strip() if len(line) >= 20 else "",
                'chain': line[21:22].strip() if len(line) >= 22 else "A",
                'res_id': int(line[22:26]) if len(line) >= 26 else 0,
                'res_insert_id': line[26:27].strip() if len(line) >= 27 else "",
                'x': float(line[30:38]) if len(line) >= 38 else 0.0,
                'y': float(line[38:46]) if len(line) >= 46 else 0.0,
                'z': float(line[46:54]) if len(line) >= 54 else 0.0,
                'occupancy': float(line[54:60]) if len(line) >= 60 else 1.0,
                'segment': line[72:76].strip() if len(line) >= 76 else "",
                'element_symb': element_symb,
                'charge': line[78:80].strip() if len(line) >= 80 else "",
            }
        elif record_type == 'HEADER':
            yield {'type': 'HEADER', 'value': line[10:].strip() if len(line) > 10 else ""}
        elif record_type == 'ENDMDL':
            break


class Protein:
    """
    Protein data structure / タンパク質データ構造
    
    Supports FASTA/PDB import/export and sequence alignment.
    FASTA/PDBのインポート/エクスポートとシーケンスアラインメントをサポート。
    """
    
    def __init__(self) -> None:
        self.name = None
        self.sequence = None
        self.residues = []
        self.all_atom = []
        self.chi_angles = None
        self.description = None
        self.kg_accession = None

    @classmethod
    def from_fasta(cls, fasta: str, name: str = None) -> Self:
        """
        Create protein from FASTA sequence.
        FASTAシーケンスからタンパク質を作成。
        """
        protein = cls()
        protein.sequence = fasta.upper()
        if name:
            protein.name = name
        else:
            protein._add_name()
        return protein

    @classmethod
    def from_pdb(cls, pdb_lines: List[str], removeHs: bool = True) -> Self:
        """
        Create protein from PDB file content.
        PDBファイル内容からタンパク質を作成。
        """
        protein = cls()
        protein._add_name()
        protein.residues = []
        protein.all_atom = []
        residues_tmp = {}
        
        for data in _enumerate_pdb_lines(pdb_lines):
            if data['type'] == 'HEADER':
                protein.description = data['value'].lower()
                continue
            if removeHs and data['element_symb'] == 'H':
                continue
            if data['res_name'] not in AA_NAME_NUMBER:
                continue
            
            atom_data = {
                "pos": np.array([data['x'], data['y'], data['z']]),
                "atom_name": data['atom_name'],
                "occupancy": data['occupancy'],
                "atomic_number": ptable.GetAtomicNumber(data['element_symb']),
                "weight": ptable.GetAtomicWeight(data['element_symb']),
                "aa_type": AA_NAME_NUMBER[data['res_name']],
                "is_backbone": data['atom_name'] in BACKBONE_NAMES,
            }
            protein.all_atom.append(atom_data)

            chain_res_id = f"{data['chain']}_{data['segment']}_{data['res_id']}_{data['res_insert_id']}"
            if chain_res_id not in residues_tmp:
                residues_tmp[chain_res_id] = Residue(
                    name=data['res_name'],
                    atoms=[len(protein.all_atom) - 1],
                    chain=data['chain'],
                    segment=data['segment'],
                    res_id=data['res_id'],
                    res_insert_id=data['res_insert_id'],
                    chain_res_id=chain_res_id
                )
            else:
                residues_tmp[chain_res_id].atoms.append(len(protein.all_atom) - 1)
        
        # Build sequence and residue properties / シーケンスと残基プロパティを構築
        protein.sequence = ""
        for residue in residues_tmp.values():
            protein.residues.append(residue)
            protein.sequence += AA_NAME_SYM.get(residue.name, 'X')
            
            # Calculate center of mass and backbone atom positions
            # 質量中心とバックボーン原子位置を計算
            sum_pos, sum_mass = np.zeros([3], dtype=np.float32), 0
            for atom_idx in residue.atoms:
                atom = protein.all_atom[atom_idx]
                sum_pos += atom["pos"] * atom["weight"]
                sum_mass += atom["weight"]
                if atom["atom_name"] in BACKBONE_NAMES:
                    setattr(protein.residues[-1], f'pos_{atom["atom_name"]}', atom["pos"])
            protein.residues[-1].center_of_mass = sum_pos / sum_mass
        
        return protein

    @classmethod
    def from_pdb_file(cls, pdb_file: str, removeHs: bool = True) -> Self:
        """
        Create protein from PDB file.
        PDBファイルからタンパク質を作成。
        """
        with open(pdb_file, "r") as f:
            protein = Protein.from_pdb(f.readlines(), removeHs=removeHs)
        protein.name = os.path.basename(pdb_file).replace(".pdb", "")
        return protein

    @classmethod
    def from_binary_file(cls, file: str) -> Self:
        """
        Load protein from binary file.
        バイナリファイルからタンパク質を読み込む。
        """
        return pickle.load(open(file, "rb"))

    def _add_name(self) -> None:
        """
        Auto-generate protein name.
        タンパク質名を自動生成。
        """
        if self.name is None:
            self.name = "protein_" + re.sub(r"[-:.]", "_", 
                                            datetime.now().isoformat(sep="_", timespec="milliseconds"))

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

    def save_pdb(self, file: Optional[str] = None, overwrite: bool = False,
                 residue_indices: Optional[List[int]] = None) -> str:
        """
        Save as PDB file.
        PDBファイルとして保存。
        """
        self._add_name()
        if file is None:
            file = f"./{self.name}.pdb"

        if not os.path.exists(file) or overwrite:
            atom_cnt = 0
            with open(file, "w") as f:
                for i, residue in enumerate(self.residues):
                    if residue_indices is not None and i not in residue_indices:
                        continue
                    for atom_id in residue.atoms:
                        atom = self.all_atom[atom_id]
                        atom_cnt += 1
                        atom_name = f"  {atom['atom_name']:<3}" if len(atom['atom_name']) <= 3 else f" {atom['atom_name']:<4}"
                        elem = atom['atom_name'][0] if atom['atom_name'] else 'C'
                        f.write(f"ATOM  {atom_cnt:5}{atom_name} {residue.name:3} {residue.chain:1}{residue.res_id:4}{residue.res_insert_id:1}   "
                               f"{atom['pos'][0]:8.3f}{atom['pos'][1]:8.3f}{atom['pos'][2]:8.3f}"
                               f"{atom.get('occupancy', 1.00):6.2f}{atom.get('temp_factor', 0.00):6.2f}           {elem}\n")
        return file

    def save_fasta(self, file: Optional[str] = None, overwrite: bool = False) -> str:
        """
        Save as FASTA file.
        FASTAファイルとして保存。
        """
        self._add_name()
        if file is None:
            file = f"./{self.name}.fasta"
        if not os.path.exists(file) or overwrite:
            with open(file, "w") as f:
                f.write(f">{self.name}\n")
                # 60 characters per line / 1行60文字
                for i in range(0, len(self.sequence), 60):
                    f.write(self.sequence[i:i+60] + "\n")
        return file

    def __len__(self) -> int:
        if self.sequence is not None:
            return len(self.sequence)
        if self.residues:
            return len(self.residues)
        return 0

    def __str__(self) -> str:
        return self.sequence if self.sequence else "Protein()"

    def __repr__(self) -> str:
        seq = self.sequence[:20] + "..." if self.sequence and len(self.sequence) > 20 else self.sequence
        return f"Protein(sequence='{seq}', length={len(self)})"


def protein_sequence_similarity(protein1: Protein, protein2: Protein) -> Tuple[float, str, str]:
    """
    Calculate sequence similarity between two proteins.
    2つのタンパク質間のシーケンス類似度を計算。
    
    Returns:
        (similarity score, aligned sequence 1, aligned sequence 2)
        (類似度スコア, アライメント済みシーケンス1, アライメント済みシーケンス2)
    """
    try:
        from Bio.Align import PairwiseAligner
    except ImportError:
        raise ImportError("Biopython required for sequence alignment. "
                         "Install with: pip install biopython")
    
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    
    alignment = aligner.align(protein1.sequence, protein2.sequence)[0]
    
    # Calculate similarity / 類似度を計算
    seq_len = max(len(protein1), len(protein2))
    similarity = alignment.score / seq_len if seq_len > 0 else 0
    
    # Format output / 出力をフォーマット
    aligned_seq1 = str(alignment[0])
    aligned_seq2 = str(alignment[1])
    
    return similarity, aligned_seq1, aligned_seq2
