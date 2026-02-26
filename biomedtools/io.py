"""
IOツールモジュール / IO Tools Module

フォーマット変換とデータ入出力
Format conversion and data import/export
"""

from typing import List, Union, Optional
import os

from .molecule import Molecule
from .protein import Protein
from .pocket import Pocket


def load_smiles_file(file_path: str) -> List[Molecule]:
    """
    Load molecules from SMILES file.
    SMILESファイルから分子を読み込む。
    
    File format: one SMILES per line / ファイル形式：1行に1つのSMILES
    """
    molecules = []
    with open(file_path, 'r') as f:
        for line in f:
            smiles = line.strip()
            if smiles and not smiles.startswith('#'):
                try:
                    mol = Molecule.from_smiles(smiles)
                    molecules.append(mol)
                except Exception:
                    continue
    return molecules


def load_sdf_file(file_path: str) -> List[Molecule]:
    """
    Load molecules from SDF file (supports multiple molecules).
    SDFファイルから分子を読み込む（複数分子対応）。
    """
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError("RDKit required for SDF loading")
    
    molecules = []
    supplier = Chem.SDMolSupplier(file_path)
    for rdmol in supplier:
        if rdmol is not None:
            mol = Molecule.from_rdmol(rdmol)
            molecules.append(mol)
    return molecules


def save_smiles_file(molecules: List[Molecule], file_path: str) -> str:
    """
    Save molecule list as SMILES file.
    分子リストをSMILESファイルとして保存。
    """
    with open(file_path, 'w') as f:
        for mol in molecules:
            if mol.smiles:
                f.write(mol.smiles + '\n')
    return os.path.abspath(file_path)


def load_fasta_file(file_path: str) -> List[Protein]:
    """
    Load protein sequences from FASTA file.
    FASTAファイルからタンパク質シーケンスを読み込む。
    
    Supports multi-sequence FASTA format.
    マルチシーケンスFASTA形式をサポート。
    """
    proteins = []
    current_name = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence / 前のシーケンスを保存
                if current_name and current_seq:
                    seq = ''.join(current_seq)
                    protein = Protein.from_fasta(seq, name=current_name)
                    proteins.append(protein)
                # Start new sequence / 新しいシーケンスを開始
                current_name = line[1:].split()[0]
                current_seq = []
            elif line:
                current_seq.append(line)
    
    # Save last sequence / 最後のシーケンスを保存
    if current_name and current_seq:
        seq = ''.join(current_seq)
        protein = Protein.from_fasta(seq, name=current_name)
        proteins.append(protein)
    
    return proteins


def apply_mutation(protein: Protein, mutation: str) -> Protein:
    """
    Apply point mutation to protein sequence.
    タンパク質シーケンスに点突然変異を適用。
    
    Args:
        protein: Original protein / 元のタンパク質
        mutation: Mutation string (e.g., "A123V" - position 123 changes from A to V)
                  変異文字列（例："A123V" - 123番目がAからVに変化）
    
    Returns:
        New mutated protein / 変異後の新しいタンパク質
    """
    if len(mutation) < 3:
        raise ValueError("Mutation string must be in format 'A123V'")
    
    # Parse mutation / 変異を解析
    original_aa = mutation[0]
    new_aa = mutation[-1]
    
    # Extract position / 位置を抽出
    pos_str = mutation[1:-1]
    try:
        position = int(pos_str) - 1  # Convert to 0-indexed / 0-indexedに変換
    except ValueError:
        raise ValueError(f"Invalid position in mutation: {mutation}")
    
    if position < 0 or position >= len(protein.sequence):
        raise ValueError(f"Position {position+1} out of range for sequence of length {len(protein)}")
    
    # Verify original amino acid / 元のアミノ酸を検証
    if protein.sequence[position] != original_aa:
        raise ValueError(f"Expected {original_aa} at position {position+1}, found {protein.sequence[position]}")
    
    # Apply mutation / 変異を適用
    new_seq = protein.sequence[:position] + new_aa + protein.sequence[position+1:]
    new_name = f"{protein.name}_{mutation}" if protein.name else f"mutated_{mutation}"
    
    return Protein.from_fasta(new_seq, name=new_name)


class MoleculeExporter:
    """
    Molecule export utility class / 分子エクスポートユーティリティクラス
    """
    
    @staticmethod
    def to_sdf(molecule: Molecule, file_path: Optional[str] = None) -> str:
        """Export to SDF format / SDF形式にエクスポート"""
        return molecule.save_sdf(file_path)
    
    @staticmethod
    def to_pkl(molecule: Molecule, file_path: Optional[str] = None) -> str:
        """Export to binary format / バイナリ形式にエクスポート"""
        return molecule.save_binary(file_path)


class ProteinExporter:
    """
    Protein export utility class / タンパク質エクスポートユーティリティクラス
    """
    
    @staticmethod
    def to_pdb(protein: Protein, file_path: Optional[str] = None) -> str:
        """Export to PDB format / PDB形式にエクスポート"""
        return protein.save_pdb(file_path)
    
    @staticmethod
    def to_fasta(protein: Protein, file_path: Optional[str] = None) -> str:
        """Export to FASTA format / FASTA形式にエクスポート"""
        return protein.save_fasta(file_path)
    
    @staticmethod
    def to_pkl(protein: Protein, file_path: Optional[str] = None) -> str:
        """Export to binary format / バイナリ形式にエクスポート"""
        return protein.save_binary(file_path)
