"""
可視化モジュール / Visualization Module

分子・タンパク質構造の可視化
Visualization of molecular and protein structures
"""

from typing import Optional
import os

from .molecule import Molecule


def _ensure_dir(file_path: str) -> str:
    """
    Ensure directory exists for file path.
    ファイルパスのディレクトリが存在することを確認。
    """
    directory = os.path.dirname(file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)
    return file_path


def draw_molecule_2d(mol: Molecule, output_file: str = "./molecule.png", 
                     size: tuple = (400, 400)) -> str:
    """
    Draw 2D molecular structure.
    2D分子構造を描画。
    
    Args:
        mol: Molecule object / 分子オブジェクト
        output_file: Output file path / 出力ファイルパス
        size: Image size (width, height) / 画像サイズ（幅、高さ）
    
    Returns:
        Absolute path to output file / 出力ファイルの絶対パス
    """
    try:
        from rdkit.Chem import Draw, rdDepictor
    except ImportError:
        raise ImportError("RDKit required for visualization")
    
    mol._add_rdmol()
    rdDepictor.Compute2DCoords(mol.rdmol)
    
    output_file = _ensure_dir(output_file)
    Draw.MolToImageFile(mol.rdmol, output_file, size=size)
    return os.path.abspath(output_file)


def draw_molecules_grid(mols: list, output_file: str = "./molecules.png",
                        mols_per_row: int = 4, size: tuple = (200, 200)) -> str:
    """
    Draw grid of multiple molecules.
    複数分子のグリッドを描画。
    
    Args:
        mols: List of molecule objects / 分子オブジェクトのリスト
        output_file: Output file path / 出力ファイルパス
        mols_per_row: Number of molecules per row / 1行あたりの分子数
        size: Size of each molecule image / 各分子画像のサイズ
    
    Returns:
        Absolute path to output file / 出力ファイルの絶対パス
    """
    try:
        from rdkit.Chem import Draw, rdDepictor
    except ImportError:
        raise ImportError("RDKit required for visualization")
    
    rdmols = []
    for mol in mols:
        mol._add_rdmol()
        rdDepictor.Compute2DCoords(mol.rdmol)
        rdmols.append(mol.rdmol)
    
    output_file = _ensure_dir(output_file)
    img = Draw.MolsToGridImage(rdmols, molsPerRow=mols_per_row, subImgSize=size)
    img.save(output_file)
    return os.path.abspath(output_file)


# 3D visualization requires PyMOL (optional) / 3D可視化にはPyMOLが必要（オプション）
def _check_pymol():
    """Check if PyMOL is available."""
    try:
        import pymol
        return True
    except ImportError:
        return False


def visualize_protein_3d(pdb_file: str, output_file: str = "./protein.png",
                         width: int = 800, height: int = 600) -> str:
    """
    Visualize protein 3D structure using PyMOL.
    PyMOLを使用してタンパク質3D構造を可視化。
    
    Args:
        pdb_file: PDB file path / PDBファイルパス
        output_file: Output image path / 出力画像パス
        width, height: Image dimensions / 画像サイズ
    
    Returns:
        Absolute path to output file / 出力ファイルの絶対パス
    """
    if not _check_pymol():
        raise ImportError("PyMOL required for 3D visualization. "
                         "Install with: conda install -c conda-forge pymol-open-source")
    
    from pymol import cmd
    
    cmd.reinitialize()
    cmd.load(pdb_file, "protein")
    cmd.hide("everything", "protein")
    cmd.show("cartoon", "protein")
    cmd.color("spectrum", "protein")
    cmd.zoom("all")
    
    output_file = _ensure_dir(output_file)
    cmd.png(output_file, width=width, height=height, dpi=300)
    return os.path.abspath(output_file)


def visualize_complex_3d(pdb_file: str, sdf_file: str, 
                         output_file: str = "./complex.png",
                         width: int = 800, height: int = 600) -> str:
    """
    Visualize protein-ligand complex using PyMOL.
    PyMOLを使用してタンパク質-リガンド複合体を可視化。
    
    Args:
        pdb_file: Protein PDB file path / タンパク質PDBファイルパス
        sdf_file: Ligand SDF file path / リガンドSDFファイルパス
        output_file: Output image path / 出力画像パス
        width, height: Image dimensions / 画像サイズ
    
    Returns:
        Absolute path to output file / 出力ファイルの絶対パス
    """
    if not _check_pymol():
        raise ImportError("PyMOL required for 3D visualization")
    
    from pymol import cmd
    
    cmd.reinitialize()
    cmd.load(pdb_file, "protein")
    cmd.load(sdf_file, "ligand")
    
    cmd.hide("everything", "protein")
    cmd.show("cartoon", "protein")
    cmd.color("grey", "protein")
    
    cmd.show("sticks", "ligand")
    cmd.color("cyan", "ligand")
    
    cmd.zoom("ligand", 10)
    
    output_file = _ensure_dir(output_file)
    cmd.png(output_file, width=width, height=height, dpi=300)
    return os.path.abspath(output_file)
