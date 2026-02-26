"""
分子処理モジュール / Molecule Processing Module

RDKitベースの分子処理・プロパティ計算ライブラリ
RDKit-based molecule processing and property calculation library

機能 / Features:
- SMILES/SELFIES/SDF形式の読み書き / Read/write SMILES/SELFIES/SDF formats
- 分子プロパティ計算 (QED, SA, LogP, Lipinski) / Molecular property calculation
- 分子フィンガープリントと類似度 / Molecular fingerprint and similarity
"""

from typing import List, Optional, Tuple
from typing_extensions import Self

import copy
from datetime import datetime
import gzip
import math
import numpy as np
import os
import pickle
import re

# RDKitインポート / RDKit import
try:
    from rdkit import Chem, DataStructs, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors, Descriptors, Lipinski
    from rdkit.Chem.AllChem import RWMol
    from rdkit.six import iteritems
    from rdkit.six.moves import cPickle
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    raise ImportError("RDKit is required. Install with: conda install -c conda-forge rdkit")

# グローバル変数：合成可及性スコアデータ / Global variable: SA score data
_fscores = None


def _load_fscores():
    """
    Load fragment score data required for SA score calculation.
    SAスコア計算に必要なフラグメントスコアデータを読み込む。
    """
    global _fscores
    if _fscores is not None:
        return _fscores
    
    # 複数の場所から読み込みを試行 / Try loading from multiple locations
    fpscores_paths = [
        os.path.join(os.path.dirname(__file__), "..", "data", "fpscores.pkl.gz"),
        "./data/fpscores.pkl.gz",
    ]
    
    data = []
    for path in fpscores_paths:
        if os.path.exists(path):
            data = cPickle.load(gzip.open(path, "rb"))
            break
    
    _fscores = {}
    for i in data:
        for j in range(1, len(i)):
            _fscores[i[j]] = float(i[0])
    
    return _fscores


class Molecule:
    """
    Molecule data structure / 分子データ構造
    
    Supports SMILES/SELFIES/RDKit Mol interconversion,
    2D/3D conformer generation, and property calculation (QED, SA, LogP, Lipinski).
    
    SMILES/SELFIES/RDKit Molの相互変換、2D/3Dコンフォーマー生成、
    プロパティ計算（QED, SA, LogP, Lipinski）をサポート。
    """
    
    def __init__(self) -> None:
        self.name = None
        self.smiles = None
        self.selfies = None
        self.rdmol = None
        self.graph = None
        self.conformer = None
        self.img = None
        self.description = None
        self.kg_accession = None

    @classmethod
    def from_smiles(cls, smiles: str) -> Self:
        """
        Create molecule from SMILES string.
        SMILES文字列から分子を作成。
        """
        molecule = cls()
        molecule.smiles = smiles
        return molecule

    @classmethod
    def from_selfies(cls, selfies: str) -> Self:
        """
        Create molecule from SELFIES string.
        SELFIES文字列から分子を作成。
        """
        try:
            import selfies as sf
        except ImportError:
            raise ImportError("selfies package required. Install with: pip install selfies")
        
        molecule = cls()
        molecule.selfies = selfies
        molecule.smiles = sf.decoder(selfies)
        return molecule

    @classmethod
    def from_rdmol(cls, rdmol: RWMol) -> Self:
        """
        Create molecule from RDKit mol object.
        RDKit molオブジェクトから分子を作成。
        """
        molecule = cls()
        molecule.rdmol = rdmol
        molecule.smiles = Chem.MolToSmiles(rdmol)
        conformer = rdmol.GetConformer()
        if conformer is not None:
            molecule.conformer = np.array(conformer.GetPositions())
        return molecule

    @classmethod
    def from_sdf_file(cls, sdf_file: str) -> Self:
        """
        Create molecule from SDF file.
        SDFファイルから分子を作成。
        """
        loader = Chem.SDMolSupplier(sdf_file)
        molecule = None
        for mol in loader:
            if mol is not None:
                molecule = Molecule.from_rdmol(mol)
                conformer = mol.GetConformer()
                molecule.conformer = np.array(conformer.GetPositions())
        if molecule is None:
            raise ValueError(f"No valid molecule found in {sdf_file}")
        molecule.name = os.path.basename(sdf_file).replace(".sdf", "")
        return molecule

    @classmethod
    def from_binary_file(cls, file: str) -> Self:
        """
        Load molecule from binary file.
        バイナリファイルから分子を読み込む。
        """
        return pickle.load(open(file, "rb"))

    def _add_name(self) -> None:
        """
        Auto-generate molecule name.
        分子名を自動生成。
        """
        if self.name is None:
            self.name = "mol_" + re.sub(r"[-:.]", "_", datetime.now().isoformat(sep="_", timespec="milliseconds"))

    def _add_rdmol(self) -> None:
        """
        Generate RDKit mol object from SMILES.
        SMILESからRDKit molオブジェクトを生成。
        """
        if self.rdmol is not None:
            return
        if self.smiles:
            self.rdmol = Chem.MolFromSmiles(self.smiles)
            if self.rdmol is None:
                raise ValueError(f"Invalid SMILES: {self.smiles}")
        if self.conformer is not None and self.rdmol is not None:
            conf = _mol_array_to_conformer(self.conformer)
            self.rdmol.AddConformer(conf)

    def _add_conformer(self, mode: str = '2D') -> None:
        """
        Generate molecular conformer.
        分子コンフォーマーを生成。
        
        mode: '2D' or '3D'
        """
        if self.conformer is None:
            self._add_rdmol()
            if mode == '2D':
                AllChem.Compute2DCoords(self.rdmol)
            elif mode == '3D':
                self.rdmol = Chem.AddHs(self.rdmol)
                AllChem.EmbedMolecule(self.rdmol)
                AllChem.MMFFOptimizeMolecule(self.rdmol)
            conformer = self.rdmol.GetConformer()
            self.conformer = np.array(conformer.GetPositions())

    # ========== Property Calculation / プロパティ計算 ==========
    
    def calc_qed(self) -> float:
        """
        Calculate QED (Quantitative Estimate of Drug-likeness).
        QED（定量的薬物相似性推定）を計算。
        
        Returns: Score between 0-1, higher is more drug-like.
                 0-1のスコア、高いほど薬物らしい。
        """
        try:
            from rdkit.Chem.QED import qed
            self._add_rdmol()
            return qed(self.rdmol)
        except Exception:
            return 0.0

    def calc_sa(self) -> float:
        """
        Calculate SA (Synthetic Accessibility) score.
        SA（合成可及性）スコアを計算。
        
        Returns: Score between 0-10, higher is easier to synthesize.
                 0-10のスコア、高いほど合成しやすい。
        """
        self._add_rdmol()
        sa = _calc_sa_score(self.rdmol)
        sa_norm = round((10 - sa) / 9, 2)
        return sa_norm

    def calc_logp(self) -> float:
        """
        Calculate LogP (partition coefficient).
        LogP（分配係数）を計算。
        
        Measures lipophilicity / 親脂性を測定。
        """
        from rdkit.Chem.Crippen import MolLogP
        self._add_rdmol()
        return MolLogP(self.rdmol)

    def calc_lipinski(self) -> int:
        """
        Calculate Lipinski's Rule of Five satisfaction count.
        Lipinskiのルールオブファイブの満足数を計算。
        
        Returns: Number of rules satisfied (0-5), typically >=4 indicates good oral bioavailability.
                 満足したルール数（0-5）、通常>=4は良好な経口生体利用率を示す。
        """
        try:
            self._add_rdmol()
            mol = copy.deepcopy(self.rdmol)
            Chem.SanitizeMol(mol)
            rule_1 = Descriptors.ExactMolWt(mol) < 500
            rule_2 = Lipinski.NumHDonors(mol) <= 5
            rule_3 = Lipinski.NumHAcceptors(mol) <= 10
            logp = self.calc_logp()
            rule_4 = (logp >= -2) & (logp <= 5)
            rule_5 = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10
            return int(sum([rule_1, rule_2, rule_3, rule_4, rule_5]))
        except Exception:
            return 0

    def calc_distance(self) -> np.ndarray:
        """
        Calculate inter-atomic distance matrix.
        原子間距離行列を計算。
        """
        self._add_conformer()
        pdist = self.conformer[None, :] - self.conformer[:, None]
        return np.sqrt(np.sum(pdist ** 2, axis=-1))

    def get_num_atoms(self) -> int:
        """
        Get number of atoms.
        原子数を取得。
        """
        self._add_rdmol()
        return self.rdmol.GetNumAtoms()

    # ========== Save Methods / 保存メソッド ==========

    def save_sdf(self, file: Optional[str] = None, overwrite: bool = False) -> str:
        """
        Save as SDF file.
        SDFファイルとして保存。
        """
        self._add_name()
        if file is None:
            file = f"./{self.name}.sdf"

        if not os.path.exists(file) or overwrite:
            writer = Chem.SDWriter(file)
            self._add_rdmol()
            self._add_conformer()
            writer.write(self.rdmol)
        return file

    def save_binary(self, file: Optional[str] = None, overwrite: bool = False) -> str:
        """
        Save as binary file.
        バイナリファイルとして保存。
        """
        self._add_name()
        if file is None:
            file = f"./{self.name}.pkl"

        if not os.path.exists(file) or overwrite:
            pickle.dump(self, open(file, "wb"))
        return file

    def __str__(self) -> str:
        return self.smiles if self.smiles else "Molecule()"

    def __repr__(self) -> str:
        return f"Molecule(smiles='{self.smiles}')"


# ========== Utility Functions / ユーティリティ関数 ==========

def molecule_fingerprint_similarity(mol1: Molecule, mol2: Molecule, 
                                    fingerprint_type: str = "morgan") -> float:
    """
    Calculate fingerprint similarity (Tanimoto coefficient) between two molecules.
    2つの分子間のフィンガープリント類似度（Tanimoto係数）を計算。
    
    Args:
        mol1, mol2: Molecule objects / 分子オブジェクト
        fingerprint_type: Type of fingerprint ("morgan", "rdkit", "maccs")
                          フィンガープリントの種類
    
    Returns:
        Similarity score (0-1) / 類似度スコア（0-1）
    """
    try:
        mol1._add_rdmol()
        mol2._add_rdmol()
        
        if fingerprint_type == "morgan":
            fp1 = AllChem.GetMorganFingerprint(mol1.rdmol, 2)
            fp2 = AllChem.GetMorganFingerprint(mol2.rdmol, 2)
            return DataStructs.TanimotoSimilarity(fp1, fp2)
        elif fingerprint_type == "rdkit":
            fp1 = Chem.RDKFingerprint(mol1.rdmol)
            fp2 = Chem.RDKFingerprint(mol2.rdmol)
        elif fingerprint_type == "maccs":
            fp1 = MACCSkeys.GenMACCSKeys(mol1.rdmol)
            fp2 = MACCSkeys.GenMACCSKeys(mol2.rdmol)
        else:
            raise ValueError(f"Unknown fingerprint type: {fingerprint_type}")
        
        return DataStructs.FingerprintSimilarity(fp1, fp2, metric=DataStructs.TanimotoSimilarity)
    except Exception:
        return 0.0


def check_identical_molecules(mol1: Molecule, mol2: Molecule) -> bool:
    """
    Check if two molecules are identical (based on InChI).
    2つの分子が同一かどうかをチェック（InChIベース）。
    """
    try:
        mol1._add_rdmol()
        mol2._add_rdmol()
        return Chem.MolToInchi(mol1.rdmol) == Chem.MolToInchi(mol2.rdmol)
    except Exception:
        return False


def calc_mol_diversity(mols: List[Molecule]) -> float:
    """
    Calculate diversity of a list of molecules.
    分子リストの多様性を計算。
    """
    dists = []
    for i in range(len(mols)):
        for j in range(i + 1, len(mols)):
            mol1 = mols[i]
            mol2 = mols[j]
            dists.append(1 - molecule_fingerprint_similarity(mol1, mol2, "rdkit"))
    return np.mean(dists) if dists else 0.0


def calc_mol_rmsd(mol1: Molecule, mol2: Molecule) -> float:
    """
    Calculate RMSD between two molecules.
    2つの分子間のRMSDを計算。
    """
    try:
        if mol1.conformer is None or mol2.conformer is None:
            raise ValueError("RMSD calculation requires 3D conformers")
        assert mol1.get_num_atoms() == mol2.get_num_atoms()
        mol1._add_rdmol()
        mol2._add_rdmol()
        return Chem.rdMolAlign.CalcRMS(mol1.rdmol, mol2.rdmol, maxMatches=30000)
    except Exception:
        return 1e4


def _mol_array_to_conformer(conf: np.ndarray) -> Chem.Conformer:
    """
    Convert numpy array to RDKit conformer object.
    numpy配列をRDKitコンフォーマーオブジェクトに変換。
    """
    new_conf = Chem.Conformer(conf.shape[0])
    for i in range(conf.shape[0]):
        new_conf.SetAtomPosition(i, tuple(conf[i]))
    return new_conf


def _calc_sa_score(molecule: Chem.RWMol) -> float:
    """
    Calculate synthetic accessibility score (based on Ertl and Schuffenhauer).
    合成可及性スコアを計算（ErtlとSchuffenhauerの手法ベース）。
    """
    fscores = _load_fscores()
    
    # Fragment score / フラグメントスコア
    fp = rdMolDescriptors.GetMorganFingerprint(molecule, 2)
    fps = fp.GetNonzeroElements()
    score1 = 0.
    nf = 0
    for bitId, v in iteritems(fps):
        nf += v
        score1 += fscores.get(bitId, -4) * v
    score1 /= nf

    # Feature score / 特徴スコア
    nAtoms = molecule.GetNumAtoms()
    nChiralCenters = len(Chem.FindMolChiralCenters(molecule, includeUnassigned=True))
    ri = molecule.GetRingInfo()
    nBridgeheads = rdMolDescriptors.CalcNumBridgeheadAtoms(molecule)
    nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(molecule)
    nMacrocycles = sum(1 for x in ri.AtomRings() if len(x) > 8)

    sizePenalty = nAtoms ** 1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1)
    spiroPenalty = math.log10(nSpiro + 1)
    bridgePenalty = math.log10(nBridgeheads + 1)
    macrocyclePenalty = math.log10(2) if nMacrocycles > 0 else 0.

    score2 = 0. - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

    # Fingerprint density correction / フィンガープリント密度補正
    score3 = 0.
    if nAtoms > len(fps):
        score3 = math.log(float(nAtoms) / len(fps)) * .5

    sascore = score1 + score2 + score3

    # Convert to 1-10 scale / 1-10スケールに変換
    min_val = -4.0
    max_val = 2.5
    sascore = 11. - (sascore - min_val + 1) / (max_val - min_val) * 9.
    if sascore > 8.:
        sascore = 8. + math.log(sascore + 1. - 9.)
    if sascore > 10.:
        sascore = 10.0
    elif sascore < 1.:
        sascore = 1.0

    return sascore
