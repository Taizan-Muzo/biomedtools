import unittest
import sys
import os

# 添加项目到路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biomedtools import Molecule, Protein, Pocket
from biomedtools import molecule_fingerprint_similarity, check_identical_molecules


class TestMolecule(unittest.TestCase):
    """测试分子功能"""
    
    def test_from_smiles(self):
        mol = Molecule.from_smiles("CCO")
        self.assertEqual(mol.smiles, "CCO")
        self.assertEqual(str(mol), "CCO")
    
    def test_num_atoms(self):
        mol = Molecule.from_smiles("CCO")  # 乙醇
        self.assertEqual(mol.get_num_atoms(), 9)  # 2C + 1O + 6H
    
    def test_properties(self):
        mol = Molecule.from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")  # 阿司匹林
        
        # 测试各项性质计算能运行（不验证具体值）
        qed = mol.calc_qed()
        self.assertIsInstance(qed, float)
        self.assertTrue(0 <= qed <= 1)
        
        sa = mol.calc_sa()
        self.assertIsInstance(sa, float)
        
        logp = mol.calc_logp()
        self.assertIsInstance(logp, float)
        
        lipinski = mol.calc_lipinski()
        self.assertIsInstance(lipinski, int)
        self.assertTrue(0 <= lipinski <= 5)
    
    def test_similarity(self):
        mol1 = Molecule.from_smiles("CCO")
        mol2 = Molecule.from_smiles("CCCO")
        
        sim = molecule_fingerprint_similarity(mol1, mol2)
        self.assertIsInstance(sim, float)
        self.assertTrue(0 <= sim <= 1)
    
    def test_identity(self):
        mol1 = Molecule.from_smiles("CCO")
        mol2 = Molecule.from_smiles("CCO")
        mol3 = Molecule.from_smiles("CCCO")
        
        self.assertTrue(check_identical_molecules(mol1, mol2))
        self.assertFalse(check_identical_molecules(mol1, mol3))


class TestProtein(unittest.TestCase):
    """测试蛋白质功能"""
    
    def test_from_fasta(self):
        seq = "MTEYKLVVV"
        protein = Protein.from_fasta(seq)
        self.assertEqual(protein.sequence, seq)
        self.assertEqual(len(protein), 9)
    
    def test_empty_protein(self):
        protein = Protein()
        self.assertEqual(len(protein), 0)


class TestPocket(unittest.TestCase):
    """测试口袋功能"""
    
    def test_basic_pocket(self):
        pocket = Pocket()
        self.assertEqual(pocket.get_num_atoms(), 0)


if __name__ == '__main__':
    unittest.main()
