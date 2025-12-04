"""
Protein structure operations module for Digital CRO.

This module provides functions for downloading, loading, cleaning, and analyzing
protein structures from PDB files.
"""

import os
import requests
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import warnings

from Bio import BiopythonWarning
from Bio.PDB import PDBParser, PDBIO, Structure, Select
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import numpy as np


# Suppress Biopython warnings for cleaner output
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', PDBConstructionWarning)


def download_pdb(pdb_id: str, output_dir: str = 'data/proteins/') -> str:
    """
    Download PDB file from RCSB database.

    Args:
        pdb_id (str): 4-letter PDB code (e.g., '1HSG', '6LU7')
        output_dir (str): Directory to save the file

    Returns:
        str: Path to downloaded PDB file

    Raises:
        ValueError: If PDB ID is invalid format
        requests.HTTPError: If download fails (PDB not found)
        IOError: If file cannot be written

    Examples:
        >>> path = download_pdb('1HSG')
        >>> print(path)
        'data/proteins/1HSG.pdb'
    """
    # Validate PDB ID format
    pdb_id = pdb_id.strip().upper()
    if len(pdb_id) != 4:
        raise ValueError(f"Invalid PDB ID '{pdb_id}'. Must be 4 characters.")

    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Construct URL and output path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = output_dir / f"{pdb_id}.pdb"

    try:
        # Download file
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        # Check if response is valid PDB file
        if 'html' in response.headers.get('Content-Type', '').lower():
            raise requests.HTTPError(f"PDB ID '{pdb_id}' not found in RCSB database")

        # Save to file
        with open(output_path, 'w') as f:
            f.write(response.text)

        return str(output_path)

    except requests.HTTPError as e:
        raise requests.HTTPError(f"Failed to download PDB {pdb_id}: {str(e)}")
    except Exception as e:
        raise IOError(f"Error saving PDB file: {str(e)}")


def load_pdb(filepath: str) -> Structure.Structure:
    """
    Load PDB file using Biopython.

    Args:
        filepath (str): Path to PDB file

    Returns:
        Bio.PDB.Structure: Biopython Structure object

    Raises:
        FileNotFoundError: If PDB file doesn't exist
        ValueError: If PDB file is malformed

    Notes:
        - If multiple models exist, all are loaded but you can access model 0
        - Warnings about discontinuous chains are suppressed

    Examples:
        >>> structure = load_pdb('data/proteins/1HSG.pdb')
        >>> print(structure.id)
        '1HSG'
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"PDB file not found: {filepath}")

    try:
        parser = PDBParser(QUIET=True)
        structure_id = Path(filepath).stem
        structure = parser.get_structure(structure_id, filepath)

        # Validate structure has at least one model
        if len(list(structure.get_models())) == 0:
            raise ValueError("PDB file contains no models")

        return structure

    except Exception as e:
        raise ValueError(f"Error parsing PDB file: {str(e)}")


class CleanStructureSelect(Select):
    """
    Custom selector for cleaning PDB structures.
    """

    def __init__(self, remove_water=True, remove_hetero=True, keep_chains=None):
        self.remove_water = remove_water
        self.remove_hetero = remove_hetero
        self.keep_chains = set(keep_chains) if keep_chains else None

    def accept_chain(self, chain):
        """Accept chain if in keep_chains list (or if None)."""
        if self.keep_chains is None:
            return True
        return chain.id in self.keep_chains

    def accept_residue(self, residue):
        """Accept residue based on filtering criteria."""
        hetero_flag, resseq, icode = residue.id

        # Check chain
        if self.keep_chains is not None:
            if residue.parent.id not in self.keep_chains:
                return False

        # Remove water molecules
        if self.remove_water:
            if residue.resname in ['HOH', 'WAT', 'H2O']:
                return False

        # Remove heteroatoms (ligands, ions, etc.)
        if self.remove_hetero:
            if hetero_flag != ' ':  # ' ' indicates standard amino acid
                return False

        return True


def clean_pdb(structure: Structure.Structure,
              remove_water: bool = True,
              remove_hetero: bool = True,
              keep_chains: Optional[List[str]] = None) -> Structure.Structure:
    """
    Clean PDB structure for docking preparation.

    Args:
        structure (Structure): Biopython Structure object
        remove_water (bool): Remove water molecules (HOH, WAT). Default: True
        remove_hetero (bool): Remove heteroatoms (ligands, ions). Default: True
        keep_chains (List[str], optional): Chain IDs to keep. None keeps all chains.

    Returns:
        Structure: Cleaned Structure object (new copy)

    Notes:
        Water residues removed: HOH, WAT, H2O
        Heteroatoms: Any residue with hetero flag != ' '

    Examples:
        >>> structure = load_pdb('1HSG.pdb')
        >>> cleaned = clean_pdb(structure, remove_water=True, keep_chains=['A'])
        >>> # Keep only chain A, remove water and ligands
    """
    # Create a deep copy to avoid modifying original
    import copy
    cleaned_structure = copy.deepcopy(structure)

    # Use custom selector to filter
    selector = CleanStructureSelect(
        remove_water=remove_water,
        remove_hetero=remove_hetero,
        keep_chains=keep_chains
    )

    # Remove unwanted residues
    for model in cleaned_structure:
        chains_to_remove = []
        for chain in model:
            residues_to_remove = []
            for residue in chain:
                if not selector.accept_residue(residue):
                    residues_to_remove.append(residue.id)

            # Remove residues
            for res_id in residues_to_remove:
                chain.detach_child(res_id)

            # Mark empty chains for removal
            if len(chain) == 0:
                chains_to_remove.append(chain.id)

        # Remove empty chains
        for chain_id in chains_to_remove:
            model.detach_child(chain_id)

    return cleaned_structure


def save_pdb(structure: Structure.Structure, output_path: str) -> None:
    """
    Save structure to PDB file.

    Args:
        structure (Structure): Biopython Structure object
        output_path (str): Path to save PDB file

    Raises:
        IOError: If file cannot be written

    Examples:
        >>> save_pdb(cleaned_structure, 'data/proteins/1HSG_clean.pdb')
    """
    try:
        # Create output directory if needed
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        io = PDBIO()
        io.set_structure(structure)
        io.save(output_path)

    except Exception as e:
        raise IOError(f"Error saving PDB file: {str(e)}")


def get_protein_info(structure: Structure.Structure) -> Dict[str, any]:
    """
    Extract metadata from protein structure.

    Args:
        structure (Structure): Biopython Structure object

    Returns:
        dict: Dictionary containing:
            - pdb_id (str): PDB identifier
            - num_chains (int): Number of chains
            - num_residues (int): Total residue count
            - chains (List[str]): List of chain IDs
            - num_atoms (int): Total atom count
            - models (int): Number of models

    Examples:
        >>> info = get_protein_info(structure)
        >>> print(f"Chains: {info['chains']}")
        Chains: ['A', 'B']
    """
    info = {
        'pdb_id': structure.id,
        'num_chains': 0,
        'num_residues': 0,
        'chains': [],
        'num_atoms': 0,
        'models': len(list(structure.get_models()))
    }

    # Use first model
    model = structure[0]

    chains = list(model.get_chains())
    info['num_chains'] = len(chains)
    info['chains'] = [chain.id for chain in chains]

    # Count residues and atoms
    residues = list(model.get_residues())
    info['num_residues'] = len(residues)

    atoms = list(model.get_atoms())
    info['num_atoms'] = len(atoms)

    return info


def calculate_center_of_mass(structure: Structure.Structure) -> Tuple[float, float, float]:
    """
    Calculate center of mass of protein structure.

    Args:
        structure (Structure): Biopython Structure object

    Returns:
        Tuple[float, float, float]: (x, y, z) coordinates in Angstroms

    Notes:
        - Uses all atoms in first model
        - Assumes uniform mass (doesn't use actual atomic masses)
        - Useful for initial binding site guess

    Examples:
        >>> com = calculate_center_of_mass(structure)
        >>> print(f"Center: {com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f}")
        Center: 2.30, 15.80, 24.10
    """
    # Get all atoms from first model
    model = structure[0]
    atoms = list(model.get_atoms())

    if len(atoms) == 0:
        raise ValueError("Structure contains no atoms")

    # Extract coordinates
    coords = np.array([atom.coord for atom in atoms])

    # Calculate center of mass (assuming uniform mass)
    center = coords.mean(axis=0)

    return tuple(center)


def get_chain_sequence(structure: Structure.Structure, chain_id: str) -> str:
    """
    Get amino acid sequence for a specific chain.

    Args:
        structure (Structure): Biopython Structure object
        chain_id (str): Chain identifier (e.g., 'A')

    Returns:
        str: Single-letter amino acid sequence

    Raises:
        ValueError: If chain not found

    Examples:
        >>> seq = get_chain_sequence(structure, 'A')
        >>> print(seq[:20])
        PQITLWQRPLVTIKIGGQLK
    """
    from Bio.PDB.Polypeptide import PPBuilder

    model = structure[0]

    if chain_id not in [c.id for c in model.get_chains()]:
        raise ValueError(f"Chain '{chain_id}' not found in structure")

    chain = model[chain_id]
    ppb = PPBuilder()
    peptides = ppb.build_peptides(chain)

    # Concatenate all peptides in chain
    sequence = ''.join([str(pp.get_sequence()) for pp in peptides])

    return sequence
