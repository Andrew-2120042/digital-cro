"""
Receptor preparation module for Digital CRO.

This module converts cleaned PDB protein files to docking-ready PDBQT format for AutoDock Vina.

Process:
    PDB → Add polar hydrogens → Assign charges → PDBQT

Key requirements:
    - Add polar hydrogens (N-H, O-H for hydrogen bonding)
    - Gasteiger partial charges
    - Proper PDBQT formatting with AD4 atom types
    - Handle protein-specific residues

Note:
    This is a simplified implementation. For production use, consider:
    - OpenBabel's -p flag for better hydrogen placement
    - reduce tool for optimal H-bond network
    - PDB2PQR for pKa-based protonation states
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
from Bio.PDB import PDBIO, Structure


def get_atom_charge(atom_name: str, residue_name: str) -> float:
    """
    Get simplified Gasteiger charge for protein atoms.

    This is a VERY simplified charge model using typical values for amino acids.
    For production, use proper Gasteiger calculation or force field charges.

    Args:
        atom_name (str): Atom name (e.g., 'CA', 'N', 'O')
        residue_name (str): Residue name (e.g., 'ALA', 'GLY')

    Returns:
        float: Partial charge estimate

    Note:
        AutoDock Vina uses its own force field, so these approximate charges
        are sufficient for the PDBQT format requirement.
    """
    # Backbone atoms - common to all residues
    if atom_name == 'N':
        return -0.47
    elif atom_name == 'CA':
        return 0.07
    elif atom_name == 'C':
        return 0.51
    elif atom_name == 'O':
        return -0.51

    # Charged residues - simplified
    elif residue_name in ['ARG', 'LYS']:  # Positive
        if atom_name in ['NZ', 'NH1', 'NH2']:
            return -0.80
    elif residue_name in ['ASP', 'GLU']:  # Negative
        if atom_name in ['OD1', 'OD2', 'OE1', 'OE2']:
            return -0.80

    # Polar sidechains
    elif atom_name.startswith('O'):
        return -0.40
    elif atom_name.startswith('N'):
        return -0.40
    elif atom_name.startswith('S'):
        return -0.20

    # Default for carbon/other
    return 0.0


def get_ad4_atom_type_protein(atom_name: str, residue_name: str, element: str) -> str:
    """
    Get AutoDock 4 atom type for protein atoms.

    AD4 protein atom types:
    - C: Aliphatic carbon
    - A: Aromatic carbon
    - N: Nitrogen
    - NA: Nitrogen acceptor
    - OA: Oxygen acceptor
    - S: Sulfur
    - SA: Sulfur acceptor
    - HD: Polar hydrogen

    Args:
        atom_name (str): Atom name (e.g., 'CA', 'N', 'O')
        residue_name (str): Residue name (e.g., 'ALA', 'GLY')
        element (str): Element symbol (e.g., 'C', 'N', 'O')

    Returns:
        str: AD4 atom type
    """
    # Hydrogen
    if element == 'H':
        # All hydrogens we add are polar (H-bond donors)
        return 'HD'

    # Backbone nitrogen
    elif atom_name == 'N':
        return 'NA'  # Acceptor

    # Backbone oxygen
    elif atom_name == 'O':
        return 'OA'  # Acceptor

    # Aromatic carbons (PHE, TYR, TRP, HIS)
    elif residue_name in ['PHE', 'TYR', 'TRP', 'HIS']:
        if atom_name in ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ',
                         'CE3', 'CZ2', 'CZ3', 'CH2']:
            return 'A'  # Aromatic

    # All other carbons
    if element == 'C':
        return 'C'

    # Nitrogen
    elif element == 'N':
        return 'NA'

    # Oxygen
    elif element == 'O':
        return 'OA'

    # Sulfur
    elif element == 'S':
        return 'SA'

    # Default
    return 'C'


def should_add_hydrogen(atom_name: str, residue_name: str) -> bool:
    """
    Determine if hydrogen should be added to this atom for docking.

    We only add POLAR hydrogens (H-bond donors):
    - Backbone N-H
    - Sidechain N-H (ARG, LYS, TRP, etc.)
    - Sidechain O-H (SER, THR, TYR)

    Args:
        atom_name (str): Atom name
        residue_name (str): Residue name

    Returns:
        bool: True if hydrogen should be added
    """
    # Backbone nitrogen (all residues except PRO)
    if atom_name == 'N' and residue_name != 'PRO':
        return True

    # Sidechain polar groups
    # Serine/Threonine OH
    if residue_name in ['SER', 'THR'] and atom_name == 'OG':
        return True
    if residue_name == 'THR' and atom_name == 'OG1':
        return True

    # Tyrosine OH
    if residue_name == 'TYR' and atom_name == 'OH':
        return True

    # Lysine NH3+
    if residue_name == 'LYS' and atom_name == 'NZ':
        return True

    # Arginine NH groups
    if residue_name == 'ARG' and atom_name in ['NE', 'NH1', 'NH2']:
        return True

    # Tryptophan NH
    if residue_name == 'TRP' and atom_name == 'NE1':
        return True

    # Histidine NH
    if residue_name == 'HIS' and atom_name in ['ND1', 'NE2']:
        return True

    # Asparagine/Glutamine NH2
    if residue_name == 'ASN' and atom_name == 'ND2':
        return True
    if residue_name == 'GLN' and atom_name == 'NE2':
        return True

    # Cysteine SH
    if residue_name == 'CYS' and atom_name == 'SG':
        return True

    return False


def add_polar_hydrogens_simple(structure: Structure.Structure) -> Structure.Structure:
    """
    Add polar hydrogens to protein structure.

    This is a SIMPLIFIED implementation that adds hydrogens at approximate
    positions. For production use:
    - OpenBabel with -p flag
    - reduce tool
    - Proper geometry optimization

    Args:
        structure (Structure): Biopython Structure object

    Returns:
        Structure: Modified structure with added hydrogens

    Note:
        Hydrogens are placed at ~1.0 Å from heavy atom in direction
        away from bonded neighbors. This is approximate but sufficient
        for rigid docking.
    """
    from Bio.PDB import Atom

    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname()

                for atom in list(residue.get_atoms()):
                    atom_name = atom.get_name()

                    if should_add_hydrogen(atom_name, res_name):
                        # Get atom position
                        atom_pos = atom.get_coord()

                        # Calculate approximate H position
                        # Direction: away from bonded neighbors
                        neighbors_pos = []

                        # Simple heuristic: find nearby atoms
                        for other_atom in residue.get_atoms():
                            if other_atom.get_name() != atom_name:
                                other_pos = other_atom.get_coord()
                                dist = np.linalg.norm(atom_pos - other_pos)
                                if dist < 2.0:  # Bonded distance threshold
                                    neighbors_pos.append(other_pos)

                        # Calculate direction away from neighbors
                        if neighbors_pos:
                            avg_neighbor_pos = np.mean(neighbors_pos, axis=0)
                            direction = atom_pos - avg_neighbor_pos
                            direction = direction / np.linalg.norm(direction)
                        else:
                            # Default direction if no neighbors found
                            direction = np.array([1.0, 0.0, 0.0])

                        # Place hydrogen ~1.0 Å away
                        h_pos = atom_pos + direction * 1.0

                        # Create hydrogen atom
                        h_name = 'H'
                        if atom_name == 'N':
                            h_name = 'H'
                        elif atom_name.startswith('O'):
                            h_name = 'HO'
                        elif atom_name.startswith('N'):
                            h_name = 'HN'
                        elif atom_name.startswith('S'):
                            h_name = 'HS'

                        # Add serial number to make unique
                        serial_number = max([a.get_serial_number() for a in residue.get_atoms()]) + 1

                        h_atom = Atom.Atom(
                            name=h_name,
                            coord=h_pos,
                            bfactor=0.0,
                            occupancy=1.0,
                            altloc=' ',
                            fullname=h_name,
                            serial_number=serial_number,
                            element='H'
                        )

                        # Add to residue
                        try:
                            residue.add(h_atom)
                        except Exception:
                            # Atom might already exist, skip
                            pass

    return structure


def write_receptor_pdbqt(
    structure: Structure.Structure,
    output_path: str,
    add_hydrogens: bool = True
) -> None:
    """
    Write protein receptor to PDBQT format for AutoDock Vina.

    PDBQT format for receptors:
    - Atomic coordinates
    - Partial charges
    - AD4 atom types
    - No ROOT/ENDROOT (those are for ligands only)

    Args:
        structure (Structure): Biopython Structure with cleaned protein
        output_path (str): Output .pdbqt file path
        add_hydrogens (bool): Add polar hydrogens. Default: True

    Raises:
        IOError: If file cannot be written

    Examples:
        >>> from scripts.utils.protein_ops import load_pdb, clean_pdb
        >>> structure = load_pdb('1HSG.pdb')
        >>> cleaned = clean_pdb(structure, remove_water=True, remove_hetero=True)
        >>> write_receptor_pdbqt(cleaned, '1HSG_receptor.pdbqt')
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)

    # Add polar hydrogens if requested
    if add_hydrogens:
        structure = add_polar_hydrogens_simple(structure)

    # Write PDBQT file
    with open(output_path, 'w') as f:
        # Write header
        f.write("REMARK  Receptor prepared for AutoDock Vina\n")

        atom_serial = 1
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_name = residue.get_resname()
                    res_id = residue.get_id()[1]

                    for atom in residue.get_atoms():
                        atom_name = atom.get_name()
                        element = atom.element
                        coord = atom.get_coord()

                        # Get charge
                        charge = get_atom_charge(atom_name, res_name)

                        # Get AD4 atom type
                        ad4_type = get_ad4_atom_type_protein(atom_name, res_name, element)

                        # PDBQT ATOM line format (standard PDB columns):
                        # ATOM    serial  name  resName chainID resSeq    X       Y       Z     occupancy tempFactor charge type
                        # Columns: 1-6    7-11  13-16  18-20   22   23-26  31-38   39-46   47-54
                        line = f"ATOM  {atom_serial:5d} {atom_name:^4s} {res_name:>3s} {chain.id}{res_id:4d}    "
                        line += f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                        line += f"  1.00  0.00    {charge:6.3f} {ad4_type}\n"

                        f.write(line)
                        atom_serial += 1


def prepare_receptor_for_docking(
    pdb_path: str,
    output_path: str,
    remove_water: bool = True,
    remove_hetero: bool = True,
    keep_chains: Optional[List[str]] = None,
    add_hydrogens: bool = True
) -> Dict:
    """
    Complete receptor preparation pipeline: PDB → PDBQT.

    Args:
        pdb_path (str): Input PDB file path
        output_path (str): Output PDBQT file path
        remove_water (bool): Remove water molecules. Default: True
        remove_hetero (bool): Remove heteroatoms (ligands, ions). Default: True
        keep_chains (List[str], optional): Specific chains to keep. None = all chains
        add_hydrogens (bool): Add polar hydrogens. Default: True

    Returns:
        dict: Result dictionary containing:
            - success (bool): Whether preparation succeeded
            - pdbqt_path (str): Path to output PDBQT file
            - pdb_path (str): Path to input PDB file
            - num_atoms (int): Number of atoms in receptor
            - num_residues (int): Number of residues
            - chains (List[str]): Chain IDs included
            - error (str): Error message if failed

    Examples:
        >>> result = prepare_receptor_for_docking('1HSG.pdb', '1HSG_receptor.pdbqt')
        >>> if result['success']:
        ...     print(f"Prepared receptor: {result['num_atoms']} atoms")
    """
    from scripts.utils.protein_ops import load_pdb, clean_pdb, get_protein_info

    result = {
        'success': False,
        'pdbqt_path': output_path,
        'pdb_path': pdb_path,
        'num_atoms': 0,
        'num_residues': 0,
        'chains': [],
        'error': None
    }

    try:
        # Step 1: Load PDB
        structure = load_pdb(pdb_path)

        if structure is None:
            result['error'] = 'Failed to load PDB file'
            return result

        # Step 2: Clean PDB
        cleaned_structure = clean_pdb(
            structure,
            remove_water=remove_water,
            remove_hetero=remove_hetero,
            keep_chains=keep_chains
        )

        # Step 3: Get protein info
        info = get_protein_info(cleaned_structure)
        result['num_residues'] = info['num_residues']
        result['chains'] = info['chains']

        # Step 4: Write PDBQT
        write_receptor_pdbqt(cleaned_structure, output_path, add_hydrogens=add_hydrogens)

        # Count final atoms (after H addition)
        atom_count = 0
        for model in cleaned_structure:
            for chain in model:
                for residue in chain:
                    atom_count += len(list(residue.get_atoms()))

        result['num_atoms'] = atom_count
        result['success'] = True

    except Exception as e:
        result['error'] = str(e)

    return result
