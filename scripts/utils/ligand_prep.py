"""
Ligand preparation module for Digital CRO.

This module converts SMILES strings to docking-ready PDBQT format for AutoDock Vina.

Process:
    SMILES → 3D RDKit Mol → PDB → PDBQT

Key requirements:
    - 3D coordinates (ETKDG conformer generation)
    - Gasteiger partial charges
    - Rotatable bond identification
    - Proper PDBQT formatting with atom types
"""

import os
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import rdPartialCharges


def smiles_to_3d_mol(smiles: str, optimize: bool = True) -> Optional[Chem.Mol]:
    """
    Convert SMILES to 3D molecule with robust conformer generation.
    Falls back to single conformer if multi-conformer fails.

    Args:
        smiles (str): SMILES string
        optimize (bool): If True, optimize geometry with MMFF force field. Default: True

    Returns:
        Chem.Mol: RDKit molecule with 3D coordinates, or None if failed

    Process:
        1. Parse SMILES
        2. Add hydrogens (essential for 3D structure)
        3. Try multi-conformer generation (10 conformers, select best)
        4. Fallback to single conformer if multi-conformer fails
        5. Optimize geometry with MMFF94 force field

    Examples:
        >>> mol = smiles_to_3d_mol('CCO')  # Ethanol
        >>> if mol:
        ...     print(f"Generated 3D structure with {mol.GetNumAtoms()} atoms")
    """
    import logging
    logger = logging.getLogger(__name__)

    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Failed to parse SMILES: {smiles}")
            return None

        # Add hydrogens
        mol = Chem.AddHs(mol)
        if mol is None:
            logger.error(f"Failed to add hydrogens for: {smiles}")
            return None

        # STRATEGY 1: Try multi-conformer generation (best quality)
        try:
            num_confs = 10
            conf_ids = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=num_confs,
                params=AllChem.ETKDGv3(),
                randomSeed=42,
                numThreads=1,  # Single thread for stability
                pruneRmsThresh=0.5  # Remove similar conformers
            )

            if conf_ids and len(conf_ids) > 0:
                # Success! Optimize each conformer
                for conf_id in conf_ids:
                    try:
                        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, maxIters=500)
                    except:
                        continue

                # Select best (lowest energy) conformer
                energies = []
                for conf_id in conf_ids:
                    try:
                        props = AllChem.MMFFGetMoleculeProperties(mol)
                        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
                        if ff:
                            energy = ff.CalcEnergy()
                            energies.append((conf_id, energy))
                    except:
                        continue

                if energies:
                    # Keep only best conformer
                    best_conf_id = min(energies, key=lambda x: x[1])[0]
                    best_energy = min(energies, key=lambda x: x[1])[1]
                    new_mol = Chem.Mol(mol)
                    new_mol.RemoveAllConformers()
                    new_mol.AddConformer(mol.GetConformer(best_conf_id))
                    logger.debug(f"Generated {len(conf_ids)} conformers, selected best (energy: {best_energy:.2f})")
                    return new_mol

        except Exception as e:
            logger.debug(f"Multi-conformer generation failed: {str(e)}, falling back to single conformer")

        # STRATEGY 2: Fallback to single conformer (original method)
        try:
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))  # Fresh molecule

            # Try ETKDGv3 first (best)
            result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

            if result != 0:  # Failed
                # Try ETKDGv2
                result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())

            if result != 0:  # Still failed
                # Try basic embedding
                result = AllChem.EmbedMolecule(mol)

            if result != 0:  # All methods failed
                logger.error(f"All embedding methods failed for SMILES: {smiles}")
                return None

            # Optimize
            if optimize:
                try:
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
                except:
                    try:
                        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
                    except:
                        pass  # Continue with unoptimized structure

            logger.debug(f"Generated single conformer (fallback method)")
            return mol

        except Exception as e:
            logger.error(f"Single conformer generation also failed: {str(e)}")
            return None

    except Exception as e:
        logger.error(f"Unexpected error in smiles_to_3d_mol: {str(e)}")
        return None


def add_gasteiger_charges(mol: Chem.Mol) -> Chem.Mol:
    """
    Add Gasteiger partial charges to molecule.

    Required for AutoDock Vina scoring function.

    Args:
        mol (Chem.Mol): RDKit molecule

    Returns:
        Chem.Mol: Molecule with '_GasteigerCharge' property on each atom

    Examples:
        >>> mol = Chem.MolFromSmiles('CCO')
        >>> mol = add_gasteiger_charges(mol)
        >>> charge = mol.GetAtomWithIdx(0).GetDoubleProp('_GasteigerCharge')
    """
    # Compute Gasteiger charges
    rdPartialCharges.ComputeGasteigerCharges(mol)
    return mol


def identify_rotatable_bonds(mol: Chem.Mol) -> List[Tuple[int, int]]:
    """
    Identify rotatable bonds for flexibility in docking.

    Rotatable bonds are:
    - Single bonds
    - Not in rings
    - Not terminal (not connected to H or degree-1 atoms)

    Args:
        mol (Chem.Mol): RDKit molecule

    Returns:
        List[Tuple[int, int]]: List of (atom_idx1, atom_idx2) for rotatable bonds

    Examples:
        >>> mol = Chem.MolFromSmiles('CCCC')  # Butane
        >>> bonds = identify_rotatable_bonds(mol)
        >>> print(f"Rotatable bonds: {len(bonds)}")
    """
    from rdkit.Chem import Lipinski

    rot_bonds = []

    # Use SMARTS pattern to find rotatable bonds
    pattern = Lipinski.RotatableBondSmarts
    matches = mol.GetSubstructMatches(pattern)

    for match in matches:
        # Each match is a tuple of atom indices
        # Find the bond between them
        for i in range(len(match) - 1):
            bond = mol.GetBondBetweenAtoms(match[i], match[i+1])
            if bond:
                rot_bonds.append((match[i], match[i+1]))
                break

    return rot_bonds


def get_ad4_atom_type(atom: Chem.Atom) -> str:
    """
    Get AutoDock 4 atom type for an atom.

    AD4 atom types:
    - C: Aliphatic carbon
    - A: Aromatic carbon
    - N: Nitrogen
    - NA: Nitrogen acceptor
    - NS: Nitrogen in S-containing ring
    - OA: Oxygen acceptor
    - OS: Oxygen in S-containing ring
    - S: Sulfur
    - SA: Sulfur acceptor
    - HD: Polar hydrogen (donor)
    - H: Non-polar hydrogen

    Args:
        atom (Chem.Atom): RDKit atom

    Returns:
        str: AD4 atom type

    Note:
        This is a simplified implementation. Full AD4 typing is complex.
    """
    atomic_num = atom.GetAtomicNum()
    is_aromatic = atom.GetIsAromatic()

    # Hydrogen
    if atomic_num == 1:
        # Check if attached to N, O, S (polar hydrogen)
        neighbors = atom.GetNeighbors()
        if neighbors:
            neighbor_num = neighbors[0].GetAtomicNum()
            if neighbor_num in [7, 8, 16]:  # N, O, S
                return 'HD'
        return 'H'

    # Carbon
    elif atomic_num == 6:
        if is_aromatic:
            return 'A'
        return 'C'

    # Nitrogen
    elif atomic_num == 7:
        # Simplified: all nitrogens as NA (acceptor)
        return 'NA'

    # Oxygen
    elif atomic_num == 8:
        # All oxygens as OA (acceptor)
        return 'OA'

    # Sulfur
    elif atomic_num == 16:
        return 'SA'

    # Phosphorus
    elif atomic_num == 15:
        return 'P'

    # Halogen s
    elif atomic_num in [9, 17, 35, 53]:  # F, Cl, Br, I
        if atomic_num == 9:
            return 'F'
        elif atomic_num == 17:
            return 'Cl'
        elif atomic_num == 35:
            return 'Br'
        elif atomic_num == 53:
            return 'I'

    # Default
    return 'C'


def write_pdbqt_file(mol: Chem.Mol, output_path: str, mol_id: str = 'LIG') -> None:
    """
    Write molecule to PDBQT format for AutoDock Vina.

    PDBQT format includes:
    - Atomic coordinates
    - Gasteiger partial charges
    - AD4 atom types
    - ROOT/ENDROOT markers
    - TORSDOF (number of rotatable bonds)

    Args:
        mol (Chem.Mol): RDKit molecule with 3D coordinates and charges
        output_path (str): Output .pdbqt file path
        mol_id (str): Molecule identifier (3 letters). Default: 'LIG'

    Raises:
        ValueError: If molecule has no 3D coordinates
        IOError: If file cannot be written

    Examples:
        >>> mol = smiles_to_3d_mol('CCO')
        >>> mol = add_gasteiger_charges(mol)
        >>> write_pdbqt_file(mol, 'ethanol.pdbqt', 'ETH')
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)

    # Check for 3D coordinates
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no 3D coordinates. Run smiles_to_3d_mol() first.")

    # Get conformer
    conf = mol.GetConformer()

    # Identify rotatable bonds
    rot_bonds = identify_rotatable_bonds(mol)

    # Write PDBQT file
    with open(output_path, 'w') as f:
        # Write header
        f.write(f"REMARK  Name = {mol_id}\n")
        f.write(f"REMARK  {mol.GetNumAtoms()} atoms\n")
        f.write("ROOT\n")

        # Ensure mol_id is exactly 3 characters for PDB format
        res_name = mol_id[:3].upper() if len(mol_id) >= 3 else mol_id.upper().ljust(3)

        # Write atoms
        for i, atom in enumerate(mol.GetAtoms()):
            atom_idx = i + 1
            atom_name = atom.GetSymbol() + str(atom_idx)
            pos = conf.GetAtomPosition(i)

            # Get charge
            try:
                charge = atom.GetDoubleProp('_GasteigerCharge')
            except KeyError:
                charge = 0.0

            # Get AD4 atom type
            ad4_type = get_ad4_atom_type(atom)

            # PDBQT ATOM line format (standard PDB columns):
            # ATOM    serial  name  resName chainID resSeq    X       Y       Z     occupancy tempFactor charge type
            # Columns: 1-6    7-11  13-16  18-20   22   23-26  31-38   39-46   47-54
            line = f"ATOM  {atom_idx:5d} {atom_name:^4s} {res_name:>3s} A   1    "
            line += f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
            line += f"  1.00  0.00    {charge:6.3f} {ad4_type}\n"

            f.write(line)

        f.write("ENDROOT\n")

        # Write torsion degrees of freedom
        f.write(f"TORSDOF {len(rot_bonds)}\n")


def prepare_ligand_for_docking(
    smiles: str,
    output_path: str,
    mol_id: Optional[str] = None
) -> Dict:
    """
    Complete ligand preparation pipeline: SMILES → PDBQT.

    Args:
        smiles (str): Input SMILES string
        output_path (str): Where to save .pdbqt file
        mol_id (str, optional): Molecule identifier. Auto-generated if None.

    Returns:
        dict: Result dictionary containing:
            - success (bool): Whether preparation succeeded
            - pdbqt_path (str): Path to output PDBQT file
            - smiles (str): Input SMILES
            - mol_id (str): Molecule identifier
            - num_atoms (int): Number of atoms
            - num_rotatable_bonds (int): Number of rotatable bonds
            - molecular_weight (float): Molecular weight
            - error (str): Error message if failed

    Examples:
        >>> result = prepare_ligand_for_docking('CCO', 'ethanol.pdbqt')
        >>> if result['success']:
        ...     print(f"Prepared ligand: {result['pdbqt_path']}")
    """
    result = {
        'success': False,
        'pdbqt_path': output_path,
        'smiles': smiles,
        'mol_id': mol_id or 'LIG',
        'num_atoms': 0,
        'num_rotatable_bonds': 0,
        'molecular_weight': 0.0,
        'error': None
    }

    import logging
    logger = logging.getLogger(__name__)

    try:
        # Step 1: Convert SMILES to 3D molecule
        mol = smiles_to_3d_mol(smiles, optimize=True)

        if mol is None:
            error_msg = f'Failed to generate 3D structure from SMILES: {smiles}'
            result['error'] = error_msg
            logger.error(f"❌ [{result['mol_id']}] {error_msg}")
            return result

        # Step 2: Add Gasteiger charges
        try:
            mol = add_gasteiger_charges(mol)
        except Exception as e:
            error_msg = f'Failed to add Gasteiger charges: {str(e)}'
            result['error'] = error_msg
            logger.error(f"❌ [{result['mol_id']}] {error_msg}")
            return result

        # Step 3: Get properties
        result['num_atoms'] = mol.GetNumAtoms()
        result['num_rotatable_bonds'] = len(identify_rotatable_bonds(mol))
        result['molecular_weight'] = Descriptors.MolWt(mol)

        # Step 4: Write PDBQT file
        try:
            write_pdbqt_file(mol, output_path, result['mol_id'])
        except Exception as e:
            error_msg = f'Failed to write PDBQT file: {str(e)}'
            result['error'] = error_msg
            logger.error(f"❌ [{result['mol_id']}] {error_msg}")
            return result

        result['success'] = True
        logger.debug(f"✓ [{result['mol_id']}] Successfully prepared ligand")

    except Exception as e:
        import traceback
        error_msg = f'Unexpected error: {str(e)}'
        result['error'] = error_msg
        logger.error(f"❌ [{result['mol_id']}] {error_msg}")
        logger.error(traceback.format_exc())

    return result


def batch_prepare_ligands(
    smiles_file: str,
    output_dir: str,
    max_workers: int = 4,
    verbose: bool = True
) -> Dict:
    """
    Prepare multiple ligands in parallel.

    Args:
        smiles_file (str): Path to SMILES file (SMILES<tab>ID format)
        output_dir (str): Directory to save .pdbqt files
        max_workers (int): Number of parallel threads. Default: 4
        verbose (bool): Show progress bar. Default: True

    Returns:
        dict: Statistics dictionary containing:
            - total (int): Total molecules processed
            - success (int): Successfully prepared
            - failed (int): Failed preparations
            - output_files (List[str]): Paths to generated PDBQT files
            - failed_smiles (List[str]): SMILES that failed
            - time_elapsed (float): Total time in seconds

    Examples:
        >>> stats = batch_prepare_ligands('library.smi', 'ligands/', max_workers=4)
        >>> print(f"Prepared {stats['success']}/{stats['total']} ligands")
    """
    start_time = time.time()

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Read SMILES file - detect format automatically
    smiles_data = []

    # First, determine if it's CSV or SMI format
    # Skip comment lines and empty lines to find the first real line
    with open(smiles_file, 'r') as f:
        first_line = None
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                first_line = line
                break

        if not first_line:
            raise ValueError("File is empty or contains only comments")

        has_comma = ',' in first_line
        has_tab = '\t' in first_line

        # Check if first line is a header
        is_header = False
        smiles_col_idx = 0
        id_col_idx = 1

        if any(keyword in first_line.lower() for keyword in ['smiles', 'id', 'name', 'mol_id']):
            is_header = True
            # Parse header to find column indices
            header_parts = first_line.lower().split(',' if has_comma else '\t')
            for idx, col_name in enumerate(header_parts):
                col_name = col_name.strip()
                if col_name in ['smiles', 'smile', 'smi', 'structure']:
                    smiles_col_idx = idx
                elif col_name in ['id', 'name', 'mol_id', 'compound_id', 'molecule_id']:
                    id_col_idx = idx
        else:
            # No header - need to detect from first data line
            # Try parsing both columns as SMILES to see which one is valid
            parts = first_line.split(',' if has_comma else '\t')
            if len(parts) >= 2:
                from rdkit import Chem
                mol_0 = Chem.MolFromSmiles(parts[0].strip())
                mol_1 = Chem.MolFromSmiles(parts[1].strip())

                # Determine which column is SMILES based on which parses successfully
                if mol_0 is not None and mol_1 is None:
                    smiles_col_idx = 0
                    id_col_idx = 1
                elif mol_1 is not None and mol_0 is None:
                    smiles_col_idx = 1
                    id_col_idx = 0
                else:
                    # Both parse or both fail - default to SMILES first
                    smiles_col_idx = 0
                    id_col_idx = 1

    # Parse the file based on detected format
    line_count = 0
    with open(smiles_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue

            line_count += 1

            # Skip header line if detected
            if line_count == 1 and is_header:
                continue

            # Split by comma if CSV, otherwise by tab/whitespace
            if has_comma:
                parts = line.split(',')
            else:
                parts = line.split('\t')

            # Handle different formats
            if len(parts) >= 2:
                smi = parts[smiles_col_idx].strip() if smiles_col_idx < len(parts) else ''
                mol_id = parts[id_col_idx].strip() if id_col_idx < len(parts) else ''
            else:
                smi = parts[0].strip() if parts else ''
                mol_id = ''

            # Skip empty, NaN, or None values
            if not smi or smi.lower() in ['nan', 'none', '', 'null']:
                continue
            if not mol_id or mol_id.lower() in ['nan', 'none', '', 'null']:
                mol_id = f'MOL_{line_count:06d}'

            smiles_data.append((smi, mol_id))

    # Prepare statistics
    stats = {
        'total': len(smiles_data),
        'success': 0,
        'failed': 0,
        'output_files': [],
        'failed_smiles': [],
        'time_elapsed': 0.0
    }

    # Worker function
    def prepare_worker(args):
        smi, mol_id = args
        output_path = os.path.join(output_dir, f'{mol_id}.pdbqt')
        result = prepare_ligand_for_docking(smi, output_path, mol_id)
        return result

    # Process in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        futures = {executor.submit(prepare_worker, args): args for args in smiles_data}

        # Process results with progress bar
        iterator = as_completed(futures)
        if verbose:
            iterator = tqdm(iterator, total=len(futures), desc="Preparing ligands")

        for future in iterator:
            result = future.result()

            if result['success']:
                stats['success'] += 1
                stats['output_files'].append(result['pdbqt_path'])
            else:
                stats['failed'] += 1
                stats['failed_smiles'].append(result['smiles'])

    stats['time_elapsed'] = time.time() - start_time

    if verbose:
        print(f"\nLigand Preparation Complete:")
        print(f"  Total: {stats['total']}")
        print(f"  Success: {stats['success']}")
        print(f"  Failed: {stats['failed']}")
        print(f"  Time: {stats['time_elapsed']:.1f} seconds")

    return stats
