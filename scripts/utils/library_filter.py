"""
Library filtering module for Digital CRO.

This module handles efficient filtering of large compound libraries,
including PAINS detection, reactive group filtering, and diversity selection.
"""

import time
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from pathlib import Path
from tqdm import tqdm

from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, FilterCatalog, AllChem
from rdkit.Chem.FilterCatalog import FilterCatalogParams

# Import from our properties module
import sys
from pathlib import Path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from scripts.utils.properties import passes_lipinski, calculate_lipinski_properties


def filter_large_library(
    input_file: str,
    output_file: str,
    filters: Optional[Dict] = None,
    chunk_size: int = 10000,
    verbose: bool = True
) -> Dict:
    """
    Filter large SMILES library efficiently with chunked processing.

    Handles files with millions of molecules by processing in chunks.

    Args:
        input_file (str): Path to input SMILES file (.smi, .csv, .txt)
        output_file (str): Path to save filtered molecules
        filters (Dict, optional): Filter criteria (see defaults below)
        chunk_size (int): Process this many molecules at once. Default: 10000
        verbose (bool): Print progress. Default: True

    Default filters:
        {
            'lipinski': True,           # Lipinski's Rule of Five
            'mw_range': (150, 500),     # Molecular weight
            'logp_range': (-0.4, 5.6),  # Lipophilicity
            'hbd_max': 5,               # H-bond donors
            'hba_max': 10,              # H-bond acceptors
            'rotatable_max': 10,        # Rotatable bonds
            'remove_pains': True,       # Remove PAINS compounds
            'remove_reactive': True,    # Remove reactive groups
            'unique_only': True         # Remove duplicates
        }

    Returns:
        dict: Statistics including:
            - input_count (int): Total input molecules
            - output_count (int): Molecules passing filters
            - filter_pass_rate (float): Percentage passing
            - filters_applied (dict): Count rejected by each filter
            - time_elapsed (float): Processing time in seconds

    Examples:
        >>> stats = filter_large_library('library.smi', 'filtered.smi')
        >>> print(f"Pass rate: {stats['filter_pass_rate']:.1f}%")
    """
    start_time = time.time()

    # Set default filters
    if filters is None:
        filters = {
            'lipinski': True,
            'mw_range': (150, 500),
            'logp_range': (-0.4, 5.6),
            'hbd_max': 5,
            'hba_max': 10,
            'rotatable_max': 10,
            'remove_pains': True,
            'remove_reactive': True,
            'unique_only': True
        }

    # Initialize statistics
    stats = {
        'input_count': 0,
        'output_count': 0,
        'filter_pass_rate': 0.0,
        'filters_applied': {
            'invalid_smiles': 0,
            'lipinski_fail': 0,
            'pains': 0,
            'reactive': 0,
            'duplicates': 0,
            'custom_filters': 0
        },
        'time_elapsed': 0.0
    }

    # Create output directory
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    # Track seen molecules for duplicate removal
    seen_smiles = set() if filters.get('unique_only', True) else None

    # Open output file
    with open(output_file, 'w') as out_f:
        # Write header
        out_f.write('smiles\tmol_id\n')

        # Read and process in chunks
        try:
            # Try to read as CSV/TSV first
            chunks = pd.read_csv(
                input_file,
                sep=None,  # Auto-detect separator
                engine='python',
                header=None,
                chunksize=chunk_size,
                on_bad_lines='skip'
            )

            chunk_iterator = tqdm(chunks, desc="Processing chunks") if verbose else chunks

            for chunk_idx, chunk in enumerate(chunk_iterator):
                # Parse SMILES column (assume first column)
                if len(chunk.columns) >= 2:
                    chunk.columns = ['smiles', 'mol_id'] + [f'col_{i}' for i in range(len(chunk.columns) - 2)]
                else:
                    chunk.columns = ['smiles']
                    chunk['mol_id'] = [f'MOL_{stats["input_count"] + i:08d}' for i in range(len(chunk))]

                # Process each molecule in chunk
                for idx, row in chunk.iterrows():
                    stats['input_count'] += 1
                    smiles = row['smiles']
                    mol_id = row['mol_id']

                    # Check for duplicates
                    if seen_smiles is not None:
                        if smiles in seen_smiles:
                            stats['filters_applied']['duplicates'] += 1
                            continue
                        seen_smiles.add(smiles)

                    # Apply filters
                    passes, reason = apply_custom_filters(smiles, filters)

                    if passes:
                        stats['output_count'] += 1
                        out_f.write(f'{smiles}\t{mol_id}\n')
                    else:
                        # Track rejection reason
                        if reason in stats['filters_applied']:
                            stats['filters_applied'][reason] += 1
                        else:
                            stats['filters_applied']['custom_filters'] += 1

        except Exception as e:
            raise IOError(f"Error processing library file: {str(e)}")

    # Calculate statistics
    if stats['input_count'] > 0:
        stats['filter_pass_rate'] = (stats['output_count'] / stats['input_count']) * 100
    stats['time_elapsed'] = time.time() - start_time

    if verbose:
        print(f"\n{'='*60}")
        print(f"Library Filtering Complete")
        print(f"{'='*60}")
        print(f"Input:  {stats['input_count']:,} molecules")
        print(f"Output: {stats['output_count']:,} molecules")
        print(f"Pass rate: {stats['filter_pass_rate']:.1f}%")
        print(f"Time: {stats['time_elapsed']:.1f} seconds")
        print(f"\nRejection reasons:")
        for reason, count in stats['filters_applied'].items():
            if count > 0:
                print(f"  {reason}: {count:,}")
        print(f"{'='*60}\n")

    return stats


def apply_custom_filters(smiles: str, filters: Dict) -> Tuple[bool, str]:
    """
    Apply all filters to a single molecule.

    Args:
        smiles (str): SMILES string
        filters (Dict): Filter criteria

    Returns:
        Tuple[bool, str]: (passes, rejection_reason)
            - passes: True if molecule passes all filters
            - rejection_reason: String explaining failure, or '' if passed

    Examples:
        >>> passes, reason = apply_custom_filters('CCO', filters)
        >>> if not passes:
        ...     print(f"Rejected: {reason}")
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, 'invalid_smiles'

    # Check PAINS
    if filters.get('remove_pains', True):
        if detect_pains_substructures(mol):
            return False, 'pains'

    # Check reactive groups
    if filters.get('remove_reactive', True):
        if detect_reactive_groups(mol):
            return False, 'reactive'

    # Check Lipinski
    if filters.get('lipinski', False):
        if not passes_lipinski(mol):
            return False, 'lipinski_fail'

    # Check MW range
    if 'mw_range' in filters:
        mw = Descriptors.MolWt(mol)
        min_mw, max_mw = filters['mw_range']
        if not (min_mw <= mw <= max_mw):
            return False, 'lipinski_fail'

    # Check LogP range
    if 'logp_range' in filters:
        logp = Descriptors.MolLogP(mol)
        min_logp, max_logp = filters['logp_range']
        if not (min_logp <= logp <= max_logp):
            return False, 'lipinski_fail'

    # Check HBD
    if 'hbd_max' in filters:
        hbd = Descriptors.NumHDonors(mol)
        if hbd > filters['hbd_max']:
            return False, 'lipinski_fail'

    # Check HBA
    if 'hba_max' in filters:
        hba = Descriptors.NumHAcceptors(mol)
        if hba > filters['hba_max']:
            return False, 'lipinski_fail'

    # Check rotatable bonds
    if 'rotatable_max' in filters:
        rotatable = Descriptors.NumRotatableBonds(mol)
        if rotatable > filters['rotatable_max']:
            return False, 'lipinski_fail'

    return True, ''


def detect_pains_substructures(mol: Chem.Mol) -> bool:
    """
    Check for PAINS (Pan-Assay Interference) patterns.

    PAINS are molecules that give false positives in biological assays.
    Uses RDKit's built-in PAINS catalog.

    Args:
        mol (Chem.Mol): RDKit molecule

    Returns:
        bool: True if contains PAINS, False if clean

    Common PAINS patterns:
        - Catechols
        - Rhodanines
        - Quinones
        - Alkyl halides
        - Hydroxyphenyl hydrazones

    Examples:
        >>> mol = Chem.MolFromSmiles('c1ccc(O)c(O)c1')  # Catechol
        >>> detect_pains_substructures(mol)
        True
    """
    if mol is None:
        return True

    try:
        # Create PAINS filter catalog
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)

        # Check for matches
        matches = catalog.GetMatches(mol)

        return len(matches) > 0

    except Exception:
        # If filtering fails, be conservative and reject
        return True


def detect_reactive_groups(mol: Chem.Mol) -> bool:
    """
    Check for reactive/unstable functional groups.

    Flags molecules with groups that may be:
    - Chemically unstable
    - Reactive with proteins (non-specifically)
    - Toxic

    Args:
        mol (Chem.Mol): RDKit molecule

    Returns:
        bool: True if contains reactive groups, False if stable

    Patterns flagged:
        - Acid chlorides: C(=O)Cl
        - Isocyanates: N=C=O
        - Thiols: [SH]
        - Peroxides: O-O
        - Diazo groups: N=[N+]=[N-]
        - Aldehydes (some)
        - Michael acceptors (some)

    Examples:
        >>> mol = Chem.MolFromSmiles('CC(=O)Cl')  # Acetyl chloride
        >>> detect_reactive_groups(mol)
        True
    """
    if mol is None:
        return True

    # Define reactive SMARTS patterns
    reactive_smarts = [
        'C(=O)Cl',           # Acid chloride
        'S(=O)(=O)Cl',       # Sulfonyl chloride
        'N=C=O',             # Isocyanate
        'N=C=S',             # Isothiocyanate
        'C(=O)O(C=O)',       # Anhydride
        '[N,O,S]=[N+]=[N-]', # Diazo
        'O-O',               # Peroxide
        'C#N',               # Nitrile (sometimes reactive)
        '[SH]',              # Thiol (can be reactive)
        'C=[N+]=[N-]',       # Diazonium
        '[F,Cl,Br,I][C,c]',  # Halides (some are reactive)
        'C(F)(F)F',          # Trifluoromethyl (sometimes toxic)
    ]

    # Check each pattern
    for smarts in reactive_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            return True

    return False


def calculate_diversity_subset(
    smiles_list: List[str],
    target_count: int,
    method: str = 'maxmin',
    fingerprint_radius: int = 2,
    fingerprint_bits: int = 2048
) -> List[str]:
    """
    Select diverse subset from large library using MaxMin algorithm.

    Args:
        smiles_list (List[str]): Full list of SMILES strings
        target_count (int): Number of molecules to select
        method (str): Selection method. Currently only 'maxmin' supported.
        fingerprint_radius (int): Morgan fingerprint radius. Default: 2
        fingerprint_bits (int): Fingerprint bit length. Default: 2048

    Algorithm (MaxMin):
        1. Convert all SMILES to Morgan fingerprints
        2. Pick random starting molecule
        3. Iteratively select molecule most dissimilar to already selected
        4. Continue until target_count reached

    Uses Tanimoto similarity on Morgan fingerprints for diversity measurement.

    Returns:
        List[str]: Selected SMILES (diverse subset)

    Examples:
        >>> diverse = calculate_diversity_subset(smiles_list, 100)
        >>> print(f"Selected {len(diverse)} diverse molecules")
    """
    if len(smiles_list) <= target_count:
        return smiles_list

    # Convert SMILES to molecules and fingerprints
    mols = []
    fps = []
    valid_smiles = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(
                mol,
                fingerprint_radius,
                nBits=fingerprint_bits
            )
            mols.append(mol)
            fps.append(fp)
            valid_smiles.append(smiles)

    if len(fps) <= target_count:
        return valid_smiles

    # MaxMin algorithm
    selected_indices = []
    remaining_indices = list(range(len(fps)))

    # Start with random molecule
    first_idx = np.random.choice(remaining_indices)
    selected_indices.append(first_idx)
    remaining_indices.remove(first_idx)

    # Iteratively select most dissimilar molecule
    for _ in range(target_count - 1):
        if not remaining_indices:
            break

        max_min_sim = -1
        best_idx = None

        for idx in remaining_indices:
            # Calculate minimum similarity to selected molecules
            min_sim = min([
                DataStructs.TanimotoSimilarity(fps[idx], fps[sel_idx])
                for sel_idx in selected_indices
            ])

            # Keep track of molecule with maximum minimum similarity
            if min_sim > max_min_sim:
                max_min_sim = min_sim
                best_idx = idx

        if best_idx is not None:
            selected_indices.append(best_idx)
            remaining_indices.remove(best_idx)

    # Return selected SMILES
    return [valid_smiles[i] for i in selected_indices]


def calculate_library_diversity(smiles_list: List[str]) -> Dict[str, float]:
    """
    Calculate diversity metrics for a compound library.

    Args:
        smiles_list (List[str]): List of SMILES strings

    Returns:
        dict: Diversity metrics:
            - avg_tanimoto (float): Average pairwise Tanimoto similarity
            - min_tanimoto (float): Minimum similarity
            - max_tanimoto (float): Maximum similarity
            - diversity_score (float): 1 - avg_tanimoto (higher = more diverse)

    Examples:
        >>> metrics = calculate_library_diversity(smiles_list)
        >>> print(f"Diversity: {metrics['diversity_score']:.2f}")
    """
    # Generate fingerprints
    fps = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fps.append(fp)

    if len(fps) < 2:
        return {
            'avg_tanimoto': 0.0,
            'min_tanimoto': 0.0,
            'max_tanimoto': 0.0,
            'diversity_score': 1.0
        }

    # Calculate pairwise similarities
    similarities = []
    for i in range(len(fps)):
        for j in range(i + 1, len(fps)):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            similarities.append(sim)

    avg_sim = np.mean(similarities)
    min_sim = np.min(similarities)
    max_sim = np.max(similarities)

    return {
        'avg_tanimoto': avg_sim,
        'min_tanimoto': min_sim,
        'max_tanimoto': max_sim,
        'diversity_score': 1.0 - avg_sim
    }


def save_filtered_library(
    smiles_list: List[str],
    output_file: str,
    include_properties: bool = False
) -> None:
    """
    Save filtered library to file with optional property calculations.

    Args:
        smiles_list (List[str]): Filtered SMILES
        output_file (str): Output file path
        include_properties (bool): Calculate and save Lipinski properties

    Examples:
        >>> save_filtered_library(filtered, 'output.smi', include_properties=True)
    """
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    if include_properties:
        # Calculate properties for each molecule
        data = []
        for idx, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                props = calculate_lipinski_properties(mol)
                props['smiles'] = smiles
                props['mol_id'] = f'MOL_{idx:08d}'
                data.append(props)

        df = pd.DataFrame(data)
        df.to_csv(output_file, index=False, sep='\t')
    else:
        # Simple SMILES file
        with open(output_file, 'w') as f:
            f.write('smiles\tmol_id\n')
            for idx, smiles in enumerate(smiles_list):
                f.write(f'{smiles}\tMOL_{idx:08d}\n')
