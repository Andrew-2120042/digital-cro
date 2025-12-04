"""
Molecular properties calculation module for Digital CRO.

This module provides functions for calculating molecular descriptors and
evaluating drug-likeness based on Lipinski's Rule of Five.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from typing import Dict, Optional


def calculate_lipinski_properties(mol: Chem.Mol) -> Dict[str, float]:
    """
    Calculate Lipinski Rule of Five properties.

    Properties calculated:
    - Molecular weight (MW)
    - LogP (octanol-water partition coefficient)
    - Number of hydrogen bond donors (HBD)
    - Number of hydrogen bond acceptors (HBA)
    - Number of rotatable bonds

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        dict: Dictionary containing property names and values

    Raises:
        ValueError: If mol is None or invalid

    Examples:
        >>> mol = Chem.MolFromSmiles('CC(=O)O')
        >>> props = calculate_lipinski_properties(mol)
        >>> print(props['molecular_weight'])
    """
    if mol is None:
        raise ValueError("Cannot calculate properties for None molecule")

    try:
        properties = {
            'molecular_weight': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'hbd': Descriptors.NumHDonors(mol),
            'hba': Descriptors.NumHAcceptors(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol)
        }
        return properties

    except Exception as e:
        raise ValueError(f"Error calculating Lipinski properties: {str(e)}")


def passes_lipinski(mol: Chem.Mol) -> bool:
    """
    Check if molecule passes Lipinski's Rule of Five.

    Lipinski criteria:
    - Molecular weight: 150-500 Da
    - LogP: -0.4 to 5.6
    - Hydrogen bond donors (HBD): ≤5
    - Hydrogen bond acceptors (HBA): ≤10

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        bool: True if molecule passes all Lipinski criteria, False otherwise

    Raises:
        ValueError: If mol is None or invalid

    Examples:
        >>> mol = Chem.MolFromSmiles('CC(=O)O')
        >>> if passes_lipinski(mol):
        ...     print("Drug-like molecule")
    """
    if mol is None:
        raise ValueError("Cannot evaluate Lipinski rules for None molecule")

    try:
        props = calculate_lipinski_properties(mol)

        # Check all Lipinski criteria
        mw_ok = 150 <= props['molecular_weight'] <= 500
        logp_ok = -0.4 <= props['logp'] <= 5.6
        hbd_ok = props['hbd'] <= 5
        hba_ok = props['hba'] <= 10

        return mw_ok and logp_ok and hbd_ok and hba_ok

    except Exception as e:
        raise ValueError(f"Error evaluating Lipinski rules: {str(e)}")


def calculate_all_properties(mol: Chem.Mol) -> Dict[str, any]:
    """
    Calculate comprehensive set of molecular properties.

    Properties calculated:
    - All Lipinski properties (MW, LogP, HBD, HBA, rotatable bonds)
    - Molecular formula
    - Number of rings
    - Number of aromatic rings
    - TPSA (topological polar surface area)
    - Number of atoms
    - Number of heavy atoms
    - Fraction of sp3 carbons
    - Number of stereocenters

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        dict: Dictionary containing all calculated properties

    Raises:
        ValueError: If mol is None or invalid

    Examples:
        >>> mol = Chem.MolFromSmiles('CC(=O)Nc1ccc(O)cc1')
        >>> props = calculate_all_properties(mol)
        >>> print(f"TPSA: {props['tpsa']}")
    """
    if mol is None:
        raise ValueError("Cannot calculate properties for None molecule")

    try:
        # Get Lipinski properties
        lipinski_props = calculate_lipinski_properties(mol)

        # Calculate additional properties
        properties = {
            **lipinski_props,  # Include all Lipinski properties
            'molecular_formula': rdMolDescriptors.CalcMolFormula(mol),
            'num_rings': Descriptors.RingCount(mol),
            'num_aromatic_rings': Descriptors.NumAromaticRings(mol),
            'tpsa': Descriptors.TPSA(mol),
            'num_atoms': mol.GetNumAtoms(),
            'num_heavy_atoms': Lipinski.HeavyAtomCount(mol),
            'fraction_csp3': Descriptors.FractionCSP3(mol),
            'num_stereocenters': len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
            'passes_lipinski': passes_lipinski(mol)
        }

        return properties

    except Exception as e:
        raise ValueError(f"Error calculating all properties: {str(e)}")


def get_property_summary(properties_list: list) -> Dict[str, Dict[str, float]]:
    """
    Calculate summary statistics for a list of property dictionaries.

    Args:
        properties_list (list): List of property dictionaries from calculate_all_properties

    Returns:
        dict: Dictionary with statistics (mean, min, max, std) for each numeric property

    Examples:
        >>> props_list = [calculate_all_properties(mol) for mol in mols]
        >>> summary = get_property_summary(props_list)
        >>> print(summary['molecular_weight']['mean'])
    """
    import numpy as np

    if not properties_list:
        return {}

    # Get numeric properties
    numeric_props = {}
    for key in properties_list[0].keys():
        if isinstance(properties_list[0][key], (int, float)):
            numeric_props[key] = [p[key] for p in properties_list if key in p]

    # Calculate statistics
    summary = {}
    for prop_name, values in numeric_props.items():
        if values:
            summary[prop_name] = {
                'mean': np.mean(values),
                'min': np.min(values),
                'max': np.max(values),
                'std': np.std(values),
                'median': np.median(values)
            }

    return summary
