"""
ADMET Predictions Module for Digital CRO.

This module predicts drug-likeness and ADMET properties using RDKit.

Key predictions:
    - Lipinski Rule of Five compliance
    - TPSA (Topological Polar Surface Area)
    - Synthetic Accessibility Score
    - QED (Quantitative Estimate of Drug-likeness)
    - Blood-Brain Barrier (BBB) penetration
    - Oral bioavailability estimation

References:
    - Lipinski et al. (2001) - Rule of Five
    - Ertl & Schuffenhauer (2009) - Synthetic Accessibility
    - Bickerton et al. (2012) - QED Score
"""

import warnings
from typing import Dict, Optional, List, Tuple
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.QED import qed


def calculate_lipinski_properties(smiles: str) -> Dict:
    """
    Calculate Lipinski Rule of Five properties.

    Lipinski's Rule of Five criteria for oral drug-likeness:
    - Molecular weight ≤ 500 Da
    - LogP ≤ 5
    - H-bond donors ≤ 5
    - H-bond acceptors ≤ 10

    Args:
        smiles (str): SMILES string of molecule

    Returns:
        dict: Dictionary containing:
            - molecular_weight (float): MW in Da
            - logp (float): Partition coefficient
            - hbd (int): H-bond donors
            - hba (int): H-bond acceptors
            - num_violations (int): Number of Lipinski violations (0-4)
            - lipinski_compliant (bool): True if ≤1 violation

    Examples:
        >>> props = calculate_lipinski_properties('CC(=O)Oc1ccccc1C(=O)O')  # Aspirin
        >>> print(f"MW: {props['molecular_weight']:.1f}, LogP: {props['logp']:.2f}")
        >>> print(f"Lipinski compliant: {props['lipinski_compliant']}")
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return {
            'molecular_weight': None,
            'logp': None,
            'hbd': None,
            'hba': None,
            'num_violations': None,
            'lipinski_compliant': False,
            'error': 'Invalid SMILES'
        }

    # Calculate properties
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    # Count violations
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1

    return {
        'molecular_weight': mw,
        'logp': logp,
        'hbd': hbd,
        'hba': hba,
        'num_violations': violations,
        'lipinski_compliant': violations <= 1  # Allow 1 violation
    }


def calculate_sa_score(mol: Chem.Mol) -> float:
    """
    Calculate Synthetic Accessibility Score.

    Score ranges from 1 (easy to synthesize) to 10 (difficult to synthesize).

    This is a simplified implementation based on:
    - Molecular complexity (number of rings, stereocenters)
    - Presence of uncommon substructures
    - Size and rotatable bonds

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        float: SA score (1-10, lower is better)

    Note:
        This is a simplified heuristic. For production, use the full
        Ertl & Schuffenhauer implementation from rdkit.Chem.SA_Score.

    Examples:
        >>> mol = Chem.MolFromSmiles('CCO')  # Ethanol
        >>> sa = calculate_sa_score(mol)
        >>> print(f"SA Score: {sa:.2f} (1=easy, 10=hard)")
    """
    # Simplified SA score based on molecular complexity

    # Factor 1: Size penalty
    num_atoms = mol.GetNumHeavyAtoms()
    size_penalty = min(num_atoms / 50.0, 1.0) * 3  # Max 3 points

    # Factor 2: Ring complexity
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    ring_penalty = min(num_rings / 5.0, 1.0) * 2  # Max 2 points

    # Factor 3: Stereo complexity
    num_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    stereo_penalty = min(num_stereo / 5.0, 1.0) * 2  # Max 2 points

    # Factor 4: Rotatable bonds (flexibility)
    num_rotatable = Lipinski.NumRotatableBonds(mol)
    rot_penalty = min(num_rotatable / 10.0, 1.0) * 1.5  # Max 1.5 points

    # Factor 5: Heteroatom complexity
    num_hetero = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
    hetero_penalty = min(num_hetero / 10.0, 1.0) * 1.5  # Max 1.5 points

    # Base score of 1 (simplest molecules)
    sa_score = 1.0 + size_penalty + ring_penalty + stereo_penalty + rot_penalty + hetero_penalty

    return min(sa_score, 10.0)


def calculate_qed(mol: Chem.Mol) -> float:
    """
    Calculate Quantitative Estimate of Drug-likeness (QED).

    QED is a composite score (0-1) based on 8 molecular properties:
    - Molecular weight
    - LogP
    - H-bond donors/acceptors
    - Polar surface area
    - Rotatable bonds
    - Aromatic rings
    - Structural alerts

    Score interpretation:
    - 0.0-0.3: Low drug-likeness
    - 0.3-0.5: Moderate drug-likeness
    - 0.5-0.7: Good drug-likeness
    - 0.7-1.0: Excellent drug-likeness

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        float: QED score (0-1, higher is better)

    Examples:
        >>> mol = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')  # Caffeine
        >>> qed_score = calculate_qed(mol)
        >>> print(f"QED: {qed_score:.3f}")
    """
    return qed(mol)


def predict_bbb_penetration(mol: Chem.Mol) -> Dict:
    """
    Predict Blood-Brain Barrier (BBB) penetration.

    Rules for BBB penetration (simplified):
    - TPSA < 90 Ų
    - Molecular weight < 450 Da
    - 1 < LogP < 3 (optimal range)
    - H-bond donors ≤ 3
    - H-bond acceptors ≤ 8

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        dict: Dictionary containing:
            - bbb_penetrant (bool): Likely to cross BBB
            - tpsa (float): Topological polar surface area
            - confidence (str): 'High', 'Medium', 'Low'
            - reasons (List[str]): Factors affecting prediction

    Examples:
        >>> mol = Chem.MolFromSmiles('CN1CCC23c4cccc5C6=C(C=Cc4C1Cc2c5O)O3')  # Morphine
        >>> bbb = predict_bbb_penetration(mol)
        >>> print(f"BBB penetrant: {bbb['bbb_penetrant']}")
    """
    # Calculate relevant properties
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    # Evaluate criteria
    criteria_met = []
    criteria_failed = []

    if tpsa < 90:
        criteria_met.append('TPSA < 90 Ų')
    else:
        criteria_failed.append(f'TPSA too high ({tpsa:.1f} Ų)')

    if mw < 450:
        criteria_met.append('MW < 450 Da')
    else:
        criteria_failed.append(f'MW too high ({mw:.1f} Da)')

    if 1 < logp < 3:
        criteria_met.append('LogP in optimal range (1-3)')
    elif logp <= 1:
        criteria_failed.append(f'LogP too low ({logp:.2f})')
    else:
        criteria_failed.append(f'LogP too high ({logp:.2f})')

    if hbd <= 3:
        criteria_met.append('HBD ≤ 3')
    else:
        criteria_failed.append(f'Too many H-bond donors ({hbd})')

    if hba <= 8:
        criteria_met.append('HBA ≤ 8')
    else:
        criteria_failed.append(f'Too many H-bond acceptors ({hba})')

    # Decision
    num_met = len(criteria_met)
    bbb_penetrant = num_met >= 4  # At least 4 of 5 criteria

    if num_met == 5:
        confidence = 'High'
    elif num_met >= 3:
        confidence = 'Medium'
    else:
        confidence = 'Low'

    return {
        'bbb_penetrant': bbb_penetrant,
        'tpsa': tpsa,
        'confidence': confidence,
        'criteria_met': criteria_met,
        'criteria_failed': criteria_failed
    }


def predict_oral_bioavailability(mol: Chem.Mol) -> Dict:
    """
    Predict oral bioavailability category.

    Combines Lipinski Rule of Five with Veber criteria:
    - Rotatable bonds ≤ 10
    - TPSA ≤ 140 Ų

    Categories:
    - High: Lipinski compliant + Veber compliant
    - Medium: 1-2 violations
    - Low: 3+ violations

    Args:
        mol (Chem.Mol): RDKit molecule object

    Returns:
        dict: Dictionary containing:
            - bioavailability (str): 'High', 'Medium', 'Low'
            - lipinski_violations (int): Number of Lipinski violations
            - veber_violations (int): Number of Veber violations
            - total_violations (int): Combined violations

    Examples:
        >>> mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')  # Aspirin
        >>> bioav = predict_oral_bioavailability(mol)
        >>> print(f"Oral bioavailability: {bioav['bioavailability']}")
    """
    # Lipinski properties
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    lipinski_violations = 0
    if mw > 500:
        lipinski_violations += 1
    if logp > 5:
        lipinski_violations += 1
    if hbd > 5:
        lipinski_violations += 1
    if hba > 10:
        lipinski_violations += 1

    # Veber criteria
    num_rotatable = Lipinski.NumRotatableBonds(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)

    veber_violations = 0
    if num_rotatable > 10:
        veber_violations += 1
    if tpsa > 140:
        veber_violations += 1

    total_violations = lipinski_violations + veber_violations

    # Category
    if total_violations == 0:
        bioavailability = 'High'
    elif total_violations <= 2:
        bioavailability = 'Medium'
    else:
        bioavailability = 'Low'

    return {
        'bioavailability': bioavailability,
        'lipinski_violations': lipinski_violations,
        'veber_violations': veber_violations,
        'total_violations': total_violations,
        'rotatable_bonds': num_rotatable,
        'tpsa': tpsa
    }


def predict_admet_properties(smiles: str) -> Dict:
    """
    Comprehensive ADMET property prediction.

    Args:
        smiles (str): SMILES string of molecule

    Returns:
        dict: Dictionary containing:
            - smiles (str): Input SMILES
            - molecular_weight (float): MW in Da
            - logp (float): Partition coefficient
            - tpsa (float): Topological polar surface area
            - hbd (int): H-bond donors
            - hba (int): H-bond acceptors
            - rotatable_bonds (int): Number of rotatable bonds
            - aromatic_rings (int): Number of aromatic rings
            - qed_score (float): Drug-likeness (0-1)
            - sa_score (float): Synthetic accessibility (1-10)
            - lipinski_compliant (bool): Passes Rule of Five
            - bbb_penetrant (bool): BBB penetration predicted
            - oral_bioavailability (str): 'High', 'Medium', 'Low'
            - error (str): Error message if failed

    Examples:
        >>> admet = predict_admet_properties('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
        >>> print(f"QED: {admet['qed_score']:.3f}")
        >>> print(f"Bioavailability: {admet['oral_bioavailability']}")
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return {
            'smiles': smiles,
            'error': 'Invalid SMILES'
        }

    # Basic properties
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    num_rotatable = Lipinski.NumRotatableBonds(mol)
    num_aromatic_rings = Lipinski.NumAromaticRings(mol)

    # Lipinski compliance
    lipinski = calculate_lipinski_properties(smiles)

    # Drug-likeness scores
    qed_score = calculate_qed(mol)
    sa_score = calculate_sa_score(mol)

    # BBB penetration
    bbb = predict_bbb_penetration(mol)

    # Oral bioavailability
    bioav = predict_oral_bioavailability(mol)

    return {
        'smiles': smiles,
        'molecular_weight': mw,
        'logp': logp,
        'tpsa': tpsa,
        'hbd': hbd,
        'hba': hba,
        'rotatable_bonds': num_rotatable,
        'aromatic_rings': num_aromatic_rings,
        'qed_score': qed_score,
        'sa_score': sa_score,
        'lipinski_compliant': lipinski['lipinski_compliant'],
        'lipinski_violations': lipinski['num_violations'],
        'bbb_penetrant': bbb['bbb_penetrant'],
        'bbb_confidence': bbb['confidence'],
        'oral_bioavailability': bioav['bioavailability'],
        'error': None
    }


def analyze_molecule_complete(smiles: str, mol_id: Optional[str] = None) -> Dict:
    """
    Complete molecular analysis including ADMET predictions.

    Args:
        smiles (str): SMILES string
        mol_id (str, optional): Molecule identifier

    Returns:
        dict: Complete analysis including:
            - Basic properties (MW, LogP, TPSA, etc.)
            - Lipinski compliance
            - Drug-likeness (QED, SA score)
            - BBB penetration prediction
            - Oral bioavailability
            - Detailed violation reports

    Examples:
        >>> analysis = analyze_molecule_complete('CC(=O)Oc1ccccc1C(=O)O', 'aspirin')
        >>> print(f"Molecule: {analysis['mol_id']}")
        >>> print(f"QED: {analysis['qed_score']:.3f}")
        >>> print(f"Bioavailability: {analysis['oral_bioavailability']}")
    """
    result = predict_admet_properties(smiles)

    if mol_id:
        result['mol_id'] = mol_id

    # Add detailed BBB analysis
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        bbb_details = predict_bbb_penetration(mol)
        result['bbb_details'] = {
            'penetrant': bbb_details['bbb_penetrant'],
            'confidence': bbb_details['confidence'],
            'criteria_met': bbb_details['criteria_met'],
            'criteria_failed': bbb_details['criteria_failed']
        }

        # Add bioavailability details
        bioav_details = predict_oral_bioavailability(mol)
        result['bioavailability_details'] = {
            'category': bioav_details['bioavailability'],
            'lipinski_violations': bioav_details['lipinski_violations'],
            'veber_violations': bioav_details['veber_violations'],
            'total_violations': bioav_details['total_violations']
        }

    return result


def batch_admet_analysis(
    df: pd.DataFrame,
    smiles_column: str = 'smiles',
    id_column: Optional[str] = None
) -> pd.DataFrame:
    """
    Batch ADMET analysis on DataFrame of molecules.

    Args:
        df (pd.DataFrame): DataFrame containing SMILES
        smiles_column (str): Name of SMILES column. Default: 'smiles'
        id_column (str, optional): Name of ID column. If None, use DataFrame index.

    Returns:
        pd.DataFrame: Input DataFrame with added ADMET property columns:
            - molecular_weight
            - logp
            - tpsa
            - qed_score
            - sa_score
            - lipinski_compliant
            - bbb_penetrant
            - oral_bioavailability

    Examples:
        >>> df = pd.DataFrame({
        ...     'smiles': ['CCO', 'CC(=O)Oc1ccccc1C(=O)O'],
        ...     'name': ['ethanol', 'aspirin']
        ... })
        >>> df_admet = batch_admet_analysis(df, id_column='name')
        >>> print(df_admet[['name', 'qed_score', 'oral_bioavailability']])
    """
    # Validate input
    if smiles_column not in df.columns:
        raise ValueError(f"Column '{smiles_column}' not found in DataFrame")

    # Initialize result columns
    admet_results = []

    for idx, row in df.iterrows():
        smiles = row[smiles_column]
        mol_id = row[id_column] if id_column and id_column in df.columns else str(idx)

        # Get ADMET predictions
        admet = predict_admet_properties(smiles)
        admet_results.append(admet)

    # Convert to DataFrame
    df_admet = pd.DataFrame(admet_results)

    # Merge with original DataFrame (keeping all original columns)
    # Drop 'smiles' from admet results to avoid duplication
    if 'smiles' in df_admet.columns:
        df_admet = df_admet.drop('smiles', axis=1)

    # Concatenate horizontally
    df_result = pd.concat([df.reset_index(drop=True), df_admet.reset_index(drop=True)], axis=1)

    return df_result
