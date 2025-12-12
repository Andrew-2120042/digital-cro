"""
Molecule operations module for Digital CRO.

This module provides functions for loading, converting, and visualizing molecules.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from typing import Optional, List, Union
import os


def load_smiles_file(filepath: str) -> pd.DataFrame:
    """
    Load SMILES file into DataFrame.

    Expected format: SMILES, ID, optional properties (comma or tab separated)

    Args:
        filepath (str): Path to the SMILES file

    Returns:
        pd.DataFrame: DataFrame with columns ['smiles', 'mol_id', ...]

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid

    Examples:
        >>> df = load_smiles_file('data/molecules/compounds.smi')
        >>> print(df.head())
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    try:
        # Try comma-separated first
        df = pd.read_csv(filepath, sep=',', header=None)

        # If only one column, try tab-separated
        if len(df.columns) == 1:
            df = pd.read_csv(filepath, sep='\t', header=None)

        # If only one column, try space-separated
        if len(df.columns) == 1:
            df = pd.read_csv(filepath, sep=' ', header=None)

        # Rename columns
        if len(df.columns) >= 2:
            df.columns = ['smiles', 'mol_id'] + [f'prop_{i}' for i in range(len(df.columns) - 2)]
        elif len(df.columns) == 1:
            df.columns = ['smiles']
            df['mol_id'] = [f'MOL_{i:04d}' for i in range(len(df))]
        else:
            raise ValueError("Invalid file format: no columns detected")

        # Validate SMILES column exists and has data
        if df.empty:
            raise ValueError("File is empty")

        if df['smiles'].isna().all():
            raise ValueError("No valid SMILES strings found")

        return df

    except Exception as e:
        raise ValueError(f"Error parsing SMILES file: {str(e)}")


def smiles_to_mol(smiles: str) -> Optional[Chem.Mol]:
    """
    Convert SMILES string to RDKit molecule object.

    Handle errors gracefully by returning None if invalid.

    Args:
        smiles (str): SMILES string representation of molecule

    Returns:
        Chem.Mol: RDKit molecule object, or None if conversion fails

    Examples:
        >>> mol = smiles_to_mol('CC(=O)O')  # Acetic acid
        >>> if mol is not None:
        ...     print("Valid molecule")
    """
    if not smiles or not isinstance(smiles, str):
        return None

    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        return mol
    except Exception:
        return None


def draw_molecule(mol: Chem.Mol, size: tuple = (300, 300), filename: Optional[str] = None) -> Optional[object]:
    """
    Draw 2D structure of molecule.

    If filename provided, save as PNG. Otherwise return image object.

    Args:
        mol (Chem.Mol): RDKit molecule object
        size (tuple): Image size as (width, height). Default: (300, 300)
        filename (str, optional): Path to save PNG file. If None, returns image object

    Returns:
        PIL.Image or None: Image object if no filename provided, otherwise None

    Raises:
        ValueError: If mol is None or invalid

    Examples:
        >>> mol = smiles_to_mol('CC(=O)O')
        >>> draw_molecule(mol, filename='data/outputs/acetic_acid.png')
        >>> img = draw_molecule(mol)  # Returns image object
    """
    if mol is None:
        raise ValueError("Cannot draw None molecule")

    try:
        img = Draw.MolToImage(mol, size=size)

        if filename:
            # Ensure directory exists
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            img.save(filename)
            return None
        else:
            return img

    except Exception as e:
        raise ValueError(f"Error drawing molecule: {str(e)}")


def draw_molecule_grid(
    mols: List[Chem.Mol],
    labels: Optional[List[str]] = None,
    mols_per_row: int = 4,
    filename: Optional[str] = None,
    sub_img_size: tuple = (300, 300)
) -> Optional[object]:
    """
    Draw grid of multiple molecules.

    Args:
        mols (List[Chem.Mol]): List of RDKit molecule objects
        labels (List[str], optional): Labels for each molecule. If None, no labels shown
        mols_per_row (int): Number of molecules per row. Default: 4
        filename (str, optional): Path to save PNG file. If None, returns image object
        sub_img_size (tuple): Size of each molecule image. Default: (300, 300)

    Returns:
        PIL.Image or None: Image object if no filename provided, otherwise None

    Raises:
        ValueError: If mols list is empty or contains invalid molecules

    Examples:
        >>> mols = [smiles_to_mol(s) for s in ['CCO', 'CC(=O)O', 'c1ccccc1']]
        >>> labels = ['Ethanol', 'Acetic Acid', 'Benzene']
        >>> draw_molecule_grid(mols, labels=labels, filename='grid.png')
    """
    if not mols:
        raise ValueError("Cannot draw grid: molecule list is empty")

    # Filter out None molecules
    valid_mols = [m for m in mols if m is not None]

    if not valid_mols:
        raise ValueError("Cannot draw grid: no valid molecules")

    # Adjust labels if needed
    if labels:
        if len(labels) != len(mols):
            raise ValueError(f"Number of labels ({len(labels)}) must match number of molecules ({len(mols)})")
        # Filter labels for valid molecules only
        valid_labels = [labels[i] for i, m in enumerate(mols) if m is not None]
    else:
        valid_labels = None

    try:
        img = Draw.MolsToGridImage(
            valid_mols,
            molsPerRow=mols_per_row,
            subImgSize=sub_img_size,
            legends=valid_labels
        )

        if filename:
            # Ensure directory exists
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            img.save(filename)
            return None
        else:
            return img

    except Exception as e:
        raise ValueError(f"Error drawing molecule grid: {str(e)}")
