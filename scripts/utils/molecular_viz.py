"""
Molecular Visualization Module for Digital CRO.

This module provides visualization capabilities for:
    - 2D structure drawing (grids, labels, annotations)
    - 3D conformer generation
    - Docking result visualization
    - ADMET radar plots
    - Property distribution plots

Uses RDKit for molecular drawing and matplotlib for plotting.
"""

import os
from pathlib import Path
from typing import List, Optional, Dict, Tuple
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import gridspec

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem import Descriptors


def draw_molecules_grid(
    smiles_list: List[str],
    labels: Optional[List[str]] = None,
    mols_per_row: int = 5,
    mol_size: Tuple[int, int] = (300, 300),
    output_path: Optional[str] = None
) -> None:
    """
    Draw grid of 2D molecule structures.

    Args:
        smiles_list (List[str]): List of SMILES strings
        labels (List[str], optional): Labels for each molecule (displayed below structure)
        mols_per_row (int): Number of molecules per row. Default: 5
        mol_size (Tuple[int, int]): Size of each molecule image (width, height). Default: (300, 300)
        output_path (str, optional): Path to save image. If None, displays interactively.

    Returns:
        None: Displays or saves image

    Examples:
        >>> smiles = ['CCO', 'CC(=O)O', 'c1ccccc1']
        >>> labels = ['Ethanol', 'Acetic acid', 'Benzene']
        >>> draw_molecules_grid(smiles, labels, mols_per_row=3, output_path='molecules.png')
    """
    # Convert SMILES to molecules
    mols = []
    valid_labels = []

    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mols.append(mol)
            if labels and i < len(labels):
                valid_labels.append(labels[i])
            else:
                valid_labels.append(f'Mol_{i+1}')
        else:
            warnings.warn(f"Invalid SMILES at index {i}: {smi}")

    if not mols:
        raise ValueError("No valid molecules to draw")

    # Use labels if provided
    legends = valid_labels if labels else None

    # Draw grid
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=mol_size,
        legends=legends,
        returnPNG=False
    )

    if output_path:
        # Ensure directory exists
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
        img.save(output_path)
        print(f"Saved molecule grid to: {output_path}")
    else:
        # Display
        img.show()


def generate_3d_conformer(smiles: str, optimize: bool = True) -> Optional[Chem.Mol]:
    """
    Generate 3D conformer from SMILES.

    Args:
        smiles (str): SMILES string
        optimize (bool): If True, optimize geometry with MMFF force field. Default: True

    Returns:
        Chem.Mol: RDKit molecule with 3D coordinates, or None if failed

    Examples:
        >>> mol_3d = generate_3d_conformer('CC(=O)Oc1ccccc1C(=O)O')  # Aspirin
        >>> if mol_3d:
        ...     print(f"Generated 3D conformer with {mol_3d.GetNumAtoms()} atoms")
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates
    params = AllChem.ETKDGv3()
    params.randomSeed = 42  # Reproducibility

    result = AllChem.EmbedMolecule(mol, params)

    if result == -1:
        # Embedding failed, try basic method
        result = AllChem.EmbedMolecule(mol)
        if result == -1:
            return None

    # Optimize geometry
    if optimize:
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            # MMFF might fail for some molecules
            pass

    return mol


def visualize_docking_results(
    df_results: pd.DataFrame,
    top_n: int = 20,
    smiles_column: str = 'smiles',
    id_column: str = 'ligand_id',
    affinity_column: str = 'binding_affinity',
    mols_per_row: int = 5,
    output_path: Optional[str] = None
) -> None:
    """
    Visualize docking results with 2D structures and binding affinities.

    Creates a grid showing top hits with:
    - 2D molecular structure
    - Molecule ID
    - Binding affinity (kcal/mol)
    - Rank

    Args:
        df_results (pd.DataFrame): Docking results DataFrame
        top_n (int): Number of top molecules to show. Default: 20
        smiles_column (str): Name of SMILES column. Default: 'smiles'
        id_column (str): Name of ID column. Default: 'ligand_id'
        affinity_column (str): Name of affinity column. Default: 'binding_affinity'
        mols_per_row (int): Molecules per row. Default: 5
        output_path (str, optional): Path to save image

    Examples:
        >>> visualize_docking_results(
        ...     df_ranked,
        ...     top_n=10,
        ...     output_path='top_hits.png'
        ... )
    """
    # Sort by affinity and take top N
    df_sorted = df_results.sort_values(affinity_column).head(top_n)

    smiles_list = df_sorted[smiles_column].tolist()
    mol_ids = df_sorted[id_column].tolist()
    affinities = df_sorted[affinity_column].tolist()

    # Create labels with rank, ID, and affinity
    labels = []
    for i, (mol_id, affinity) in enumerate(zip(mol_ids, affinities), 1):
        label = f"#{i} {mol_id}\n{affinity:.2f} kcal/mol"
        labels.append(label)

    # Draw grid
    draw_molecules_grid(
        smiles_list,
        labels=labels,
        mols_per_row=mols_per_row,
        output_path=output_path
    )


def create_admet_radar_plot(
    df: pd.DataFrame,
    molecule_id: str,
    id_column: str = 'ligand_id',
    output_path: Optional[str] = None
) -> None:
    """
    Create radar plot showing ADMET profile for a molecule.

    Visualizes 6 normalized properties:
    - QED (drug-likeness)
    - Lipinski compliance
    - Oral bioavailability
    - BBB penetration
    - Synthetic accessibility (inverted: 10-SA so higher is better)
    - Molecular weight (normalized: 1 - MW/1000)

    Args:
        df (pd.DataFrame): DataFrame with ADMET properties
        molecule_id (str): ID of molecule to visualize
        id_column (str): Name of ID column. Default: 'ligand_id'
        output_path (str, optional): Path to save plot

    Examples:
        >>> create_admet_radar_plot(df_admet, 'aspirin', output_path='aspirin_admet.png')
    """
    # Find molecule
    mol_data = df[df[id_column] == molecule_id]

    if mol_data.empty:
        raise ValueError(f"Molecule '{molecule_id}' not found in DataFrame")

    mol_data = mol_data.iloc[0]

    # Define properties and their values
    # Normalize all to 0-1 scale where 1 is "good"
    properties = {
        'QED\n(Drug-likeness)': mol_data.get('qed_score', 0),
        'Lipinski\nCompliance': 1.0 if mol_data.get('lipinski_compliant', False) else 0.0,
        'Oral\nBioavailability': {
            'High': 1.0,
            'Medium': 0.5,
            'Low': 0.0
        }.get(mol_data.get('oral_bioavailability', 'Low'), 0.0),
        'BBB\nPenetration': 1.0 if mol_data.get('bbb_penetrant', False) else 0.0,
        'Synthetic\nAccessibility': max(0, (10 - mol_data.get('sa_score', 10)) / 10),  # Invert: lower SA is better
        'Molecular\nWeight': max(0, 1 - mol_data.get('molecular_weight', 500) / 1000)  # Normalize to 0-1
    }

    # Number of properties
    num_vars = len(properties)

    # Compute angle for each axis
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

    # Complete the circle
    values = list(properties.values())
    values += values[:1]
    angles += angles[:1]

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))

    # Draw polygon
    ax.plot(angles, values, 'o-', linewidth=2, label=molecule_id, color='#2E86AB')
    ax.fill(angles, values, alpha=0.25, color='#2E86AB')

    # Fix axis to go in the right order
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    # Draw axis lines for each property
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(properties.keys(), size=10)

    # Set y-axis limits
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], size=8, color='gray')

    # Add grid
    ax.grid(True, linestyle='--', alpha=0.5)

    # Title
    plt.title(f'ADMET Profile: {molecule_id}', size=16, pad=20, fontweight='bold')

    # Legend with property details
    legend_text = f"MW: {mol_data.get('molecular_weight', 0):.1f} Da\n"
    legend_text += f"LogP: {mol_data.get('logp', 0):.2f}\n"
    legend_text += f"QED: {mol_data.get('qed_score', 0):.3f}\n"
    legend_text += f"SA: {mol_data.get('sa_score', 0):.1f}"

    ax.text(1.3, 0.5, legend_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='center',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved ADMET radar plot to: {output_path}")
    else:
        plt.show()

    plt.close()


def plot_property_distribution(
    df: pd.DataFrame,
    property_name: str,
    bins: int = 30,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    output_path: Optional[str] = None
) -> None:
    """
    Plot distribution of a molecular property.

    Args:
        df (pd.DataFrame): DataFrame with molecular properties
        property_name (str): Name of property column to plot
        bins (int): Number of histogram bins. Default: 30
        title (str, optional): Plot title. Auto-generated if None.
        xlabel (str, optional): X-axis label. Uses property_name if None.
        output_path (str, optional): Path to save plot

    Examples:
        >>> plot_property_distribution(
        ...     df_admet,
        ...     'qed_score',
        ...     title='QED Distribution',
        ...     output_path='qed_dist.png'
        ... )
    """
    if property_name not in df.columns:
        raise ValueError(f"Column '{property_name}' not found in DataFrame")

    # Get data (remove NaN)
    data = df[property_name].dropna()

    if len(data) == 0:
        raise ValueError(f"No valid data for property '{property_name}'")

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Histogram
    n, bins_edges, patches = ax.hist(data, bins=bins, alpha=0.7, color='#2E86AB', edgecolor='black')

    # Add mean line
    mean_val = data.mean()
    ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.2f}')

    # Add median line
    median_val = data.median()
    ax.axvline(median_val, color='orange', linestyle='--', linewidth=2, label=f'Median: {median_val:.2f}')

    # Labels
    ax.set_xlabel(xlabel or property_name, fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title(title or f'Distribution of {property_name}', fontsize=14, fontweight='bold')

    # Legend
    ax.legend()

    # Grid
    ax.grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved distribution plot to: {output_path}")
    else:
        plt.show()

    plt.close()


def plot_admet_summary(
    df: pd.DataFrame,
    output_path: Optional[str] = None
) -> None:
    """
    Create comprehensive ADMET summary visualization.

    Creates a multi-panel figure showing:
    - QED distribution
    - LogP distribution
    - Molecular weight distribution
    - Lipinski compliance pie chart
    - BBB penetration pie chart
    - Oral bioavailability bar chart

    Args:
        df (pd.DataFrame): DataFrame with ADMET properties
        output_path (str, optional): Path to save figure

    Examples:
        >>> plot_admet_summary(df_admet, output_path='admet_summary.png')
    """
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)

    # 1. QED Distribution
    ax1 = fig.add_subplot(gs[0, 0])
    if 'qed_score' in df.columns:
        qed_data = df['qed_score'].dropna()
        ax1.hist(qed_data, bins=20, alpha=0.7, color='#2E86AB', edgecolor='black')
        ax1.axvline(qed_data.mean(), color='red', linestyle='--', linewidth=2)
        ax1.set_xlabel('QED Score', fontsize=10)
        ax1.set_ylabel('Count', fontsize=10)
        ax1.set_title('QED Distribution', fontsize=11, fontweight='bold')
        ax1.grid(True, alpha=0.3)

    # 2. LogP Distribution
    ax2 = fig.add_subplot(gs[0, 1])
    if 'logp' in df.columns:
        logp_data = df['logp'].dropna()
        ax2.hist(logp_data, bins=20, alpha=0.7, color='#A23B72', edgecolor='black')
        ax2.axvline(logp_data.mean(), color='red', linestyle='--', linewidth=2)
        ax2.axvline(5, color='orange', linestyle=':', linewidth=2, label='Lipinski limit')
        ax2.set_xlabel('LogP', fontsize=10)
        ax2.set_ylabel('Count', fontsize=10)
        ax2.set_title('LogP Distribution', fontsize=11, fontweight='bold')
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

    # 3. Molecular Weight Distribution
    ax3 = fig.add_subplot(gs[0, 2])
    if 'molecular_weight' in df.columns:
        mw_data = df['molecular_weight'].dropna()
        ax3.hist(mw_data, bins=20, alpha=0.7, color='#F18F01', edgecolor='black')
        ax3.axvline(mw_data.mean(), color='red', linestyle='--', linewidth=2)
        ax3.axvline(500, color='orange', linestyle=':', linewidth=2, label='Lipinski limit')
        ax3.set_xlabel('Molecular Weight (Da)', fontsize=10)
        ax3.set_ylabel('Count', fontsize=10)
        ax3.set_title('MW Distribution', fontsize=11, fontweight='bold')
        ax3.legend(fontsize=8)
        ax3.grid(True, alpha=0.3)

    # 4. Lipinski Compliance Pie Chart
    ax4 = fig.add_subplot(gs[1, 0])
    if 'lipinski_compliant' in df.columns:
        lipinski_counts = df['lipinski_compliant'].value_counts()
        # Create labels based on actual data
        label_map = {True: 'Compliant', False: 'Non-compliant'}
        colors_map = {True: '#2E86AB', False: '#E63946'}
        labels = [label_map[val] for val in lipinski_counts.index]
        colors = [colors_map[val] for val in lipinski_counts.index]
        ax4.pie(lipinski_counts, labels=labels, autopct='%1.1f%%', colors=colors, startangle=90)
        ax4.set_title('Lipinski Compliance', fontsize=11, fontweight='bold')

    # 5. BBB Penetration Pie Chart
    ax5 = fig.add_subplot(gs[1, 1])
    if 'bbb_penetrant' in df.columns:
        bbb_counts = df['bbb_penetrant'].value_counts()
        # Create labels based on actual data
        label_map = {True: 'Penetrant', False: 'Non-penetrant'}
        colors_map = {True: '#06A77D', False: '#D62246'}
        labels = [label_map[val] for val in bbb_counts.index]
        colors = [colors_map[val] for val in bbb_counts.index]
        ax5.pie(bbb_counts, labels=labels, autopct='%1.1f%%', colors=colors, startangle=90)
        ax5.set_title('BBB Penetration', fontsize=11, fontweight='bold')

    # 6. Oral Bioavailability Bar Chart
    ax6 = fig.add_subplot(gs[1, 2])
    if 'oral_bioavailability' in df.columns:
        bioav_counts = df['oral_bioavailability'].value_counts()
        colors_map = {'High': '#2E86AB', 'Medium': '#F18F01', 'Low': '#E63946'}
        colors = [colors_map.get(cat, 'gray') for cat in bioav_counts.index]
        ax6.bar(bioav_counts.index, bioav_counts.values, color=colors, edgecolor='black', alpha=0.7)
        ax6.set_xlabel('Category', fontsize=10)
        ax6.set_ylabel('Count', fontsize=10)
        ax6.set_title('Oral Bioavailability', fontsize=11, fontweight='bold')
        ax6.grid(True, alpha=0.3, axis='y')

    # 7. TPSA Distribution
    ax7 = fig.add_subplot(gs[2, 0])
    if 'tpsa' in df.columns:
        tpsa_data = df['tpsa'].dropna()
        ax7.hist(tpsa_data, bins=20, alpha=0.7, color='#118AB2', edgecolor='black')
        ax7.axvline(tpsa_data.mean(), color='red', linestyle='--', linewidth=2)
        ax7.axvline(140, color='orange', linestyle=':', linewidth=2, label='Veber limit')
        ax7.set_xlabel('TPSA (Ų)', fontsize=10)
        ax7.set_ylabel('Count', fontsize=10)
        ax7.set_title('TPSA Distribution', fontsize=11, fontweight='bold')
        ax7.legend(fontsize=8)
        ax7.grid(True, alpha=0.3)

    # 8. Synthetic Accessibility
    ax8 = fig.add_subplot(gs[2, 1])
    if 'sa_score' in df.columns:
        sa_data = df['sa_score'].dropna()
        ax8.hist(sa_data, bins=20, alpha=0.7, color='#073B4C', edgecolor='black')
        ax8.axvline(sa_data.mean(), color='red', linestyle='--', linewidth=2)
        ax8.set_xlabel('SA Score (1=easy, 10=hard)', fontsize=10)
        ax8.set_ylabel('Count', fontsize=10)
        ax8.set_title('Synthetic Accessibility', fontsize=11, fontweight='bold')
        ax8.grid(True, alpha=0.3)

    # 9. Summary Statistics
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.axis('off')

    summary_text = f"ADMET Summary Statistics\n\n"
    summary_text += f"Total molecules: {len(df)}\n\n"

    if 'lipinski_compliant' in df.columns:
        lipinski_pct = (df['lipinski_compliant'].sum() / len(df)) * 100
        summary_text += f"Lipinski compliant: {lipinski_pct:.1f}%\n"

    if 'bbb_penetrant' in df.columns:
        bbb_pct = (df['bbb_penetrant'].sum() / len(df)) * 100
        summary_text += f"BBB penetrant: {bbb_pct:.1f}%\n"

    if 'oral_bioavailability' in df.columns:
        high_bioav = (df['oral_bioavailability'] == 'High').sum()
        high_bioav_pct = (high_bioav / len(df)) * 100
        summary_text += f"High bioavailability: {high_bioav_pct:.1f}%\n\n"

    if 'qed_score' in df.columns:
        summary_text += f"Mean QED: {df['qed_score'].mean():.3f}\n"

    if 'sa_score' in df.columns:
        summary_text += f"Mean SA: {df['sa_score'].mean():.2f}\n"

    if 'molecular_weight' in df.columns:
        summary_text += f"Mean MW: {df['molecular_weight'].mean():.1f} Da\n"

    if 'logp' in df.columns:
        summary_text += f"Mean LogP: {df['logp'].mean():.2f}\n"

    ax9.text(0.1, 0.9, summary_text, transform=ax9.transAxes,
             fontsize=10, verticalalignment='top',
             fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Main title
    fig.suptitle('ADMET Analysis Summary', fontsize=16, fontweight='bold', y=0.98)

    if output_path:
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved ADMET summary to: {output_path}")
    else:
        plt.show()

    plt.close()
