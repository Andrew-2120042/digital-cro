#!/usr/bin/env python3
"""
Phase 1 Validation Test Script for Digital CRO

This script validates the basic infrastructure by:
1. Loading test molecules from SMILES
2. Calculating molecular properties
3. Filtering based on Lipinski rules
4. Generating visualizations
5. Saving results to data/outputs/phase1_test/

Usage:
    python scripts/tests/test_phase1.py
"""

import sys
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scripts.utils.molecule_ops import smiles_to_mol, draw_molecule, draw_molecule_grid
from scripts.utils.properties import calculate_all_properties, passes_lipinski, get_property_summary


# Test molecules (known drugs)
TEST_MOLECULES = {
    'Aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
    'Caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'Ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
    'Paracetamol': 'CC(=O)Nc1ccc(O)cc1',
    'Penicillin': 'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O',
    'Atorvastatin': 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O',
    'Metformin': 'CN(C)C(=N)NC(=N)N',
    'Warfarin': 'CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O',
    'Sildenafil': 'CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12',
    'Lisinopril': 'NCCCC[C@@H](C(=O)N1CCC[C@H]1C(=O)O)NC(=O)[C@H](CC2=CC=CC=C2)N'
}


def create_output_directory():
    """Create output directory for test results."""
    output_dir = project_root / 'data' / 'outputs' / 'phase1_test'
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def load_test_molecules():
    """
    Load and convert test molecules from SMILES.

    Returns:
        pd.DataFrame: DataFrame with molecule names, SMILES, and RDKit mol objects
    """
    print("=" * 80)
    print("PHASE 1 VALIDATION TEST - Digital CRO")
    print("=" * 80)
    print("\n[1] Loading test molecules...")

    data = []
    for name, smiles in TEST_MOLECULES.items():
        mol = smiles_to_mol(smiles)
        data.append({
            'name': name,
            'smiles': smiles,
            'mol': mol,
            'valid': mol is not None
        })

    df = pd.DataFrame(data)
    valid_count = df['valid'].sum()
    print(f"    ✓ Loaded {len(df)} molecules")
    print(f"    ✓ {valid_count} valid molecules, {len(df) - valid_count} failed")

    return df


def calculate_properties(df):
    """
    Calculate properties for all molecules.

    Args:
        df (pd.DataFrame): DataFrame with molecules

    Returns:
        pd.DataFrame: Updated DataFrame with properties
    """
    print("\n[2] Calculating molecular properties...")

    properties_list = []
    for idx, row in df.iterrows():
        if row['valid']:
            try:
                props = calculate_all_properties(row['mol'])
                props['name'] = row['name']
                props['smiles'] = row['smiles']
                properties_list.append(props)
            except Exception as e:
                print(f"    ✗ Error calculating properties for {row['name']}: {e}")
                properties_list.append({'name': row['name'], 'error': str(e)})
        else:
            properties_list.append({'name': row['name'], 'error': 'Invalid SMILES'})

    props_df = pd.DataFrame(properties_list)
    print(f"    ✓ Calculated properties for {len(properties_list)} molecules")

    return props_df


def filter_lipinski(props_df):
    """
    Filter molecules based on Lipinski rules.

    Args:
        props_df (pd.DataFrame): DataFrame with molecular properties

    Returns:
        tuple: (passing molecules DataFrame, failing molecules DataFrame)
    """
    print("\n[3] Filtering by Lipinski's Rule of Five...")

    if 'passes_lipinski' in props_df.columns:
        passing = props_df[props_df['passes_lipinski'] == True]
        failing = props_df[props_df['passes_lipinski'] == False]

        print(f"    ✓ {len(passing)} molecules pass Lipinski rules")
        print(f"    ✓ {len(failing)} molecules fail Lipinski rules")

        if len(passing) > 0:
            print("\n    Passing molecules:")
            for name in passing['name'].values:
                print(f"      • {name}")

        if len(failing) > 0:
            print("\n    Failing molecules:")
            for name in failing['name'].values:
                print(f"      • {name}")

        return passing, failing
    else:
        print("    ✗ No Lipinski data available")
        return pd.DataFrame(), pd.DataFrame()


def generate_individual_images(df, output_dir):
    """
    Generate individual PNG images for each molecule.

    Args:
        df (pd.DataFrame): DataFrame with molecules
        output_dir (Path): Output directory path
    """
    print("\n[4] Generating individual molecule images...")

    individual_dir = output_dir / 'individual_molecules'
    individual_dir.mkdir(exist_ok=True)

    count = 0
    for idx, row in df.iterrows():
        if row['valid']:
            try:
                filename = individual_dir / f"{row['name'].replace(' ', '_')}.png"
                draw_molecule(row['mol'], size=(400, 400), filename=str(filename))
                count += 1
            except Exception as e:
                print(f"    ✗ Error drawing {row['name']}: {e}")

    print(f"    ✓ Generated {count} individual molecule images")
    print(f"    ✓ Saved to: {individual_dir}")


def generate_molecule_grid(df, output_dir):
    """
    Generate grid image of all molecules.

    Args:
        df (pd.DataFrame): DataFrame with molecules
        output_dir (Path): Output directory path
    """
    print("\n[5] Generating molecule grid...")

    valid_mols = df[df['valid']]
    if len(valid_mols) > 0:
        try:
            mols = valid_mols['mol'].tolist()
            labels = valid_mols['name'].tolist()
            filename = output_dir / 'molecule_grid.png'

            draw_molecule_grid(
                mols,
                labels=labels,
                mols_per_row=5,
                filename=str(filename),
                sub_img_size=(250, 250)
            )

            print(f"    ✓ Generated grid with {len(mols)} molecules")
            print(f"    ✓ Saved to: {filename}")
        except Exception as e:
            print(f"    ✗ Error generating grid: {e}")
    else:
        print("    ✗ No valid molecules to display")


def save_properties_csv(props_df, output_dir):
    """
    Save properties to CSV file.

    Args:
        props_df (pd.DataFrame): DataFrame with properties
        output_dir (Path): Output directory path
    """
    print("\n[6] Saving properties to CSV...")

    filename = output_dir / 'properties.csv'
    props_df.to_csv(filename, index=False)
    print(f"    ✓ Saved properties for {len(props_df)} molecules")
    print(f"    ✓ Saved to: {filename}")


def generate_property_distributions(props_df, output_dir):
    """
    Generate property distribution plots.

    Args:
        props_df (pd.DataFrame): DataFrame with properties
        output_dir (Path): Output directory path
    """
    print("\n[7] Generating property distribution plots...")

    # Filter valid data
    valid_props = props_df[~props_df['molecular_weight'].isna()]

    if len(valid_props) == 0:
        print("    ✗ No valid property data for plotting")
        return

    # Create figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Molecular Property Distributions - Phase 1 Test', fontsize=16, fontweight='bold')

    # 1. Molecular Weight
    axes[0, 0].hist(valid_props['molecular_weight'], bins=10, color='steelblue', edgecolor='black', alpha=0.7)
    axes[0, 0].axvline(500, color='red', linestyle='--', linewidth=2, label='Lipinski limit')
    axes[0, 0].set_xlabel('Molecular Weight (Da)', fontweight='bold')
    axes[0, 0].set_ylabel('Frequency', fontweight='bold')
    axes[0, 0].set_title('Molecular Weight Distribution')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # 2. LogP
    axes[0, 1].hist(valid_props['logp'], bins=10, color='forestgreen', edgecolor='black', alpha=0.7)
    axes[0, 1].axvline(5.6, color='red', linestyle='--', linewidth=2, label='Lipinski limit')
    axes[0, 1].set_xlabel('LogP', fontweight='bold')
    axes[0, 1].set_ylabel('Frequency', fontweight='bold')
    axes[0, 1].set_title('LogP Distribution')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # 3. H-bond Donors
    axes[0, 2].hist(valid_props['hbd'], bins=range(int(valid_props['hbd'].max()) + 2),
                    color='coral', edgecolor='black', alpha=0.7, align='left')
    axes[0, 2].axvline(5.5, color='red', linestyle='--', linewidth=2, label='Lipinski limit')
    axes[0, 2].set_xlabel('H-bond Donors', fontweight='bold')
    axes[0, 2].set_ylabel('Frequency', fontweight='bold')
    axes[0, 2].set_title('H-bond Donor Distribution')
    axes[0, 2].legend()
    axes[0, 2].grid(True, alpha=0.3)

    # 4. H-bond Acceptors
    axes[1, 0].hist(valid_props['hba'], bins=range(int(valid_props['hba'].max()) + 2),
                    color='mediumpurple', edgecolor='black', alpha=0.7, align='left')
    axes[1, 0].axvline(10.5, color='red', linestyle='--', linewidth=2, label='Lipinski limit')
    axes[1, 0].set_xlabel('H-bond Acceptors', fontweight='bold')
    axes[1, 0].set_ylabel('Frequency', fontweight='bold')
    axes[1, 0].set_title('H-bond Acceptor Distribution')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    # 5. TPSA
    axes[1, 1].hist(valid_props['tpsa'], bins=10, color='gold', edgecolor='black', alpha=0.7)
    axes[1, 1].set_xlabel('TPSA (Ų)', fontweight='bold')
    axes[1, 1].set_ylabel('Frequency', fontweight='bold')
    axes[1, 1].set_title('TPSA Distribution')
    axes[1, 1].grid(True, alpha=0.3)

    # 6. Lipinski Pass/Fail
    lipinski_counts = valid_props['passes_lipinski'].value_counts()
    axes[1, 2].bar(['Pass', 'Fail'],
                   [lipinski_counts.get(True, 0), lipinski_counts.get(False, 0)],
                   color=['green', 'red'], alpha=0.7, edgecolor='black')
    axes[1, 2].set_ylabel('Count', fontweight='bold')
    axes[1, 2].set_title('Lipinski Rule of Five')
    axes[1, 2].grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    filename = output_dir / 'property_distributions.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    ✓ Generated property distribution plots")
    print(f"    ✓ Saved to: {filename}")


def print_summary_statistics(props_df):
    """
    Print summary statistics for molecular properties.

    Args:
        props_df (pd.DataFrame): DataFrame with properties
    """
    print("\n[8] Summary Statistics:")
    print("=" * 80)

    valid_props = props_df[~props_df['molecular_weight'].isna()]

    if len(valid_props) == 0:
        print("    ✗ No valid property data")
        return

    # Key properties
    properties_to_summarize = [
        ('molecular_weight', 'Molecular Weight (Da)'),
        ('logp', 'LogP'),
        ('hbd', 'H-bond Donors'),
        ('hba', 'H-bond Acceptors'),
        ('tpsa', 'TPSA (Ų)'),
        ('rotatable_bonds', 'Rotatable Bonds')
    ]

    for prop_key, prop_name in properties_to_summarize:
        if prop_key in valid_props.columns:
            values = valid_props[prop_key]
            print(f"\n{prop_name}:")
            print(f"  Mean:   {values.mean():.2f}")
            print(f"  Median: {values.median():.2f}")
            print(f"  Min:    {values.min():.2f}")
            print(f"  Max:    {values.max():.2f}")
            print(f"  Std:    {values.std():.2f}")

    # Lipinski summary
    if 'passes_lipinski' in valid_props.columns:
        pass_count = valid_props['passes_lipinski'].sum()
        pass_rate = (pass_count / len(valid_props)) * 100
        print(f"\nLipinski Rule of Five:")
        print(f"  Passing: {pass_count}/{len(valid_props)} ({pass_rate:.1f}%)")

    print("=" * 80)


def main():
    """Main execution function."""
    try:
        # Create output directory
        output_dir = create_output_directory()

        # Load test molecules
        df = load_test_molecules()

        # Calculate properties
        props_df = calculate_properties(df)

        # Filter by Lipinski rules
        passing, failing = filter_lipinski(props_df)

        # Generate visualizations
        generate_individual_images(df, output_dir)
        generate_molecule_grid(df, output_dir)
        generate_property_distributions(props_df, output_dir)

        # Save properties
        save_properties_csv(props_df, output_dir)

        # Print summary
        print_summary_statistics(props_df)

        print("\n" + "=" * 80)
        print("PHASE 1 VALIDATION COMPLETE")
        print("=" * 80)
        print(f"\nAll results saved to: {output_dir}")
        print("\nGenerated files:")
        print(f"  • individual_molecules/*.png - Individual molecule images")
        print(f"  • molecule_grid.png - Grid of all molecules")
        print(f"  • properties.csv - All calculated properties")
        print(f"  • property_distributions.png - Property distribution plots")
        print("\n✓ Phase 1 infrastructure validated successfully!")
        print("=" * 80)

        return 0

    except Exception as e:
        print(f"\n✗ ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
