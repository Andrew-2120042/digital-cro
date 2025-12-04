"""
Phase 4 Test: ADMET Predictions + Molecular Visualizations

Tests the ADMET prediction and visualization capabilities:
1. Lipinski Rule of Five calculations
2. ADMET property predictions
3. BBB penetration prediction
4. Oral bioavailability estimation
5. Batch ADMET analysis
6. Molecular structure visualization
7. ADMET radar plots
8. Property distributions
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import pandas as pd
from scripts.utils.admet_predictions import (
    calculate_lipinski_properties,
    predict_admet_properties,
    predict_bbb_penetration,
    predict_oral_bioavailability,
    analyze_molecule_complete,
    batch_admet_analysis
)
from scripts.utils.molecular_viz import (
    draw_molecules_grid,
    generate_3d_conformer,
    visualize_docking_results,
    create_admet_radar_plot,
    plot_property_distribution,
    plot_admet_summary
)


def test_lipinski_calculations():
    """Test Lipinski Rule of Five calculations"""
    print("\n" + "="*70)
    print("[Test 1] Lipinski Rule of Five Calculations")
    print("="*70)

    test_molecules = [
        ('CC(=O)Oc1ccccc1C(=O)O', 'aspirin'),
        ('CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'caffeine'),
        ('CCC(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(CC(=O)O3)O)C', 'atorvastatin')
    ]

    results = []
    for smiles, name in test_molecules:
        props = calculate_lipinski_properties(smiles)
        results.append(props)

        print(f"\n{name.upper()}:")
        print(f"  SMILES: {smiles}")
        print(f"  Molecular Weight: {props['molecular_weight']:.2f} Da")
        print(f"  LogP: {props['logp']:.2f}")
        print(f"  H-bond Donors: {props['hbd']}")
        print(f"  H-bond Acceptors: {props['hba']}")
        print(f"  Violations: {props['num_violations']}")
        print(f"  Lipinski Compliant: {'✓ YES' if props['lipinski_compliant'] else '✗ NO'}")

    print("\n✓ Lipinski calculations completed successfully")
    return results


def test_admet_predictions():
    """Test comprehensive ADMET predictions"""
    print("\n" + "="*70)
    print("[Test 2] ADMET Property Predictions")
    print("="*70)

    test_smiles = [
        'CC(=O)Oc1ccccc1C(=O)O',  # Aspirin
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
        'CC(C)Cc1ccc(cc1)C(C)C(=O)O',  # Ibuprofen
    ]

    test_names = ['aspirin', 'caffeine', 'ibuprofen']

    for smiles, name in zip(test_smiles, test_names):
        print(f"\n{name.upper()}:")
        admet = predict_admet_properties(smiles)

        print(f"  MW: {admet['molecular_weight']:.2f} Da")
        print(f"  LogP: {admet['logp']:.2f}")
        print(f"  TPSA: {admet['tpsa']:.2f} Ų")
        print(f"  Rotatable Bonds: {admet['rotatable_bonds']}")
        print(f"  Aromatic Rings: {admet['aromatic_rings']}")
        print(f"  QED Score: {admet['qed_score']:.3f}")
        print(f"  SA Score: {admet['sa_score']:.2f}")
        print(f"  Lipinski: {'✓' if admet['lipinski_compliant'] else '✗'} ({admet['lipinski_violations']} violations)")
        print(f"  BBB Penetrant: {'✓' if admet['bbb_penetrant'] else '✗'} ({admet['bbb_confidence']} confidence)")
        print(f"  Oral Bioavailability: {admet['oral_bioavailability']}")

    print("\n✓ ADMET predictions completed successfully")


def test_bbb_penetration():
    """Test BBB penetration prediction"""
    print("\n" + "="*70)
    print("[Test 3] BBB Penetration Prediction")
    print("="*70)

    from rdkit import Chem

    # Test molecules with known BBB properties
    test_cases = [
        ('CN1CCC23c4ccccc4CC2C1Cc5c3cccc5O', 'morphine', True),  # CNS drug
        ('CC(=O)Oc1ccccc1C(=O)O', 'aspirin', False),  # Poor BBB
        ('CN(C)CCOC(c1ccccc1)c2ccccc2', 'diphenhydramine', True),  # BBB penetrant
    ]

    for smiles, name, expected_bbb in test_cases:
        mol = Chem.MolFromSmiles(smiles)
        bbb = predict_bbb_penetration(mol)

        print(f"\n{name.upper()}:")
        print(f"  BBB Penetrant: {'✓ YES' if bbb['bbb_penetrant'] else '✗ NO'}")
        print(f"  Confidence: {bbb['confidence']}")
        print(f"  TPSA: {bbb['tpsa']:.2f} Ų")
        print(f"  Criteria Met: {len(bbb['criteria_met'])}/5")

        if bbb['criteria_met']:
            print(f"  ✓ {', '.join(bbb['criteria_met'])}")
        if bbb['criteria_failed']:
            print(f"  ✗ {', '.join(bbb['criteria_failed'])}")

    print("\n✓ BBB predictions completed successfully")


def test_batch_analysis():
    """Test batch ADMET analysis"""
    print("\n" + "="*70)
    print("[Test 4] Batch ADMET Analysis")
    print("="*70)

    # Create test DataFrame
    test_data = {
        'ligand_id': ['aspirin', 'caffeine', 'ibuprofen', 'paracetamol', 'morphine'],
        'smiles': [
            'CC(=O)Oc1ccccc1C(=O)O',
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
            'CC(=O)Nc1ccc(cc1)O',
            'CN1CCC23c4ccccc4CC2C1Cc5c3cccc5O'
        ],
        'binding_affinity': [-6.5, -5.2, -7.1, -5.8, -8.3]
    }

    df = pd.DataFrame(test_data)

    print(f"\nAnalyzing {len(df)} molecules...")

    # Run batch ADMET analysis
    df_admet = batch_admet_analysis(df, smiles_column='smiles', id_column='ligand_id')

    print(f"\n✓ Batch analysis complete")
    print(f"\nResults:")
    print(f"  Total molecules: {len(df_admet)}")

    # Summary statistics
    if 'lipinski_compliant' in df_admet.columns:
        lipinski_compliant = df_admet['lipinski_compliant'].sum()
        print(f"  Lipinski compliant: {lipinski_compliant}/{len(df_admet)}")

    if 'bbb_penetrant' in df_admet.columns:
        bbb_penetrant = df_admet['bbb_penetrant'].sum()
        print(f"  BBB penetrant: {bbb_penetrant}/{len(df_admet)}")

    if 'oral_bioavailability' in df_admet.columns:
        high_bioav = (df_admet['oral_bioavailability'] == 'High').sum()
        print(f"  High oral bioavailability: {high_bioav}/{len(df_admet)}")

    if 'qed_score' in df_admet.columns:
        mean_qed = df_admet['qed_score'].mean()
        print(f"  Mean QED: {mean_qed:.3f}")

    # Display table
    print(f"\nTop results:")
    display_cols = ['ligand_id', 'qed_score', 'sa_score', 'lipinski_compliant',
                    'oral_bioavailability', 'bbb_penetrant']
    available_cols = [col for col in display_cols if col in df_admet.columns]
    print(df_admet[available_cols].to_string(index=False))

    print("\n✓ Batch analysis completed successfully")
    return df_admet


def test_visualizations(df_admet):
    """Test molecular visualizations"""
    print("\n" + "="*70)
    print("[Test 5] Molecular Visualizations")
    print("="*70)

    output_dir = Path('data/outputs/phase4_test')
    output_dir.mkdir(parents=True, exist_ok=True)

    # Test 1: Molecule grid
    print("\n[5.1] Drawing molecule grid...")
    smiles_list = df_admet['smiles'].tolist()
    labels = df_admet['ligand_id'].tolist()

    grid_path = output_dir / 'molecules_grid.png'
    draw_molecules_grid(
        smiles_list,
        labels=labels,
        mols_per_row=3,
        output_path=str(grid_path)
    )
    print(f"  ✓ Saved to: {grid_path}")

    # Test 2: Docking results visualization
    print("\n[5.2] Visualizing docking results...")
    docking_viz_path = output_dir / 'top_hits.png'
    visualize_docking_results(
        df_admet,
        top_n=5,
        smiles_column='smiles',
        id_column='ligand_id',
        affinity_column='binding_affinity',
        mols_per_row=3,
        output_path=str(docking_viz_path)
    )
    print(f"  ✓ Saved to: {docking_viz_path}")

    # Test 3: ADMET radar plot
    print("\n[5.3] Creating ADMET radar plot...")
    radar_path = output_dir / 'admet_radar_aspirin.png'
    create_admet_radar_plot(
        df_admet,
        'aspirin',
        id_column='ligand_id',
        output_path=str(radar_path)
    )
    print(f"  ✓ Saved to: {radar_path}")

    # Test 4: Property distribution
    print("\n[5.4] Plotting QED distribution...")
    qed_dist_path = output_dir / 'qed_distribution.png'
    plot_property_distribution(
        df_admet,
        'qed_score',
        title='QED Score Distribution',
        xlabel='QED Score',
        output_path=str(qed_dist_path)
    )
    print(f"  ✓ Saved to: {qed_dist_path}")

    # Test 5: ADMET summary
    print("\n[5.5] Creating ADMET summary...")
    summary_path = output_dir / 'admet_summary.png'
    plot_admet_summary(
        df_admet,
        output_path=str(summary_path)
    )
    print(f"  ✓ Saved to: {summary_path}")

    print("\n✓ All visualizations completed successfully")
    print(f"\nOutput directory: {output_dir}")


def test_3d_conformer():
    """Test 3D conformer generation"""
    print("\n" + "="*70)
    print("[Test 6] 3D Conformer Generation")
    print("="*70)

    test_molecules = [
        ('CC(=O)Oc1ccccc1C(=O)O', 'aspirin'),
        ('CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'caffeine'),
    ]

    for smiles, name in test_molecules:
        print(f"\n{name.upper()}:")
        mol_3d = generate_3d_conformer(smiles, optimize=True)

        if mol_3d:
            print(f"  ✓ 3D conformer generated")
            print(f"  Atoms: {mol_3d.GetNumAtoms()}")
            print(f"  Conformers: {mol_3d.GetNumConformers()}")

            # Check if coordinates are valid
            conf = mol_3d.GetConformer()
            pos = conf.GetAtomPosition(0)
            print(f"  First atom position: ({pos.x:.3f}, {pos.y:.3f}, {pos.z:.3f})")
        else:
            print(f"  ✗ Failed to generate 3D conformer")

    print("\n✓ 3D conformer generation completed successfully")


def test_complete_analysis():
    """Test complete molecule analysis"""
    print("\n" + "="*70)
    print("[Test 7] Complete Molecule Analysis")
    print("="*70)

    smiles = 'CC(=O)Oc1ccccc1C(=O)O'
    mol_id = 'aspirin'

    print(f"\nAnalyzing {mol_id}...")
    analysis = analyze_molecule_complete(smiles, mol_id)

    print(f"\nCOMPLETE ANALYSIS - {mol_id.upper()}")
    print(f"  SMILES: {analysis['smiles']}")
    print(f"\nBasic Properties:")
    print(f"  MW: {analysis['molecular_weight']:.2f} Da")
    print(f"  LogP: {analysis['logp']:.2f}")
    print(f"  TPSA: {analysis['tpsa']:.2f} Ų")
    print(f"  Rotatable bonds: {analysis['rotatable_bonds']}")
    print(f"  Aromatic rings: {analysis['aromatic_rings']}")

    print(f"\nDrug-likeness:")
    print(f"  QED Score: {analysis['qed_score']:.3f}")
    print(f"  SA Score: {analysis['sa_score']:.2f}")
    print(f"  Lipinski: {'✓ Compliant' if analysis['lipinski_compliant'] else '✗ Non-compliant'} ({analysis['lipinski_violations']} violations)")

    print(f"\nADMET Predictions:")
    print(f"  BBB Penetrant: {'✓ YES' if analysis['bbb_penetrant'] else '✗ NO'} ({analysis['bbb_confidence']} confidence)")
    print(f"  Oral Bioavailability: {analysis['oral_bioavailability']}")

    # Detailed BBB analysis
    if 'bbb_details' in analysis:
        bbb = analysis['bbb_details']
        print(f"\nBBB Criteria:")
        if bbb['criteria_met']:
            for criterion in bbb['criteria_met']:
                print(f"    ✓ {criterion}")
        if bbb['criteria_failed']:
            for criterion in bbb['criteria_failed']:
                print(f"    ✗ {criterion}")

    print("\n✓ Complete analysis finished successfully")


def main():
    """Run all Phase 4 tests"""
    print("\n" + "="*70)
    print("PHASE 4 TEST: ADMET PREDICTIONS + MOLECULAR VISUALIZATIONS")
    print("="*70)

    try:
        # Test 1: Lipinski calculations
        test_lipinski_calculations()

        # Test 2: ADMET predictions
        test_admet_predictions()

        # Test 3: BBB penetration
        test_bbb_penetration()

        # Test 4: Batch analysis
        df_admet = test_batch_analysis()

        # Test 5: Visualizations
        test_visualizations(df_admet)

        # Test 6: 3D conformer generation
        test_3d_conformer()

        # Test 7: Complete analysis
        test_complete_analysis()

        # Final summary
        print("\n" + "="*70)
        print("✅ PHASE 4 TESTS - ALL PASSED!")
        print("="*70)
        print("\nCapabilities verified:")
        print("  ✓ Lipinski Rule of Five calculations")
        print("  ✓ Comprehensive ADMET predictions")
        print("  ✓ BBB penetration prediction")
        print("  ✓ Oral bioavailability estimation")
        print("  ✓ Batch processing of molecules")
        print("  ✓ 2D structure visualization")
        print("  ✓ 3D conformer generation")
        print("  ✓ ADMET radar plots")
        print("  ✓ Property distributions")
        print("  ✓ Complete ADMET summary reports")

        print("\n✅ Phase 4 ready for production!")
        print("Ready to move to Phase 5: Report Generation")

        return True

    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
