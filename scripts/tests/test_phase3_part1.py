"""
PHASE 3.2 VALIDATION TEST: AutoDock Vina Integration

This test validates:
1. AutoDock Vina installation
2. Receptor preparation (PDB → PDBQT)
3. Ligand preparation (SMILES → PDBQT)
4. Single molecule docking
5. Result parsing

Test protein: 1HSG (HIV-1 Protease)
Test ligand: Aspirin (simple molecule)

Expected results:
- Receptor PDBQT file created with proper format
- Ligand PDBQT file created with proper format
- Vina docking completes successfully
- Binding affinity extracted from results
"""

import os
import sys
import json
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from scripts.utils import (
    # Phase 2: Protein and pocket detection
    download_pdb,
    load_pdb,
    clean_pdb,
    detect_pockets_auto,
    # Phase 3: Ligand preparation
    prepare_ligand_for_docking,
    # Phase 3: Receptor preparation
    prepare_receptor_for_docking,
    # Phase 3: Vina wrapper
    check_vina_installed,
    dock_ligand_with_pocket
)


def main():
    """Run Phase 3.2 validation tests."""

    print("="*80)
    print("PHASE 3.2 VALIDATION TEST: AutoDock Vina Integration")
    print("="*80)
    print()

    # Setup paths
    base_dir = Path(__file__).parent.parent.parent
    output_dir = base_dir / 'data' / 'outputs' / 'phase3_test'
    output_dir.mkdir(parents=True, exist_ok=True)

    proteins_dir = output_dir / 'proteins'
    ligands_dir = output_dir / 'ligands'
    docking_dir = output_dir / 'docking'

    proteins_dir.mkdir(exist_ok=True)
    ligands_dir.mkdir(exist_ok=True)
    docking_dir.mkdir(exist_ok=True)

    # Test results
    test_results = {
        'vina_installed': False,
        'receptor_prepared': False,
        'ligand_prepared': False,
        'docking_completed': False,
        'errors': []
    }

    #
    # TEST 1: Check Vina Installation
    #
    print("-" * 80)
    print("TEST 1: Checking AutoDock Vina Installation")
    print("-" * 80)

    vina_status = check_vina_installed()

    if vina_status['installed']:
        print(f"✓ AutoDock Vina installed")
        print(f"  Version: {vina_status['version']}")
        print(f"  Path: {vina_status['path']}")
        test_results['vina_installed'] = True
    else:
        print(f"✗ AutoDock Vina NOT installed")
        print(f"  Error: {vina_status['error']}")
        print()
        print("To install Vina:")
        print("  conda install -c conda-forge autodock-vina")
        print()
        test_results['errors'].append(f"Vina not installed: {vina_status['error']}")

        # Cannot proceed without Vina
        print("\nCANNOT PROCEED: Vina is required for docking tests")
        print_summary(test_results)
        return

    print()

    #
    # TEST 2: Prepare Receptor (1HSG → PDBQT)
    #
    print("-" * 80)
    print("TEST 2: Receptor Preparation (1HSG HIV-1 Protease)")
    print("-" * 80)

    pdb_id = '1HSG'

    # Download PDB if needed
    pdb_file = proteins_dir / f'{pdb_id}.pdb'
    if not pdb_file.exists():
        print(f"Downloading {pdb_id}...")
        download_pdb(pdb_id, str(proteins_dir))
    else:
        print(f"Using existing PDB: {pdb_file}")

    # Prepare receptor
    receptor_pdbqt = proteins_dir / f'{pdb_id}_receptor.pdbqt'

    print("Preparing receptor (PDB → PDBQT)...")
    receptor_result = prepare_receptor_for_docking(
        pdb_path=str(pdb_file),
        output_path=str(receptor_pdbqt),
        remove_water=True,
        remove_hetero=True,
        add_hydrogens=True
    )

    if receptor_result['success']:
        print(f"✓ Receptor prepared successfully")
        print(f"  Output: {receptor_pdbqt}")
        print(f"  Atoms: {receptor_result['num_atoms']}")
        print(f"  Residues: {receptor_result['num_residues']}")
        print(f"  Chains: {', '.join(receptor_result['chains'])}")
        test_results['receptor_prepared'] = True
    else:
        print(f"✗ Receptor preparation FAILED")
        print(f"  Error: {receptor_result['error']}")
        test_results['errors'].append(f"Receptor prep failed: {receptor_result['error']}")
        print_summary(test_results)
        return

    print()

    #
    # TEST 3: Detect Binding Pocket
    #
    print("-" * 80)
    print("TEST 3: Binding Pocket Detection")
    print("-" * 80)

    # Load PDB structure
    structure = load_pdb(str(pdb_file))

    # Detect pockets (will use ligand-based method for 1HSG)
    print("Detecting binding pocket...")
    pockets = detect_pockets_auto(structure, prefer_ligand=True)

    if pockets:
        pocket = pockets[0]
        print(f"✓ Pocket detected")
        print(f"  Method: {pocket.get('method', 'unknown')}")
        print(f"  Center: ({pocket['center'][0]:.2f}, {pocket['center'][1]:.2f}, {pocket['center'][2]:.2f})")
        print(f"  Volume: {pocket['volume']:.1f} Ų")
        print(f"  Score: {pocket.get('score', 'N/A')}")
    else:
        print(f"✗ No pockets detected")
        test_results['errors'].append("No pockets detected")
        print_summary(test_results)
        return

    print()

    #
    # TEST 4: Prepare Ligand (Aspirin SMILES → PDBQT)
    #
    print("-" * 80)
    print("TEST 4: Ligand Preparation (Aspirin)")
    print("-" * 80)

    # Aspirin SMILES
    aspirin_smiles = 'CC(=O)Oc1ccccc1C(=O)O'
    ligand_pdbqt = ligands_dir / 'aspirin.pdbqt'

    print(f"Preparing ligand: {aspirin_smiles}")
    ligand_result = prepare_ligand_for_docking(
        smiles=aspirin_smiles,
        output_path=str(ligand_pdbqt),
        mol_id='ASA'
    )

    if ligand_result['success']:
        print(f"✓ Ligand prepared successfully")
        print(f"  Output: {ligand_pdbqt}")
        print(f"  Atoms: {ligand_result['num_atoms']}")
        print(f"  Rotatable bonds: {ligand_result['num_rotatable_bonds']}")
        print(f"  Molecular weight: {ligand_result['molecular_weight']:.1f} Da")
        test_results['ligand_prepared'] = True
    else:
        print(f"✗ Ligand preparation FAILED")
        print(f"  Error: {ligand_result['error']}")
        test_results['errors'].append(f"Ligand prep failed: {ligand_result['error']}")
        print_summary(test_results)
        return

    print()

    #
    # TEST 5: Run Docking
    #
    print("-" * 80)
    print("TEST 5: AutoDock Vina Docking")
    print("-" * 80)

    output_pdbqt = docking_dir / 'aspirin_docked.pdbqt'

    print(f"Running Vina docking...")
    print(f"  Receptor: {receptor_pdbqt.name}")
    print(f"  Ligand: {ligand_pdbqt.name}")
    print(f"  Search box center: ({pocket['center'][0]:.1f}, {pocket['center'][1]:.1f}, {pocket['center'][2]:.1f})")
    print(f"  This may take 1-2 minutes...")
    print()

    docking_result = dock_ligand_with_pocket(
        receptor_pdbqt=str(receptor_pdbqt),
        ligand_pdbqt=str(ligand_pdbqt),
        output_pdbqt=str(output_pdbqt),
        pocket=pocket,
        padding=5.0,
        exhaustiveness=8,
        num_modes=9
    )

    if docking_result['success']:
        print(f"✓ Docking completed successfully")
        print(f"  Output: {output_pdbqt}")
        print(f"  Best affinity: {docking_result['best_affinity']:.2f} kcal/mol")
        print()
        print(f"  All binding modes ({len(docking_result['modes'])}):")
        for mode in docking_result['modes']:
            print(f"    Mode {mode['mode']}: {mode['affinity']:7.2f} kcal/mol  "
                  f"(RMSD l.b.={mode['rmsd_lb']:.3f}, u.b.={mode['rmsd_ub']:.3f})")
        test_results['docking_completed'] = True
    else:
        print(f"✗ Docking FAILED")
        print(f"  Error: {docking_result['error']}")
        if docking_result['log']:
            print(f"\n  Vina log:")
            print("  " + "\n  ".join(docking_result['log'].split('\n')[:20]))
        test_results['errors'].append(f"Docking failed: {docking_result['error']}")

    print()

    #
    # Save Results
    #
    print("-" * 80)
    print("Saving Results")
    print("-" * 80)

    results_file = output_dir / 'phase3_test_results.json'

    results_data = {
        'vina_status': vina_status,
        'receptor': receptor_result,
        'pocket': {
            'center': [float(x) for x in pocket['center']],
            'volume': float(pocket['volume']),
            'method': pocket.get('method', 'unknown')
        },
        'ligand': ligand_result,
        'docking': {
            'success': docking_result['success'],
            'best_affinity': docking_result.get('best_affinity'),
            'num_modes': len(docking_result.get('modes', [])),
            'modes': docking_result.get('modes', []),
            'error': docking_result.get('error')
        }
    }

    with open(results_file, 'w') as f:
        json.dump(results_data, f, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print()

    #
    # Summary
    #
    print_summary(test_results)


def print_summary(test_results):
    """Print test summary."""
    print("="*80)
    print("PHASE 3.2 TEST SUMMARY")
    print("="*80)
    print()

    all_passed = all([
        test_results['vina_installed'],
        test_results['receptor_prepared'],
        test_results['ligand_prepared'],
        test_results['docking_completed']
    ])

    print(f"1. Vina Installation:    {'✓ PASS' if test_results['vina_installed'] else '✗ FAIL'}")
    print(f"2. Receptor Preparation: {'✓ PASS' if test_results['receptor_prepared'] else '✗ FAIL'}")
    print(f"3. Ligand Preparation:   {'✓ PASS' if test_results['ligand_prepared'] else '✗ FAIL'}")
    print(f"4. Vina Docking:         {'✓ PASS' if test_results['docking_completed'] else '✗ FAIL'}")
    print()

    if all_passed:
        print("="*80)
        print("✓ ALL TESTS PASSED - Phase 3.2 Complete!")
        print("="*80)
        print()
        print("Ready for Phase 3.3: Batch Docking + Result Analysis")
    else:
        print("="*80)
        print("✗ SOME TESTS FAILED")
        print("="*80)
        if test_results['errors']:
            print("\nErrors:")
            for error in test_results['errors']:
                print(f"  - {error}")

    print()


if __name__ == '__main__':
    main()
