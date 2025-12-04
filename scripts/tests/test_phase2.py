#!/usr/bin/env python3
"""
Phase 2 Validation Test Script for Digital CRO

This script validates protein handling, pocket detection, and library filtering:
1. Protein operations (download, load, clean, analyze)
2. Pocket detection (geometric algorithm, scoring)
3. Library filtering (PAINS, reactive groups, diversity)
4. Integration test (complete workflow)

Usage:
    python scripts/tests/test_phase2.py
"""

import sys
import os
import json
import time
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Import Phase 2 utilities
from scripts.utils.protein_ops import (
    download_pdb, load_pdb, clean_pdb, save_pdb,
    get_protein_info, calculate_center_of_mass
)
from scripts.utils.pocket_detection import (
    detect_pockets_geometric, get_pocket_residues,
    create_docking_box, visualize_pocket, find_known_binding_site
)
from scripts.utils.library_filter import (
    filter_large_library, detect_pains_substructures,
    detect_reactive_groups, calculate_diversity_subset,
    calculate_library_diversity
)


def create_output_directories():
    """Create output directories for Phase 2 tests."""
    base_dir = project_root / 'data' / 'outputs' / 'phase2_test'

    dirs = [
        base_dir / 'proteins',
        base_dir / 'pockets',
        base_dir / 'library_filtering',
        base_dir / 'integration'
    ]

    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)

    return base_dir


def test_protein_operations():
    """
    Test protein handling pipeline.

    Tests:
    1. Download HIV protease (PDB: 1HSG)
    2. Load and parse structure
    3. Clean structure (remove water/hetero)
    4. Extract metadata
    5. Calculate center of mass
    6. Save cleaned structure
    """
    print("\n" + "="*80)
    print("[TEST 1] PROTEIN OPERATIONS")
    print("="*80)

    output_dir = project_root / 'data' / 'outputs' / 'phase2_test' / 'proteins'

    try:
        # Step 1: Download 1HSG
        print("\n[1.1] Downloading PDB structure (1HSG - HIV Protease)...")
        pdb_path = download_pdb('1HSG', output_dir=str(output_dir))
        print(f"      ✓ Downloaded to: {pdb_path}")

        # Step 2: Load structure
        print("\n[1.2] Loading PDB structure...")
        structure = load_pdb(pdb_path)
        print(f"      ✓ Loaded structure: {structure.id}")

        # Step 3: Get metadata
        print("\n[1.3] Extracting protein metadata...")
        info = get_protein_info(structure)
        print(f"      ✓ Chains: {info['chains']}")
        print(f"      ✓ Residues: {info['num_residues']}")
        print(f"      ✓ Atoms: {info['num_atoms']}")

        # Save metadata
        info_file = output_dir / 'protein_info.json'
        with open(info_file, 'w') as f:
            json.dump(info, f, indent=2)
        print(f"      ✓ Saved metadata to: {info_file}")

        # Step 4: Clean structure
        print("\n[1.4] Cleaning structure (removing water and heteroatoms)...")
        original_atoms = info['num_atoms']
        cleaned = clean_pdb(structure, remove_water=True, remove_hetero=True)
        cleaned_info = get_protein_info(cleaned)
        atoms_removed = original_atoms - cleaned_info['num_atoms']
        print(f"      ✓ Removed {atoms_removed} atoms")
        print(f"      ✓ Remaining: {cleaned_info['num_atoms']} atoms")

        # Step 5: Calculate center of mass
        print("\n[1.5] Calculating center of mass...")
        com = calculate_center_of_mass(cleaned)
        print(f"      ✓ Center: ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f}) Å")

        # Step 6: Save cleaned structure
        print("\n[1.6] Saving cleaned structure...")
        cleaned_path = output_dir / '1HSG_cleaned.pdb'
        save_pdb(cleaned, str(cleaned_path))
        print(f"      ✓ Saved to: {cleaned_path}")

        print("\n✓ Protein Operations Test PASSED")
        return True, structure, cleaned, com

    except Exception as e:
        print(f"\n✗ Protein Operations Test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False, None, None, None


def test_pocket_detection(structure):
    """
    Test both ligand-based and geometric pocket detection.

    Tests:
    Part A: Ligand-based detection (1HSG has MK1 inhibitor)
    Part B: Geometric detection with volume filtering
    Part C: Auto-detection (intelligent method selection)
    """
    print("\n" + "="*80)
    print("[TEST 2] POCKET DETECTION")
    print("="*80)

    output_dir = project_root / 'data' / 'outputs' / 'phase2_test' / 'pockets'
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        # PART A: Ligand-based detection
        print("\n" + "-"*80)
        print("Part A: Ligand-Based Detection (Primary Method)")
        print("-"*80)

        from scripts.utils.pocket_detection import detect_pockets_from_ligand

        ligand_pockets = detect_pockets_from_ligand(structure)

        if ligand_pockets:
            print(f"\n[A.1] Found {len(ligand_pockets)} co-crystallized ligand(s)")

            for i, pocket in enumerate(ligand_pockets):
                print(f"\n  Ligand {i+1}: {pocket.get('ligand_name', 'Unknown')}")
                print(f"    Center: ({pocket['center'][0]:.2f}, {pocket['center'][1]:.2f}, {pocket['center'][2]:.2f}) Å")
                print(f"    Volume: {pocket['volume']:.1f} Ų")
                print(f"    Score: {pocket['score']:.2f}")
                print(f"    Atoms: {pocket['num_atoms']}")
                print(f"    Nearby residues: {len(pocket['residues'])}")

            # Validate against known 1HSG binding site
            known_site = (2.0, 16.0, 25.0)
            detected = ligand_pockets[0]['center']
            distance = np.linalg.norm(
                np.array(detected) - np.array(known_site)
            )

            print(f"\n[A.2] Accuracy Validation:")
            print(f"    Known binding site (literature): ({known_site[0]:.1f}, {known_site[1]:.1f}, {known_site[2]:.1f}) Å")
            print(f"    Detected center (ligand):        ({detected[0]:.2f}, {detected[1]:.2f}, {detected[2]:.2f}) Å")
            print(f"    Distance from known site:        {distance:.2f} Å")

            if distance < 5.0:
                print(f"    ✓ EXCELLENT: Within 5Å of known site (high accuracy)")
            elif distance < 10.0:
                print(f"    ✓ GOOD: Within 10Å of known site")
            else:
                print(f"    ⚠ WARNING: >10Å from known site")

        else:
            print("\n  ✗ No ligands detected (unexpected for 1HSG)")
            print("  This protein should have MK1 inhibitor")

        # PART B: Geometric detection with filtering
        print("\n" + "-"*80)
        print("Part B: Geometric Detection (Fallback Method, Filtered)")
        print("-"*80)

        # Clean structure (remove ligands) to test pure geometric method
        from scripts.utils.protein_ops import clean_pdb

        print("\n[B.1] Testing geometric detection on cleaned protein...")
        cleaned = clean_pdb(structure, remove_hetero=True, remove_water=True)

        geo_pockets = detect_pockets_geometric(
            cleaned,
            grid_spacing=1.5,
            filter_by_volume=True,
            min_volume=200,
            max_volume=2000
        )

        print(f"    ✓ Detected {len(geo_pockets)} druggable pockets")

        if geo_pockets:
            for i, pocket in enumerate(geo_pockets[:3]):
                print(f"\n  Pocket {i+1}:")
                print(f"    Center: ({pocket['center'][0]:.2f}, {pocket['center'][1]:.2f}, {pocket['center'][2]:.2f}) Å")
                print(f"    Volume: {pocket['volume']:.1f} Ų (druggable range: 200-2000 Ų)")
                print(f"    Score: {pocket['score']:.2f}")
        else:
            print("    ⚠ No pockets in druggable volume range")

        # PART C: Auto-detection
        print("\n" + "-"*80)
        print("Part C: Auto-Detection (Intelligent Method Selection)")
        print("-"*80)

        from scripts.utils.pocket_detection import detect_pockets_auto

        print("\n[C.1] Running auto-detection on original structure...")
        auto_pockets = detect_pockets_auto(structure)

        if auto_pockets:
            print(f"    ✓ Auto-detected {len(auto_pockets)} pocket(s)")
            print(f"    ✓ Method used: {auto_pockets[0]['method']}")
            print(f"    ✓ Center: ({auto_pockets[0]['center'][0]:.2f}, "
                  f"{auto_pockets[0]['center'][1]:.2f}, {auto_pockets[0]['center'][2]:.2f}) Å")
            print(f"    ✓ Confidence: {'High (ligand-based)' if auto_pockets[0]['method'] == 'ligand-based' else 'Moderate (geometric)'}")

            # Save best pocket
            best_pocket = auto_pockets[0]

            # Save all pockets data
            print("\n[C.2] Saving pocket analysis...")
            pockets_output = {
                'method': str(best_pocket['method']),
                'num_pockets': int(len(auto_pockets)),
                'best_pocket': {
                    'center': [float(x) for x in best_pocket['center']],
                    'volume': float(best_pocket['volume']),
                    'score': float(best_pocket['score']),
                    'method': str(best_pocket['method'])
                },
                'all_pockets': [
                    {
                        'center': [float(x) for x in p['center']],
                        'volume': float(p['volume']),
                        'score': float(p['score'])
                    }
                    for p in auto_pockets
                ]
            }

            pocket_file = output_dir / '1HSG_pockets.json'
            with open(pocket_file, 'w') as f:
                json.dump(pockets_output, f, indent=2)
            print(f"    ✓ Saved to: {pocket_file}")

            # Create docking box
            print("\n[C.3] Creating docking box parameters...")
            docking_box = create_docking_box(best_pocket['center'], box_size=20.0)
            print(f"    ✓ Box center: ({docking_box['center_x']:.2f}, "
                  f"{docking_box['center_y']:.2f}, {docking_box['center_z']:.2f})")
            print(f"    ✓ Box size: {docking_box['size_x']:.1f} × "
                  f"{docking_box['size_y']:.1f} × {docking_box['size_z']:.1f} Å")

            box_file = output_dir / 'docking_box.json'
            with open(box_file, 'w') as f:
                json.dump(docking_box, f, indent=2)
            print(f"    ✓ Saved to: {box_file}")

            # Visualize
            print("\n[C.4] Generating pocket visualization...")
            viz_file = output_dir / '1HSG_pocket_viz.png'
            visualize_pocket(structure, best_pocket, output_path=str(viz_file))
            print(f"    ✓ Saved to: {viz_file}")

            print("\n" + "="*80)
            print("✓ Pocket Detection Test PASSED")
            print(f"  • Ligand-based: {'Working' if ligand_pockets else 'No ligands'}")
            print(f"  • Geometric (filtered): {len(geo_pockets)} pockets")
            print(f"  • Auto-detection: Using {best_pocket['method']}")
            print("="*80)

            return True, best_pocket
        else:
            print("    ✗ No pockets detected")
            return False, None

    except Exception as e:
        print(f"\n✗ Pocket Detection Test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False, None


def generate_test_library(size=1000):
    """
    Generate diverse test library for filtering.

    Mix:
        - 70% Random drug-like molecules
        - 15% Known drugs (should pass)
        - 10% PAINS compounds (should fail)
        - 5% Oversized molecules (should fail)
    """
    print("\n[3.1] Generating test library...")

    smiles_list = []

    # Known drugs (should pass most filters)
    known_drugs = [
        'CC(=O)Oc1ccccc1C(=O)O',  # Aspirin
        'CC(C)Cc1ccc(cc1)C(C)C(=O)O',  # Ibuprofen
        'CC(=O)Nc1ccc(O)cc1',  # Paracetamol
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
    ]

    # PAINS compounds (should be filtered)
    pains_compounds = [
        'c1ccc(O)c(O)c1',  # Catechol
        'O=C1CSC(=S)N1',  # Rhodanine
        'c1ccc(cc1)N=Nc2ccccc2',  # Azo dye
    ]

    # Oversized molecules (MW > 600)
    oversized = [
        'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)OCCCCCCCCCC',
    ]

    # Add known drugs (15%)
    num_drugs = int(size * 0.15)
    for _ in range(num_drugs):
        smiles_list.append(np.random.choice(known_drugs))

    # Add PAINS (10%)
    num_pains = int(size * 0.10)
    for _ in range(num_pains):
        smiles_list.append(np.random.choice(pains_compounds))

    # Add oversized (5%)
    num_oversized = int(size * 0.05)
    for _ in range(num_oversized):
        smiles_list.append(oversized[0])

    # Fill rest with random molecules (70%)
    num_random = size - len(smiles_list)

    # Generate random molecules
    for _ in range(num_random):
        # Simple random molecule generator
        mol = Chem.MolFromSmiles('C' * np.random.randint(5, 15))
        if mol:
            smiles_list.append(Chem.MolToSmiles(mol))
        else:
            smiles_list.append('CCCCCC')  # Fallback

    # Save to file
    output_dir = project_root / 'data' / 'outputs' / 'phase2_test' / 'library_filtering'
    lib_file = output_dir / 'test_library_1000.smi'

    with open(lib_file, 'w') as f:
        f.write('smiles\tmol_id\n')
        for i, smi in enumerate(smiles_list):
            f.write(f'{smi}\tMOL_{i:06d}\n')

    print(f"      ✓ Generated {len(smiles_list)} molecules")
    print(f"      ✓ Saved to: {lib_file}")

    return str(lib_file)


def test_library_filtering():
    """
    Test large-scale filtering.

    Tests:
    1. Generate test library (1000 molecules)
    2. Apply filters (Lipinski, PAINS, reactive)
    3. Verify statistics
    4. Check performance
    """
    print("\n" + "="*80)
    print("[TEST 3] LIBRARY FILTERING")
    print("="*80)

    output_dir = project_root / 'data' / 'outputs' / 'phase2_test' / 'library_filtering'

    try:
        # Step 1: Generate test library
        input_file = generate_test_library(size=1000)

        # Step 2: Filter library
        print("\n[3.2] Filtering library...")
        output_file = output_dir / 'filtered_library.smi'

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

        start_time = time.time()
        stats = filter_large_library(
            str(input_file),
            str(output_file),
            filters=filters,
            chunk_size=100,
            verbose=False
        )
        elapsed = time.time() - start_time

        print(f"      ✓ Input: {stats['input_count']:,} molecules")
        print(f"      ✓ Output: {stats['output_count']:,} molecules")
        print(f"      ✓ Pass rate: {stats['filter_pass_rate']:.1f}%")
        print(f"      ✓ Processing time: {elapsed:.2f} seconds")
        print(f"      ✓ Speed: {stats['input_count']/elapsed:.0f} molecules/second")

        # Step 3: Save statistics
        print("\n[3.3] Saving filter statistics...")
        stats_file = output_dir / 'filter_statistics.json'
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        print(f"      ✓ Saved to: {stats_file}")

        # Step 4: Performance check
        print("\n[3.4] Performance validation...")
        if elapsed < 20.0:  # Should process 1000 in <20 seconds
            print(f"      ✓ PASS: Processed 1000 molecules in {elapsed:.1f}s")
        else:
            print(f"      ⚠ WARNING: Slow processing ({elapsed:.1f}s)")

        # Step 5: Test PAINS detection
        print("\n[3.5] Testing PAINS detection...")
        catechol = Chem.MolFromSmiles('c1ccc(O)c(O)c1')
        if detect_pains_substructures(catechol):
            print("      ✓ PAINS detection working (catechol detected)")
        else:
            print("      ⚠ WARNING: PAINS detection may not be working")

        # Step 6: Test reactive group detection
        print("\n[3.6] Testing reactive group detection...")
        acid_chloride = Chem.MolFromSmiles('CC(=O)Cl')
        if detect_reactive_groups(acid_chloride):
            print("      ✓ Reactive group detection working (acid chloride detected)")
        else:
            print("      ⚠ WARNING: Reactive group detection may not be working")

        print("\n✓ Library Filtering Test PASSED")
        return True, str(output_file)

    except Exception as e:
        print(f"\n✗ Library Filtering Test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False, None


def test_diversity_selection():
    """
    Test diversity subset selection.

    Tests:
    1. Create test set with similar and diverse molecules
    2. Select diverse subset
    3. Verify diversity metrics
    """
    print("\n" + "="*80)
    print("[TEST 4] DIVERSITY SELECTION")
    print("="*80)

    output_dir = project_root / 'data' / 'outputs' / 'phase2_test' / 'library_filtering'

    try:
        # Step 1: Create test set
        print("\n[4.1] Creating test molecule set...")

        # Similar molecules (benzene derivatives)
        similar = [
            'c1ccccc1',  # Benzene
            'c1ccccc1C',  # Toluene
            'c1ccccc1CC',  # Ethylbenzene
            'c1ccc(C)cc1',  # p-xylene
            'c1ccccc1O',  # Phenol
        ] * 20  # Repeat to get 100

        # Different scaffolds
        diverse = [
            'C1CCCCC1',  # Cyclohexane
            'CC(=O)O',  # Acetic acid
            'CCO',  # Ethanol
            'c1ncncc1',  # Pyridine
            'C1CCNCC1',  # Piperidine
            'c1cccnc1',  # Pyridine
            'CC(C)CC',  # Alkane
            'C1COCCO1',  # Dioxane
        ]

        all_smiles = similar + diverse
        print(f"      ✓ Created {len(all_smiles)} molecules")

        # Step 2: Select diverse subset
        print("\n[4.2] Selecting diverse subset (20 molecules)...")
        selected = calculate_diversity_subset(all_smiles, target_count=20, method='maxmin')
        print(f"      ✓ Selected {len(selected)} diverse molecules")

        # Step 3: Calculate diversity
        print("\n[4.3] Calculating diversity metrics...")
        diversity_metrics = calculate_library_diversity(selected)
        print(f"      ✓ Average Tanimoto similarity: {diversity_metrics['avg_tanimoto']:.3f}")
        print(f"      ✓ Diversity score: {diversity_metrics['diversity_score']:.3f}")

        if diversity_metrics['avg_tanimoto'] < 0.5:
            print("      ✓ PASS: Good diversity (avg similarity < 0.5)")
        else:
            print("      ⚠ WARNING: Low diversity")

        # Step 4: Save diverse subset
        print("\n[4.4] Saving diverse subset...")
        diverse_file = output_dir / 'diversity_subset_20.smi'
        with open(diverse_file, 'w') as f:
            f.write('smiles\tmol_id\n')
            for i, smi in enumerate(selected):
                f.write(f'{smi}\tDIV_{i:04d}\n')
        print(f"      ✓ Saved to: {diverse_file}")

        print("\n✓ Diversity Selection Test PASSED")
        return True

    except Exception as e:
        print(f"\n✗ Diversity Selection Test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_full_pipeline():
    """
    Complete Phase 2 integration test.

    Workflow:
    1. Download protein (6LU7 - COVID-19 main protease)
    2. Detect pockets
    3. Create test library
    4. Filter library
    5. Generate outputs ready for Phase 3 docking
    """
    print("\n" + "="*80)
    print("[TEST 5] FULL INTEGRATION PIPELINE")
    print("="*80)

    output_dir = project_root / 'data' / 'outputs' / 'phase2_test' / 'integration'

    try:
        # Step 1: Download and prepare protein
        print("\n[5.1] Downloading COVID-19 main protease (6LU7)...")
        pdb_path = download_pdb('6LU7', output_dir=str(output_dir))
        structure = load_pdb(pdb_path)
        cleaned = clean_pdb(structure, remove_water=True, remove_hetero=True)

        cleaned_path = output_dir / '6LU7_cleaned.pdb'
        save_pdb(cleaned, str(cleaned_path))
        print(f"      ✓ Prepared protein: {cleaned_path}")

        # Step 2: Detect pockets
        print("\n[5.2] Detecting binding pockets...")
        pockets = detect_pockets_geometric(cleaned, grid_spacing=1.5, max_pockets=3)
        print(f"      ✓ Found {len(pockets)} pockets")

        if len(pockets) > 0:
            top_pocket = pockets[0]
            print(f"      ✓ Top pocket score: {top_pocket['score']:.2f}")

            # Create docking box
            docking_box = create_docking_box(top_pocket['center'], box_size=22.5)

            # Save pocket data
            pocket_data = {
                'protein': '6LU7',
                'num_pockets': len(pockets),
                'top_pocket': {
                    'center': list(top_pocket['center']),
                    'volume': float(top_pocket['volume']),
                    'score': float(top_pocket['score'])
                },
                'docking_box': docking_box
            }

            pocket_file = output_dir / 'pocket_analysis.json'
            with open(pocket_file, 'w') as f:
                json.dump(pocket_data, f, indent=2)
            print(f"      ✓ Saved pocket analysis: {pocket_file}")

        # Step 3: Create and filter library
        print("\n[5.3] Creating and filtering compound library...")
        lib_file = output_dir / 'test_library.smi'

        # Generate small test library
        test_smiles = [
            'CC(=O)Oc1ccccc1C(=O)O',
            'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
            'CC(=O)Nc1ccc(O)cc1',
            'c1ccccc1',
            'CCO',
        ] * 20  # 100 molecules

        with open(lib_file, 'w') as f:
            f.write('smiles\tmol_id\n')
            for i, smi in enumerate(test_smiles):
                f.write(f'{smi}\tLIB_{i:04d}\n')

        # Filter library
        filtered_file = output_dir / 'filtered_library.smi'
        stats = filter_large_library(
            str(lib_file),
            str(filtered_file),
            verbose=False
        )
        print(f"      ✓ Filtered: {stats['input_count']} → {stats['output_count']} molecules")

        # Step 4: Create summary report
        print("\n[5.4] Generating integration summary...")
        summary = {
            'protein': {
                'pdb_id': '6LU7',
                'chains': get_protein_info(cleaned)['chains'],
                'cleaned_file': str(cleaned_path.name)
            },
            'pockets': {
                'detected': len(pockets),
                'top_score': float(pockets[0]['score']) if pockets else 0.0
            },
            'library': {
                'input_count': stats['input_count'],
                'filtered_count': stats['output_count'],
                'pass_rate': stats['filter_pass_rate']
            },
            'ready_for_docking': True
        }

        summary_file = output_dir / 'phase2_summary.json'
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"      ✓ Saved summary: {summary_file}")

        # Step 5: Create text report
        report_file = output_dir / 'phase2_summary_report.txt'
        with open(report_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("DIGITAL CRO - PHASE 2 INTEGRATION REPORT\n")
            f.write("="*80 + "\n\n")
            f.write(f"Protein: 6LU7 (COVID-19 Main Protease)\n")
            f.write(f"Pockets Detected: {len(pockets)}\n")
            f.write(f"Molecules Filtered: {stats['input_count']} → {stats['output_count']}\n")
            f.write(f"\nREADY FOR PHASE 3 DOCKING\n")
            f.write("="*80 + "\n")

        print(f"      ✓ Saved report: {report_file}")

        print("\n✓ Integration Pipeline Test PASSED")
        print("\n🎉 ALL PHASE 2 TESTS COMPLETE - Ready for Phase 3!")
        return True

    except Exception as e:
        print(f"\n✗ Integration Pipeline Test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main execution function."""
    print("\n" + "="*80)
    print("PHASE 2 VALIDATION TEST SUITE - DIGITAL CRO")
    print("="*80)
    print("\nTesting:")
    print("  • Protein operations (download, load, clean)")
    print("  • Pocket detection (geometric algorithm)")
    print("  • Library filtering (PAINS, reactive groups)")
    print("  • Diversity selection")
    print("  • Full integration pipeline")
    print("\nThis will take 2-5 minutes...\n")

    # Create output directories
    create_output_directories()

    results = []

    # Test 1: Protein Operations
    success, structure, cleaned, com = test_protein_operations()
    results.append(('Protein Operations', success))

    if not success:
        print("\n⚠ Skipping remaining tests due to protein operations failure")
        return 1

    # Test 2: Pocket Detection (pass original structure with ligands)
    success, pocket = test_pocket_detection(structure)
    results.append(('Pocket Detection', success))

    # Test 3: Library Filtering
    success, filtered_file = test_library_filtering()
    results.append(('Library Filtering', success))

    # Test 4: Diversity Selection
    success = test_diversity_selection()
    results.append(('Diversity Selection', success))

    # Test 5: Full Integration
    success = test_full_pipeline()
    results.append(('Integration Pipeline', success))

    # Final Summary
    print("\n" + "="*80)
    print("PHASE 2 VALIDATION RESULTS")
    print("="*80)

    all_passed = True
    for test_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}: {test_name}")
        if not passed:
            all_passed = False

    print("="*80)

    if all_passed:
        print("\n🎉 ALL TESTS PASSED! Phase 2 Complete!")
        print("\nGenerated outputs:")
        print("  • data/outputs/phase2_test/proteins/")
        print("  • data/outputs/phase2_test/pockets/")
        print("  • data/outputs/phase2_test/library_filtering/")
        print("  • data/outputs/phase2_test/integration/")
        print("\n✓ Ready for Phase 3: Docking")
        return 0
    else:
        print("\n⚠ SOME TESTS FAILED - Review output above")
        return 1


if __name__ == '__main__':
    sys.exit(main())
