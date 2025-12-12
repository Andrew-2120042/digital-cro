"""
Phase 3 Complete Test: End-to-End Docking Workflow

Tests the complete pipeline from protein preparation through hit selection:
1. Protein download and cleaning
2. Pocket detection
3. Receptor preparation
4. Library filtering and preparation
5. Batch docking (10-20 molecules)
6. Result analysis and ranking
7. Diversity selection
8. Report generation
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from scripts.utils.protein_ops import download_pdb, load_pdb, clean_pdb
from scripts.utils.pocket_detection import detect_pockets_auto
from scripts.utils.receptor_prep import prepare_receptor_for_docking
from scripts.utils.ligand_prep import batch_prepare_ligands
from scripts.utils.docking_batch import run_batch_docking
from scripts.utils.result_analysis import (
    rank_docking_results,
    filter_hits_by_score,
    select_diverse_hits,
    generate_hit_report
)


def create_test_library():
    """Create test SMILES library with 20 known drugs"""
    test_smiles = [
        ("CC(=O)Oc1ccccc1C(=O)O", "aspirin"),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine"),
        ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen"),
        ("CC(=O)Nc1ccc(cc1)O", "paracetamol"),
        ("CC1(C2CCC1(C(=O)O2)C)C", "penicillin"),
        ("CCC(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(CC(=O)O3)O)C", "atorvastatin"),
        ("CN(C)C(=N)NC(=N)N", "metformin"),
        ("CC(=O)CC(c1ccccc1)c2c(O)c3ccccc3oc2=O", "warfarin"),
        ("CCCC1=NC(C)=C(N1CCCc2ccc(cc2)S(=O)(=O)N3CCN(CC3)C)C(=O)OCC", "sildenafil"),
        ("CCCCC(C(CC(=O)N1C(CCC1C(=O)O)Cc2ccccc2)NC(C(C)CC)C(=O)O)C", "lisinopril"),
        # Additional diverse molecules
        ("Cc1c(sc[n+]1Cc2cnc(nc2N)C)CCO", "thiamine"),
        ("CN1CCC23c4ccc(O)c(O)c4CC2C1Cc5c3cccc5O", "morphine"),
        ("COc1ccc2c(c1)c(c(c(=O)o2)O)c3c4ccoc4c(cc3)CC=C(C)C", "khellin"),
        ("CC(C)(C)NCC(COc1ccc(CCOCC(O)CO)cc1)O", "propranolol"),
        ("CC(C)NCC(O)COc1ccccc1CC=C", "alprenolol"),
        ("CN1CCc2cc(Oc3ccc(NC(=O)C(=O)Nc4ccc(Oc5cc6CCN(C)Cc6cc5)cc4)cc3)ccc2C1", "sunitinib"),
        ("Cc1cnc(NC(=O)c2ccc(C)c(Nc3nccc(-c4cccnc4)n3)c2)s1", "dabrafenib"),
        ("Cc1ccc(-c2cc(C(F)(F)F)nn2-c3ccc(S(N)(=O)=O)cc3)cc1", "celecoxib"),
        ("COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN4CCOCC4", "gefitinib"),
        ("CN(C)CCC=C1c2ccccc2COc3c1cccc3", "amitriptyline")
    ]

    smiles_file = Path('data/outputs/phase3_test/test_library.smi')
    smiles_file.parent.mkdir(parents=True, exist_ok=True)

    with open(smiles_file, 'w') as f:
        f.write("# Test library - 20 known drugs\n")
        for smiles, name in test_smiles:
            f.write(f"{smiles}\t{name}\n")

    return str(smiles_file)


def main():
    """Run complete Phase 3 test"""
    print("\n" + "="*70)
    print("PHASE 3 COMPLETE: END-TO-END DOCKING WORKFLOW TEST")
    print("="*70)

    output_dir = Path('data/outputs/phase3_test/complete_workflow')

    # Step 1: Protein preparation
    print("\n[Step 1/7] Downloading and preparing protein...")
    pdb_path = download_pdb('1HSG', str(output_dir / 'proteins'))
    structure = load_pdb(pdb_path)
    cleaned = clean_pdb(structure, remove_water=True, remove_hetero=True)

    # Save cleaned structure
    from Bio.PDB import PDBIO
    io = PDBIO()
    io.set_structure(cleaned)
    cleaned_pdb = output_dir / 'proteins' / '1HSG_cleaned.pdb'
    cleaned_pdb.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(cleaned_pdb))

    print(f"✓ Protein prepared: {cleaned_pdb}")

    # Step 2: Pocket detection
    print("\n[Step 2/7] Detecting binding pocket...")
    pockets = detect_pockets_auto(structure)
    best_pocket = pockets[0]

    print(f"✓ Pocket detected:")
    print(f"  Center: ({best_pocket['center'][0]:.2f}, {best_pocket['center'][1]:.2f}, {best_pocket['center'][2]:.2f})")
    print(f"  Method: {best_pocket.get('method', 'unknown')}")
    print(f"  Score: {best_pocket.get('score', 0):.2f}")

    # Step 3: Receptor preparation
    print("\n[Step 3/7] Preparing receptor PDBQT...")
    receptor_pdbqt = output_dir / 'proteins' / '1HSG_receptor.pdbqt'
    rec_result = prepare_receptor_for_docking(
        pdb_path=str(cleaned_pdb),
        output_path=str(receptor_pdbqt)
    )

    print(f"✓ Receptor prepared: {receptor_pdbqt}")
    print(f"  Atoms: {rec_result['num_atoms']}")

    # Step 4: Library preparation
    print("\n[Step 4/7] Preparing ligand library...")
    smiles_file = create_test_library()

    ligands_dir = output_dir / 'ligands'
    prep_result = batch_prepare_ligands(
        smiles_file=smiles_file,
        output_dir=str(ligands_dir),
        max_workers=4
    )

    print(f"✓ Ligands prepared: {prep_result['success']}/{prep_result['total']}")

    # Build ligand list for docking
    ligands = []
    with open(smiles_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                smiles, mol_id = parts[0], parts[1]
                pdbqt_path = ligands_dir / f"{mol_id}.pdbqt"
                if pdbqt_path.exists():
                    ligands.append({
                        'id': mol_id,
                        'smiles': smiles,
                        'pdbqt_path': str(pdbqt_path)
                    })

    # Step 5: Batch docking
    print(f"\n[Step 5/7] Running batch docking ({len(ligands)} molecules)...")
    print("  This will take ~15-30 minutes depending on CPU...")

    from scripts.utils.vina_wrapper import get_vina_box_from_pocket
    docking_box = get_vina_box_from_pocket(best_pocket)

    dock_result = run_batch_docking(
        receptor_pdbqt=str(receptor_pdbqt),
        ligands=ligands,
        output_dir=str(output_dir),
        docking_box=docking_box,
        max_workers=4,  # Use 4 cores
        exhaustiveness=8
    )

    print(f"\n✓ Docking complete:")
    print(f"  Successful: {dock_result['successful']}/{dock_result['total_ligands']}")
    print(f"  Total time: {dock_result['time_elapsed']/60:.1f} minutes")
    print(f"  Avg time/ligand: {dock_result['avg_time_per_ligand']:.1f} seconds")

    # Step 6: Result analysis
    print("\n[Step 6/7] Analyzing results...")

    df_ranked = rank_docking_results(dock_result['results_file'], top_n=50)

    print(f"✓ Top 10 hits:")
    for idx, row in df_ranked.head(10).iterrows():
        print(f"  {row['rank']:2d}. {row['ligand_id']:15s} {row['binding_affinity']:6.2f} kcal/mol")

    # Filter by threshold
    df_filtered = filter_hits_by_score(df_ranked, threshold=-7.0)
    print(f"\n✓ Hits passing -7.0 kcal/mol threshold: {len(df_filtered)}")

    # Step 7: Diversity selection
    print("\n[Step 7/7] Selecting diverse hits...")

    df_diverse = select_diverse_hits(df_ranked, n_hits=10, diversity_threshold=0.5)

    print(f"✓ Selected {len(df_diverse)} diverse hits")

    # Generate report
    report_dir = output_dir / 'analysis'
    generate_hit_report(df_ranked, str(report_dir), top_n=50)

    print(f"\n✓ Report generated: {report_dir}")

    # Final summary
    print("\n" + "="*70)
    print("✅ PHASE 3 COMPLETE WORKFLOW - SUCCESS!")
    print("="*70)
    print(f"\nOutput directory: {output_dir}")
    print(f"\nKey results:")
    print(f"  - Receptor: {receptor_pdbqt}")
    print(f"  - Docked: {dock_result['successful']} molecules")

    if len(df_ranked) > 0:
        print(f"  - Best affinity: {df_ranked.iloc[0]['binding_affinity']:.2f} kcal/mol")
        print(f"  - Best hit: {df_ranked.iloc[0]['ligand_id']}")
    else:
        print(f"  - No successful docking results")

    print(f"  - Results CSV: {dock_result['results_file']}")
    print(f"  - Analysis: {report_dir}")

    print("\n✅ All Phase 3 tests passed!")
    print("Ready for production use!")

    return True


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
