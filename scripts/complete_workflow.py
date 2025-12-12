"""
Complete Drug Discovery Workflow

End-to-end pipeline from protein PDB ID to final PDF report.

Usage:
    python scripts/complete_workflow.py --pdb 1HSG --library data/test_library.smi --output results/
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
from datetime import datetime
import logging

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.utils.protein_ops import download_pdb, load_pdb, clean_pdb
from scripts.utils.pocket_detection import detect_pockets_auto
from scripts.utils.receptor_prep import prepare_receptor_for_docking
from scripts.utils.ligand_prep import batch_prepare_ligands
from scripts.utils.vina_wrapper import get_vina_box_from_pocket
from scripts.utils.docking_batch import run_batch_docking
from scripts.utils.result_analysis import rank_docking_results, filter_hits_by_score
from scripts.utils.admet_predictions import batch_admet_analysis
from scripts.utils.molecular_viz import (
    visualize_docking_results,
    plot_admet_summary,
    create_admet_radar_plot
)
from scripts.utils.report_generator import generate_complete_report

from Bio.PDB import PDBIO

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_complete_workflow(
    pdb_id: str,
    ligand_library_path: str,
    output_dir: str,
    project_name: str = "Drug Discovery Project",
    client_name: str = "Client",
    affinity_threshold: float = -7.0,
    max_workers: int = 4,
    exhaustiveness: int = 8,
    use_consensus: bool = False
):
    """
    Run complete drug discovery workflow.

    Steps:
        1. Download and prepare protein
        2. Detect binding pocket
        3. Prepare receptor for docking
        4. Prepare ligand library
        5. Run batch docking
        6. Analyze results and rank hits
        7. Run ADMET predictions
        8. Generate visualizations
        9. Create PDF report

    Args:
        pdb_id: PDB ID of target protein
        ligand_library_path: SMILES file with tab-separated format (SMILES<tab>ID)
        output_dir: Output directory
        project_name: Project name for report
        client_name: Client name for report
        affinity_threshold: Binding affinity threshold for hits
        max_workers: Number of parallel workers for docking
        exhaustiveness: Vina exhaustiveness parameter

    Returns:
        Path to generated PDF report
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create subdirectories
    protein_dir = output_dir / 'proteins'
    ligand_dir = output_dir / 'ligands'
    docking_dir = output_dir / 'docking'
    viz_dir = output_dir / 'visualizations'

    for d in [protein_dir, ligand_dir, docking_dir, viz_dir]:
        d.mkdir(exist_ok=True)

    logger.info("="*70)
    logger.info(f"DIGITAL CRO - COMPLETE WORKFLOW")
    logger.info(f"Target: {pdb_id}")
    logger.info(f"Output: {output_dir}")
    logger.info("="*70)

    # ========================================================================
    # STEP 1: Download and prepare protein
    # ========================================================================
    logger.info("\n[1/9] Downloading protein structure...")

    raw_pdb_path = protein_dir / f'{pdb_id}.pdb'
    clean_pdb_path = protein_dir / f'{pdb_id}_cleaned.pdb'

    # Download
    download_pdb(pdb_id, str(protein_dir))

    # Load and clean
    structure = load_pdb(str(raw_pdb_path))
    cleaned = clean_pdb(structure, remove_water=True, remove_hetero=True)

    # Save cleaned structure
    io = PDBIO()
    io.set_structure(cleaned)
    io.save(str(clean_pdb_path))

    logger.info(f"✓ Protein prepared: {clean_pdb_path}")

    # ========================================================================
    # STEP 2: Detect binding pocket
    # ========================================================================
    logger.info("\n[2/9] Detecting binding pocket...")

    pockets = detect_pockets_auto(structure)

    if not pockets or len(pockets) == 0:
        raise ValueError("Failed to detect binding pocket")

    best_pocket = pockets[0]
    center = best_pocket['center']
    volume = best_pocket['volume']

    logger.info(f"✓ Pocket detected:")
    logger.info(f"  Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
    logger.info(f"  Volume: {volume:.1f} ų")

    # ========================================================================
    # STEP 3: Prepare receptor
    # ========================================================================
    logger.info("\n[3/9] Preparing receptor for docking...")

    receptor_pdbqt = protein_dir / f'{pdb_id}_receptor.pdbqt'

    rec_result = prepare_receptor_for_docking(
        pdb_path=str(clean_pdb_path),
        output_path=str(receptor_pdbqt)
    )

    logger.info(f"✓ Receptor prepared: {receptor_pdbqt}")
    logger.info(f"  Atoms: {rec_result['num_atoms']}")

    # ========================================================================
    # STEP 4: Prepare ligand library
    # ========================================================================
    logger.info("\n[4/9] Preparing ligand library...")

    # Prepare ligands using batch function
    prep_result = batch_prepare_ligands(
        smiles_file=ligand_library_path,
        output_dir=str(ligand_dir),
        max_workers=max_workers,
        verbose=True
    )

    logger.info(f"✓ Prepared {prep_result['success']}/{prep_result['total']} ligands")

    # Build ligand list for docking
    ligands = []

    # Detect format (CSV or SMI) and column indices
    with open(ligand_library_path, 'r') as f:
        first_line = f.readline().strip()
        has_comma = ',' in first_line
        is_header = False
        smiles_col_idx = 0
        id_col_idx = 1

        if any(keyword in first_line.lower() for keyword in ['smiles', 'id', 'name', 'mol_id']):
            is_header = True
            # Parse header to find column indices
            header_parts = first_line.lower().split(',' if has_comma else '\t')
            for idx, col_name in enumerate(header_parts):
                col_name = col_name.strip()
                if col_name in ['smiles', 'smile', 'smi', 'structure']:
                    smiles_col_idx = idx
                elif col_name in ['id', 'name', 'mol_id', 'compound_id', 'molecule_id']:
                    id_col_idx = idx
        else:
            # No header - detect from first data line
            parts = first_line.split(',' if has_comma else '\t')
            if len(parts) >= 2:
                from rdkit import Chem
                mol_0 = Chem.MolFromSmiles(parts[0].strip())
                mol_1 = Chem.MolFromSmiles(parts[1].strip())

                if mol_0 is not None and mol_1 is None:
                    smiles_col_idx = 0
                    id_col_idx = 1
                elif mol_1 is not None and mol_0 is None:
                    smiles_col_idx = 1
                    id_col_idx = 0

    with open(ligand_library_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#') or not line.strip():
                continue

            # Skip header line
            if line_num == 1 and is_header:
                continue

            # Split by comma if CSV, otherwise by tab
            parts = line.strip().split(',' if has_comma else '\t')

            if len(parts) >= 2:
                smiles = parts[smiles_col_idx].strip()
                mol_id = parts[id_col_idx].strip()

                # Skip empty, NaN, or None values
                if not smiles or smiles.lower() in ['nan', 'none', '', 'null']:
                    continue
                if not mol_id or mol_id.lower() in ['nan', 'none', '', 'null']:
                    continue

                pdbqt_path = ligand_dir / f"{mol_id}.pdbqt"
                if pdbqt_path.exists():
                    ligands.append({
                        'id': mol_id,
                        'smiles': smiles,
                        'pdbqt_path': str(pdbqt_path)
                    })

    logger.info(f"  {len(ligands)} ligands ready for docking")

    # Check if any ligands were prepared
    if len(ligands) == 0:
        error_msg = (
            f"\n❌ ERROR: Failed to prepare any ligands from library.\n\n"
            f"This usually means:\n"
            f"1. Invalid SMILES format in library file\n"
            f"2. Molecules too small/simple (< 3 heavy atoms)\n"
            f"3. RDKit cannot generate 3D coordinates\n"
            f"4. Library file format is incorrect\n\n"
            f"Expected format: SMILES<tab>ID (no headers)\n"
            f"Example: CCO\tethanol\n\n"
            f"Check the error messages above for specific failures.\n"
            f"Library file: {ligand_library_path}\n"
        )
        logger.error(error_msg)
        raise ValueError(error_msg)

    # ========================================================================
    # STEP 5: Run batch docking
    # ========================================================================

    # Get docking box from pocket
    docking_box = get_vina_box_from_pocket(best_pocket)

    if use_consensus:
        logger.info("\n[5/9] Running consensus molecular docking (multi-method)...")
        logger.info(f"  This will take 2-3x longer for {len(ligands)} molecules...")

        from scripts.utils.consensus_docking import ConsensusDocking
        from scripts.utils.consensus_viz import (
            plot_method_agreement,
            plot_confidence_distribution,
            plot_consensus_summary
        )

        consensus = ConsensusDocking()

        # Prepare ligand data for consensus docking
        ligand_pdbqts = [lig['pdbqt_path'] for lig in ligands]
        ligand_ids = [lig['id'] for lig in ligands]

        # Run consensus docking
        import time
        start_time = time.time()

        df_consensus = consensus.batch_consensus_dock(
            receptor_pdbqt=str(receptor_pdbqt),
            ligand_pdbqts=ligand_pdbqts,
            ligand_ids=ligand_ids,
            center=tuple(docking_box['center']),
            size=tuple(docking_box['size']),
            output_dir=str(docking_dir / "consensus"),
            max_workers=max_workers
        )

        elapsed_time = time.time() - start_time

        logger.info(f"\n✓ Consensus docking complete:")
        logger.info(f"  Successful: {len(df_consensus)}/{len(ligands)}")
        logger.info(f"  Total time: {elapsed_time/60:.1f} minutes")
        logger.info(f"  Best consensus score: {df_consensus['consensus_score'].min():.2f} kcal/mol")

        # Analyze consensus
        consensus_analysis = consensus.analyze_consensus(df_consensus)
        logger.info(f"  High confidence hits: {consensus_analysis['high_confidence']}")
        logger.info(f"  Low confidence hits: {consensus_analysis['low_confidence']}")

        # Merge with SMILES data
        library_df = pd.DataFrame(ligands)
        library_df = library_df.rename(columns={'id': 'ligand_id'})

        df_docking = library_df.merge(df_consensus, on='ligand_id', how='inner')

        # Rename consensus_score to binding_affinity for compatibility
        df_docking['binding_affinity'] = df_docking['consensus_score']
        df_docking['success'] = True

        # Generate consensus visualizations
        logger.info("  Generating consensus visualizations...")

        # Get available methods
        method_cols = [col.replace('_score', '') for col in df_consensus.columns if col.endswith('_score') and col != 'consensus_score']

        if len(method_cols) >= 2:
            # Method agreement plot
            agreement_plot = viz_dir / 'consensus_method_agreement.png'
            plot_method_agreement(
                df_consensus,
                methods=method_cols,
                output_path=str(agreement_plot),
                top_n=20
            )
            logger.info(f"    ✓ Method agreement: {agreement_plot}")

        # Confidence distribution
        confidence_plot = viz_dir / 'consensus_confidence.png'
        plot_confidence_distribution(
            df_consensus,
            output_path=str(confidence_plot)
        )
        logger.info(f"    ✓ Confidence analysis: {confidence_plot}")

        # Summary plot
        summary_plot = viz_dir / 'consensus_summary.png'
        plot_consensus_summary(
            df_consensus,
            consensus_analysis,
            output_path=str(summary_plot)
        )
        logger.info(f"    ✓ Consensus summary: {summary_plot}")

        # Save consensus results
        consensus_csv = docking_dir / 'consensus_results.csv'
        df_consensus.to_csv(consensus_csv, index=False)
        logger.info(f"    ✓ Consensus data: {consensus_csv}")

    else:
        logger.info("\n[5/9] Running batch molecular docking...")
        logger.info(f"  This may take several minutes for {len(ligands)} molecules...")

        # Run batch docking
        dock_result = run_batch_docking(
            receptor_pdbqt=str(receptor_pdbqt),
            ligands=ligands,
            output_dir=str(docking_dir),
            docking_box=docking_box,
            max_workers=max_workers,
            exhaustiveness=exhaustiveness
        )

        logger.info(f"\n✓ Docking complete:")
        logger.info(f"  Successful: {dock_result['successful']}/{dock_result['total_ligands']}")
        logger.info(f"  Total time: {dock_result['time_elapsed']/60:.1f} minutes")

        # Load results
        df_docking = pd.read_csv(dock_result['results_file'])

        # Filter out failed dockings
        df_docking = df_docking[df_docking['success'] == True].copy()

        if len(df_docking) == 0:
            raise ValueError("No successful docking results")

        logger.info(f"  Best affinity: {df_docking['binding_affinity'].min():.2f} kcal/mol")

    # ========================================================================
    # STEP 6: Rank and filter results
    # ========================================================================
    logger.info("\n[6/9] Ranking and filtering results...")

    # Rank results
    df_ranked = rank_docking_results(dock_result['results_file'], top_n=100)

    # Filter hits
    df_hits = filter_hits_by_score(df_ranked, threshold=affinity_threshold)

    logger.info(f"✓ {len(df_hits)} hits identified (threshold: {affinity_threshold} kcal/mol)")

    # ========================================================================
    # STEP 7: ADMET predictions
    # ========================================================================
    logger.info("\n[7/9] Running ADMET predictions...")

    df_final = batch_admet_analysis(df_ranked, smiles_column='smiles', id_column='ligand_id')

    lipinski_pass = df_final['lipinski_compliant'].sum()
    logger.info(f"✓ ADMET analysis complete")
    logger.info(f"  Lipinski compliant: {lipinski_pass}/{len(df_final)}")

    # Save final results
    final_output = output_dir / 'final_results.csv'
    df_final.to_csv(final_output, index=False)
    logger.info(f"  Results saved: {final_output}")

    # ========================================================================
    # STEP 8: Generate visualizations
    # ========================================================================
    logger.info("\n[8/9] Generating visualizations...")

    # Top hits visualization
    top_hits_path = viz_dir / 'top_hits.png'
    visualize_docking_results(
        df_final.sort_values('binding_affinity').head(20),
        output_path=str(top_hits_path)
    )
    logger.info(f"  ✓ Top hits: {top_hits_path}")

    # ADMET summary
    admet_summary_path = viz_dir / 'admet_summary.png'
    plot_admet_summary(df_final, output_path=str(admet_summary_path))
    logger.info(f"  ✓ ADMET summary: {admet_summary_path}")

    # Radar plot for top molecule
    top_molecule = df_final.sort_values('binding_affinity').iloc[0]
    radar_path = viz_dir / f"admet_radar_{top_molecule['ligand_id']}.png"
    create_admet_radar_plot(
        df_final,
        top_molecule['ligand_id'],
        id_column='ligand_id',
        output_path=str(radar_path)
    )
    logger.info(f"  ✓ Radar plot: {radar_path}")

    # ========================================================================
    # STEP 9: Generate PDF report
    # ========================================================================
    logger.info("\n[9/9] Generating PDF report...")

    report_path = output_dir / f'{pdb_id}_report.pdf'

    pdf_path = generate_complete_report(
        df_results=df_final.sort_values('binding_affinity'),
        output_path=str(report_path),
        project_name=project_name,
        client_name=client_name,
        target_protein=f"PDB {pdb_id}",
        target_pdb_id=pdb_id,
        visualization_dir=str(viz_dir),
        affinity_threshold=affinity_threshold
    )

    logger.info(f"✓ PDF report generated: {pdf_path}")

    # ========================================================================
    # Summary
    # ========================================================================
    logger.info("\n" + "="*70)
    logger.info("✅ WORKFLOW COMPLETE!")
    logger.info("="*70)
    logger.info(f"Molecules screened: {len(df_final)}")
    logger.info(f"Hits identified: {len(df_hits)}")
    logger.info(f"Best binding affinity: {df_final['binding_affinity'].min():.2f} kcal/mol")
    logger.info(f"Lipinski compliant: {lipinski_pass}/{len(df_final)}")
    logger.info(f"\nFinal report: {pdf_path}")
    logger.info(f"Full results: {final_output}")
    logger.info("="*70)

    return pdf_path


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description='Complete Drug Discovery Workflow - Digital CRO Platform',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python scripts/complete_workflow.py --pdb 1HSG --library data/test_library.smi

  # With custom settings
  python scripts/complete_workflow.py \\
    --pdb 1HSG \\
    --library data/molecules.smi \\
    --output results/project1 \\
    --project "HIV Protease Inhibitors" \\
    --client "Biotech Corp" \\
    --threshold -7.5 \\
    --workers 8

  # Run full workflow with all options
  python scripts/complete_workflow.py \\
    --pdb 1HSG \\
    --library data/large_library.smi \\
    --output results/hiv_screen \\
    --project "HIV-1 Protease Screening Campaign" \\
    --client "PharmaCo Research" \\
    --threshold -7.0 \\
    --workers 8 \\
    --exhaustiveness 16
        """
    )

    parser.add_argument(
        '--pdb',
        type=str,
        required=True,
        help='PDB ID of target protein (e.g., 1HSG, 3CLN, etc.)'
    )

    parser.add_argument(
        '--library',
        type=str,
        required=True,
        help='Path to ligand library file (SMILES<tab>ID format)'
    )

    parser.add_argument(
        '--output',
        type=str,
        default='data/outputs/workflow',
        help='Output directory (default: data/outputs/workflow)'
    )

    parser.add_argument(
        '--project',
        type=str,
        default='Drug Discovery Project',
        help='Project name for report (default: Drug Discovery Project)'
    )

    parser.add_argument(
        '--client',
        type=str,
        default='Client',
        help='Client name for report (default: Client)'
    )

    parser.add_argument(
        '--threshold',
        type=float,
        default=-7.0,
        help='Binding affinity threshold for hits in kcal/mol (default: -7.0)'
    )

    parser.add_argument(
        '--workers',
        type=int,
        default=4,
        help='Number of parallel workers for docking (default: 4)'
    )

    parser.add_argument(
        '--exhaustiveness',
        type=int,
        default=8,
        help='Vina exhaustiveness parameter (default: 8, range: 1-32)'
    )

    parser.add_argument(
        '--consensus',
        action='store_true',
        help='Use consensus docking with multiple methods (Vina + Smina) for higher confidence'
    )

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.library).exists():
        logger.error(f"Ligand library not found: {args.library}")
        sys.exit(1)

    # Run workflow
    try:
        run_complete_workflow(
            pdb_id=args.pdb,
            ligand_library_path=args.library,
            output_dir=args.output,
            project_name=args.project,
            client_name=args.client,
            affinity_threshold=args.threshold,
            max_workers=args.workers,
            exhaustiveness=args.exhaustiveness,
            use_consensus=args.consensus
        )
    except Exception as e:
        logger.error(f"Workflow failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
