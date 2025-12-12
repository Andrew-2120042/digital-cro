"""
Batch docking module for Digital CRO.

Processes multiple ligands in parallel using AutoDock Vina.

Key features:
    - Parallel execution (ThreadPoolExecutor)
    - Progress tracking (tqdm)
    - Robust error handling
    - Result aggregation
    - Memory-efficient processing
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
import pandas as pd
from datetime import datetime
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import logging

from .vina_wrapper import dock_single_ligand
from .ligand_prep import prepare_ligand_for_docking

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def dock_ligand_worker(args: Tuple) -> Dict:
    """
    Worker function for parallel docking.

    Args:
        args: Tuple of (receptor_pdbqt, ligand_info, output_dir, box_params, vina_params)

    Returns:
        dict with docking results for this ligand

    This is called by ThreadPoolExecutor for each ligand.
    """
    receptor_pdbqt, ligand_info, output_dir, box_params, vina_params = args

    ligand_id = ligand_info.get('id', 'unknown')
    smiles = ligand_info.get('smiles', '')
    ligand_pdbqt = ligand_info.get('pdbqt_path', '')

    try:
        # Prepare output path
        output_pdbqt = Path(output_dir) / 'results' / f'{ligand_id}_docked.pdbqt'
        output_pdbqt.parent.mkdir(parents=True, exist_ok=True)

        # Run docking
        start_time = time.time()

        result = dock_single_ligand(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            output_pdbqt=str(output_pdbqt),
            center=(box_params['center_x'], box_params['center_y'], box_params['center_z']),
            size=(box_params['size_x'], box_params['size_y'], box_params['size_z']),
            exhaustiveness=vina_params.get('exhaustiveness', 8),
            num_modes=vina_params.get('num_modes', 9)
        )

        runtime = time.time() - start_time

        if result['success']:
            return {
                'ligand_id': ligand_id,
                'smiles': smiles,
                'binding_affinity': result['best_affinity'],
                'num_modes': len(result.get('modes', [])),
                'all_modes': result.get('modes', []),
                'output_file': str(output_pdbqt),
                'runtime': runtime,
                'success': True,
                'error': None
            }
        else:
            return {
                'ligand_id': ligand_id,
                'smiles': smiles,
                'binding_affinity': None,
                'num_modes': 0,
                'success': False,
                'error': result.get('error', 'Unknown error'),
                'runtime': runtime
            }

    except Exception as e:
        logger.error(f"Error docking {ligand_id}: {e}")
        return {
            'ligand_id': ligand_id,
            'smiles': smiles,
            'binding_affinity': None,
            'num_modes': 0,
            'success': False,
            'error': str(e),
            'runtime': 0
        }


def run_batch_docking(
    receptor_pdbqt: str,
    ligands: List[Dict],
    output_dir: str,
    docking_box: Dict,
    max_workers: int = 4,
    exhaustiveness: int = 8,
    num_modes: int = 9
) -> Dict:
    """
    Dock multiple ligands in parallel.

    Args:
        receptor_pdbqt: Path to prepared receptor .pdbqt
        ligands: List of dicts with ligand info:
                 [{'id': 'mol1', 'smiles': '...', 'pdbqt_path': '...'},  ...]
        output_dir: Base directory for all outputs
        docking_box: Dict with center_x, center_y, center_z, size_x, size_y, size_z
        max_workers: Number of parallel Vina processes (CPU cores)
        exhaustiveness: Vina search parameter
        num_modes: Number of binding modes per ligand

    Returns:
        dict with:
        {
            'total_ligands': int,
            'successful': int,
            'failed': int,
            'results_file': str (path to CSV),
            'results': list of all individual results,
            'time_elapsed': float,
            'avg_time_per_ligand': float
        }

    Example:
        ligands = [
            {'id': 'aspirin', 'smiles': 'CC(=O)Oc1ccccc1C(=O)O',
             'pdbqt_path': 'ligands/aspirin.pdbqt'},
            ...
        ]

        results = run_batch_docking(
            receptor_pdbqt='receptor.pdbqt',
            ligands=ligands,
            output_dir='project1',
            docking_box={'center_x': 10, 'center_y': 20, 'center_z': 15,
                        'size_x': 20, 'size_y': 20, 'size_z': 20},
            max_workers=8
        )
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Starting batch docking: {len(ligands)} ligands")
    logger.info(f"Parallel workers: {max_workers}")
    logger.info(f"Exhaustiveness: {exhaustiveness}")

    start_time = time.time()

    # Prepare arguments for workers
    vina_params = {
        'exhaustiveness': exhaustiveness,
        'num_modes': num_modes
    }

    worker_args = [
        (receptor_pdbqt, ligand, output_dir, docking_box, vina_params)
        for ligand in ligands
    ]

    # Run parallel docking
    results = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        futures = {executor.submit(dock_ligand_worker, args): args[1]['id']
                  for args in worker_args}

        # Process as they complete
        with tqdm(total=len(ligands), desc="Docking ligands") as pbar:
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                pbar.update(1)

                # Show current best
                successful_affinities = [r['binding_affinity'] for r in results if r['success'] and r['binding_affinity'] is not None]
                if successful_affinities:
                    pbar.set_postfix({'best': f"{min(successful_affinities):.2f} kcal/mol"})

    # Calculate statistics
    successful = [r for r in results if r['success']]
    failed = [r for r in results if not r['success']]

    total_time = time.time() - start_time
    avg_time = total_time / len(ligands) if ligands else 0

    logger.info(f"Batch docking complete:")
    logger.info(f"  Successful: {len(successful)}/{len(ligands)}")
    logger.info(f"  Failed: {len(failed)}/{len(ligands)}")
    logger.info(f"  Total time: {total_time/60:.1f} minutes")
    logger.info(f"  Avg time/ligand: {avg_time:.1f} seconds")

    # Save results to CSV
    results_csv = output_dir / 'docking_results.csv'
    create_docking_summary(results, str(results_csv))

    # Save detailed JSON
    results_json = output_dir / 'docking_results_detailed.json'
    with open(results_json, 'w') as f:
        json.dump({
            'metadata': {
                'receptor': receptor_pdbqt,
                'total_ligands': len(ligands),
                'successful': len(successful),
                'failed': len(failed),
                'docking_box': docking_box,
                'vina_params': vina_params,
                'timestamp': datetime.now().isoformat(),
                'runtime_seconds': total_time
            },
            'results': results
        }, f, indent=2)

    return {
        'total_ligands': len(ligands),
        'successful': len(successful),
        'failed': len(failed),
        'results_file': str(results_csv),
        'results_json': str(results_json),
        'results': results,
        'time_elapsed': total_time,
        'avg_time_per_ligand': avg_time
    }


def create_docking_summary(results: List[Dict], output_csv: str):
    """
    Create summary CSV of docking results.

    Args:
        results: List of result dicts from batch docking
        output_csv: Where to save CSV

    CSV columns:
        - ligand_id
        - smiles
        - binding_affinity (kcal/mol)
        - num_modes
        - runtime (seconds)
        - success
        - error (if failed)

    Sorted by binding affinity (best first)
    """
    df = pd.DataFrame(results)

    # Handle empty results
    if df.empty or 'success' not in df.columns:
        logger.warning("No docking results to summarize")
        # Create empty CSV with correct columns
        empty_df = pd.DataFrame(columns=[
            'ligand_id', 'smiles', 'binding_affinity', 'num_modes',
            'runtime', 'success', 'error'
        ])
        empty_df.to_csv(output_csv, index=False)
        logger.info(f"Empty results file saved to: {output_csv}")
        return empty_df

    # Sort by binding affinity (lowest = best)
    df_success = df[df['success'] == True].copy()
    df_failed = df[df['success'] == False].copy()

    if not df_success.empty:
        df_success = df_success.sort_values('binding_affinity')

    # Combine (successful first, then failed)
    df_final = pd.concat([df_success, df_failed], ignore_index=True)

    # Select columns for CSV
    columns = ['ligand_id', 'smiles', 'binding_affinity', 'num_modes',
               'runtime', 'success', 'error']
    df_final = df_final[columns]

    # Save
    df_final.to_csv(output_csv, index=False)
    logger.info(f"Results saved to: {output_csv}")

    return df_final


def prepare_and_dock_from_smiles(
    receptor_pdbqt: str,
    smiles_file: str,
    output_dir: str,
    docking_box: Dict,
    max_workers: int = 4
) -> Dict:
    """
    Complete workflow: prepare ligands from SMILES file, then dock all.

    Args:
        receptor_pdbqt: Prepared receptor
        smiles_file: File with SMILES strings (format: "SMILES ID" per line)
        output_dir: Output directory
        docking_box: Docking box parameters
        max_workers: Parallel workers

    Returns:
        Batch docking results dict

    This combines ligand preparation (Phase 3.1) with batch docking (Phase 3.3).
    """
    from .ligand_prep import batch_prepare_ligands

    output_dir = Path(output_dir)
    ligands_dir = output_dir / 'ligands'
    ligands_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Step 1: Preparing ligands from SMILES")

    # Prepare all ligands
    prep_results = batch_prepare_ligands(
        smiles_file=smiles_file,
        output_dir=str(ligands_dir),
        max_workers=max_workers
    )

    logger.info(f"Ligand preparation: {prep_results['success']}/{prep_results['total']} successful")

    # Build ligand list for docking
    ligands = []

    # Read SMILES file to get IDs
    with open(smiles_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) >= 2:
                smiles, mol_id = parts[0], parts[1]
            else:
                smiles = parts[0]
                mol_id = f"mol_{len(ligands)}"

            pdbqt_path = ligands_dir / f"{mol_id}.pdbqt"

            if pdbqt_path.exists():
                ligands.append({
                    'id': mol_id,
                    'smiles': smiles,
                    'pdbqt_path': str(pdbqt_path)
                })

    logger.info(f"Step 2: Docking {len(ligands)} prepared ligands")

    # Run batch docking
    docking_results = run_batch_docking(
        receptor_pdbqt=receptor_pdbqt,
        ligands=ligands,
        output_dir=output_dir,
        docking_box=docking_box,
        max_workers=max_workers
    )

    return docking_results
