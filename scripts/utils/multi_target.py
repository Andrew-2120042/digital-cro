"""
Multi-Target Screening Module

Enables screening one molecule library against multiple protein targets.
Generates comparative analysis and selectivity reports.
"""

from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd
import logging
from datetime import datetime
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from scripts.complete_workflow import run_complete_workflow

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def run_multi_target_screening(
    target_proteins: List[str],
    ligand_library_path: str,
    output_dir: str,
    project_name: str = "Multi-Target Screening",
    client_name: str = "Client",
    affinity_threshold: float = -7.0,
    max_workers_per_target: int = 4,
    run_targets_parallel: bool = False
) -> Dict:
    """
    Run screening against multiple protein targets.

    Args:
        target_proteins: List of PDB IDs (e.g., ['1HSG', '3CL5', '1M17'])
        ligand_library_path: Path to molecule library CSV or SMI
        output_dir: Output directory for all results
        project_name: Project name for reports
        client_name: Client name for reports
        affinity_threshold: Binding affinity threshold
        max_workers_per_target: Workers for parallel docking per target
        run_targets_parallel: If True, run targets in parallel (uses more RAM)

    Returns:
        Dict with results summary and paths
    """
    from rdkit import Chem

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("="*70)
    logger.info(f"MULTI-TARGET SCREENING")
    logger.info(f"Targets: {len(target_proteins)}")
    logger.info(f"Library: {ligand_library_path}")
    logger.info("="*70)

    # ========== CRITICAL FIX: Proper CSV/Library Handling ==========

    logger.info(f"\n[1/3] Reading and validating library file...")
    logger.info(f"  Source: {ligand_library_path}")

    library_path_obj = Path(ligand_library_path)

    if library_path_obj.suffix.lower() == '.csv':
        # Read CSV with proper column parsing
        logger.info("  Format: CSV - parsing columns...")

        try:
            library_df = pd.read_csv(library_path_obj, dtype=str)

            # Normalize column names to lowercase
            library_df.columns = library_df.columns.str.lower().str.strip()

            logger.info(f"  Original columns: {library_df.columns.tolist()}")

            # Map common column name variations to standard names
            column_mapping = {}

            # Find SMILES column
            smiles_candidates = ['smiles', 'smile', 'smi', 'structure', 'mol']
            for col in library_df.columns:
                if col in smiles_candidates:
                    column_mapping[col] = 'smiles'
                    break

            # Find ID column
            id_candidates = ['id', 'name', 'mol_id', 'compound_id', 'molecule_id']
            for col in library_df.columns:
                if col in id_candidates and col not in column_mapping:
                    column_mapping[col] = 'id'
                    break

            # Rename columns
            if column_mapping:
                library_df = library_df.rename(columns=column_mapping)

            # Validate required columns exist
            if 'smiles' not in library_df.columns:
                # If only 2 columns and no header match, assume first is ID, second is SMILES
                if len(library_df.columns) == 2:
                    library_df.columns = ['id', 'smiles']
                    logger.info("  Applied default column mapping: first=id, second=smiles")
                else:
                    raise ValueError(f"CSV must have 'smiles' column. Found columns: {library_df.columns.tolist()}")

            if 'id' not in library_df.columns:
                # Generate IDs if missing
                library_df['id'] = [f'mol_{i:06d}' for i in range(len(library_df))]
                logger.info("  Generated molecule IDs (column missing)")

            # Keep only required columns
            library_df = library_df[['id', 'smiles']].copy()

            # Remove any rows with missing values
            library_df = library_df.dropna()

            logger.info(f"  ✓ Parsed {len(library_df)} molecules from CSV")

        except Exception as e:
            logger.error(f"  ✗ Failed to parse CSV: {e}")
            raise ValueError(f"Failed to read CSV library: {e}")

    elif library_path_obj.suffix.lower() in ['.smi', '.txt']:
        # Read SMI file (tab-separated or space-separated: SMILES<tab>ID)
        logger.info("  Format: SMI - reading tab-separated format...")

        try:
            lines = library_path_obj.read_text().strip().split('\n')

            smiles_list = []
            id_list = []

            for idx, line in enumerate(lines):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                if '\t' in line:
                    parts = line.split('\t')
                    smiles_list.append(parts[0].strip())
                    id_list.append(parts[1].strip() if len(parts) > 1 else f'mol_{idx:06d}')
                elif ' ' in line:
                    parts = line.split(None, 1)  # Split on first whitespace
                    smiles_list.append(parts[0].strip())
                    id_list.append(parts[1].strip() if len(parts) > 1 else f'mol_{idx:06d}')
                else:
                    # Just SMILES, no ID
                    smiles_list.append(line)
                    id_list.append(f'mol_{idx:06d}')

            library_df = pd.DataFrame({'id': id_list, 'smiles': smiles_list})

            logger.info(f"  ✓ Read {len(library_df)} molecules from SMI file")

        except Exception as e:
            logger.error(f"  ✗ Failed to read SMI file: {e}")
            raise ValueError(f"Failed to read SMI library: {e}")

    else:
        raise ValueError(f"Unsupported library format: {library_path_obj.suffix}. Use .csv, .smi, or .txt")

    # Validate SMILES
    logger.info(f"\n[2/3] Validating SMILES structures...")

    valid_indices = []
    invalid_count = 0

    for idx, row in library_df.iterrows():
        try:
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol is not None:
                valid_indices.append(idx)
            else:
                logger.warning(f"  ⚠️  Invalid SMILES skipped: {row['id']} - {row['smiles'][:50]}")
                invalid_count += 1
        except Exception as e:
            logger.warning(f"  ⚠️  Error parsing {row['id']}: {e}")
            invalid_count += 1

    if not valid_indices:
        raise ValueError("No valid molecules found in library! All SMILES are invalid.")

    library_df = library_df.loc[valid_indices].reset_index(drop=True)

    logger.info(f"  ✓ Validated {len(library_df)} molecules")
    if invalid_count > 0:
        logger.info(f"  ⚠️  Skipped {invalid_count} invalid SMILES")

    # Save cleaned library in .smi format (what the workflow expects)
    cleaned_library_path = output_dir / 'library_cleaned.smi'

    logger.info(f"\n[3/3] Saving cleaned library...")
    logger.info(f"  Format: .smi (tab-separated: SMILES<tab>ID)")

    with open(cleaned_library_path, 'w') as f:
        for _, row in library_df.iterrows():
            f.write(f"{row['smiles']}\t{row['id']}\n")

    logger.info(f"  ✓ Saved to: {cleaned_library_path}")

    # DEBUG: Verify library data
    logger.info("\n" + "="*70)
    logger.info("LIBRARY VERIFICATION:")
    logger.info(f"  Shape: {library_df.shape}")
    logger.info(f"  Columns: {library_df.columns.tolist()}")
    logger.info(f"  First 3 molecules:")
    for idx, row in library_df.head(3).iterrows():
        logger.info(f"    [{idx}] ID: '{row['id']}' | SMILES: '{row['smiles'][:50]}...'")
        mol = Chem.MolFromSmiles(row['smiles'])
        logger.info(f"         Valid: {mol is not None}")
    logger.info("="*70 + "\n")

    # Update library path to use cleaned version
    ligand_library_path = str(cleaned_library_path)

    # ========== END OF CRITICAL FIX ==========

    # Results storage
    target_results = {}

    # Process each target
    if run_targets_parallel:
        # Run targets in parallel (faster but uses more RAM)
        target_results = _run_targets_parallel(
            target_proteins=target_proteins,
            ligand_library_path=ligand_library_path,
            output_dir=output_dir,
            project_name=project_name,
            client_name=client_name,
            affinity_threshold=affinity_threshold,
            max_workers_per_target=max_workers_per_target
        )
    else:
        # Run targets sequentially (safer, less RAM)
        target_results = _run_targets_sequential(
            target_proteins=target_proteins,
            ligand_library_path=ligand_library_path,
            output_dir=output_dir,
            project_name=project_name,
            client_name=client_name,
            affinity_threshold=affinity_threshold,
            max_workers_per_target=max_workers_per_target
        )

    # Create comparative analysis
    logger.info("\nGenerating comparative analysis...")

    comparative_df = create_comparative_table(target_results)

    # Save comparative results
    comparative_path = output_dir / "multi_target_results.csv"
    comparative_df.to_csv(comparative_path, index=False)

    logger.info(f"✓ Comparative results saved: {comparative_path}")

    # Generate selectivity analysis
    selectivity_df = analyze_selectivity(comparative_df, target_proteins, affinity_threshold)

    selectivity_path = output_dir / "selectivity_analysis.csv"
    selectivity_df.to_csv(selectivity_path, index=False)

    logger.info(f"✓ Selectivity analysis saved: {selectivity_path}")

    # Create summary
    summary = create_multi_target_summary(target_results, comparative_df, selectivity_df)

    summary_path = output_dir / "multi_target_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info(f"✓ Summary saved: {summary_path}")

    # Print summary
    print_multi_target_summary(summary, target_proteins)

    return {
        'target_results': target_results,
        'comparative_df': comparative_df,
        'selectivity_df': selectivity_df,
        'summary': summary,
        'output_dir': str(output_dir)
    }


def _run_targets_sequential(
    target_proteins: List[str],
    ligand_library_path: str,
    output_dir: Path,
    project_name: str,
    client_name: str,
    affinity_threshold: float,
    max_workers_per_target: int
) -> Dict:
    """Run targets one by one (safer, less RAM)"""

    target_results = {}

    for idx, pdb_id in enumerate(target_proteins, 1):
        logger.info(f"\n{'='*70}")
        logger.info(f"TARGET {idx}/{len(target_proteins)}: {pdb_id}")
        logger.info(f"{'='*70}")

        # Create target-specific output directory
        target_dir = output_dir / f"target_{pdb_id}"

        try:
            # Run complete workflow for this target
            pdf_path = run_complete_workflow(
                pdb_id=pdb_id,
                ligand_library_path=ligand_library_path,
                output_dir=str(target_dir),
                project_name=f"{project_name} - {pdb_id}",
                client_name=client_name,
                affinity_threshold=affinity_threshold,
                max_workers=max_workers_per_target
            )

            # Load results
            results_csv = target_dir / "final_results.csv"

            if results_csv.exists():
                df = pd.read_csv(results_csv)

                target_results[pdb_id] = {
                    'success': True,
                    'results_df': df,
                    'results_path': str(results_csv),
                    'pdf_path': pdf_path,
                    'output_dir': str(target_dir),
                    'num_molecules': len(df),
                    'num_hits': len(df[df['binding_affinity'] <= affinity_threshold]),
                    'best_affinity': df['binding_affinity'].min()
                }

                logger.info(f"✓ {pdb_id} complete: {len(df)} molecules, {target_results[pdb_id]['num_hits']} hits")
            else:
                target_results[pdb_id] = {'success': False, 'error': 'Results file not found'}
                logger.error(f"✗ {pdb_id} failed: Results file not found")

        except Exception as e:
            target_results[pdb_id] = {'success': False, 'error': str(e)}
            logger.error(f"✗ {pdb_id} failed: {e}")

    return target_results


def _run_targets_parallel(
    target_proteins: List[str],
    ligand_library_path: str,
    output_dir: Path,
    project_name: str,
    client_name: str,
    affinity_threshold: float,
    max_workers_per_target: int
) -> Dict:
    """Run targets in parallel (faster but uses more RAM)"""

    target_results = {}

    # Use ThreadPoolExecutor to run targets in parallel
    max_parallel_targets = min(len(target_proteins), 3)  # Limit to 3 parallel targets

    with ThreadPoolExecutor(max_workers=max_parallel_targets) as executor:
        # Submit all targets
        future_to_target = {}

        for pdb_id in target_proteins:
            target_dir = output_dir / f"target_{pdb_id}"

            future = executor.submit(
                run_complete_workflow,
                pdb_id=pdb_id,
                ligand_library_path=ligand_library_path,
                output_dir=str(target_dir),
                project_name=f"{project_name} - {pdb_id}",
                client_name=client_name,
                affinity_threshold=affinity_threshold,
                max_workers=max_workers_per_target
            )

            future_to_target[future] = pdb_id

        # Collect results as they complete
        for future in as_completed(future_to_target):
            pdb_id = future_to_target[future]

            try:
                pdf_path = future.result()
                target_dir = output_dir / f"target_{pdb_id}"
                results_csv = target_dir / "final_results.csv"

                if results_csv.exists():
                    df = pd.read_csv(results_csv)

                    target_results[pdb_id] = {
                        'success': True,
                        'results_df': df,
                        'results_path': str(results_csv),
                        'pdf_path': pdf_path,
                        'output_dir': str(target_dir),
                        'num_molecules': len(df),
                        'num_hits': len(df[df['binding_affinity'] <= affinity_threshold]),
                        'best_affinity': df['binding_affinity'].min()
                    }

                    logger.info(f"✓ {pdb_id} complete")
                else:
                    target_results[pdb_id] = {'success': False, 'error': 'Results not found'}

            except Exception as e:
                target_results[pdb_id] = {'success': False, 'error': str(e)}
                logger.error(f"✗ {pdb_id} failed: {e}")

    return target_results


def create_comparative_table(target_results: Dict) -> pd.DataFrame:
    """
    Create comparative table showing all molecules across all targets.

    Returns DataFrame with columns:
    - ligand_id
    - target_X_affinity (for each target)
    - target_X_rank (for each target)
    - best_target
    - best_affinity
    - num_targets_hit
    """
    # Get all successful targets
    successful_targets = {
        pdb_id: data
        for pdb_id, data in target_results.items()
        if data.get('success', False)
    }

    if not successful_targets:
        logger.warning("No successful targets to compare")
        return pd.DataFrame()

    # Start with first target's molecules
    first_target = list(successful_targets.keys())[0]
    df_combined = successful_targets[first_target]['results_df'][['ligand_id', 'smiles']].copy()

    # Add affinity columns for each target
    for pdb_id, data in successful_targets.items():
        df_target = data['results_df'][['ligand_id', 'binding_affinity']].copy()
        df_target = df_target.rename(columns={'binding_affinity': f'{pdb_id}_affinity'})

        df_combined = df_combined.merge(df_target, on='ligand_id', how='outer')

    # Calculate best target and affinity for each molecule
    affinity_cols = [col for col in df_combined.columns if col.endswith('_affinity')]

    df_combined['best_affinity'] = df_combined[affinity_cols].min(axis=1)
    df_combined['best_target'] = df_combined[affinity_cols].idxmin(axis=1).str.replace('_affinity', '')

    # Count how many targets each molecule hits (affinity <= -7.0)
    df_combined['num_targets_hit'] = (df_combined[affinity_cols] <= -7.0).sum(axis=1)

    # Sort by best affinity
    df_combined = df_combined.sort_values('best_affinity')

    return df_combined


def analyze_selectivity(
    comparative_df: pd.DataFrame,
    target_proteins: List[str],
    affinity_threshold: float = -7.0
) -> pd.DataFrame:
    """
    Analyze selectivity - which molecules are selective for specific targets.

    Returns DataFrame with:
    - ligand_id
    - selective_for (target name or "Non-selective")
    - selectivity_score (difference between best and second-best)
    """
    selectivity_data = []

    # CRITICAL FIX: Only use affinity columns that actually exist in the dataframe
    # Some targets may have failed, so not all columns will be present
    affinity_cols = [f'{pdb_id}_affinity' for pdb_id in target_proteins]

    # Filter to only columns that exist in the dataframe
    existing_affinity_cols = [col for col in affinity_cols if col in comparative_df.columns]

    if len(existing_affinity_cols) < 2:
        logger.warning(f"Not enough successful targets for selectivity analysis. Need at least 2, found {len(existing_affinity_cols)}")
        # Return empty dataframe with expected structure
        return pd.DataFrame(columns=['ligand_id', 'best_target', 'best_affinity',
                                    'second_best_affinity', 'selectivity_score',
                                    'classification', 'selective_for'])

    for _, row in comparative_df.iterrows():
        affinities = row[existing_affinity_cols].values

        # Filter out NaN values (molecules that failed docking for some targets)
        valid_affinities = [a for a in affinities if pd.notna(a)]

        if len(valid_affinities) < 2:
            # Not enough valid affinities to calculate selectivity
            logger.warning(f"Skipping {row.get('ligand_id')} - insufficient valid affinities: {valid_affinities}")
            continue

        # Sort affinities
        sorted_affinities = sorted(valid_affinities)
        best = sorted_affinities[0]
        second_best = sorted_affinities[1] if len(sorted_affinities) > 1 else 0

        # Selectivity score (larger = more selective)
        selectivity_score = second_best - best

        # Determine if selective
        if selectivity_score >= 2.0:  # 2 kcal/mol difference
            selective_for = row['best_target']
            classification = "Selective"
        else:
            selective_for = "Non-selective"
            classification = "Promiscuous"

        selectivity_data.append({
            'ligand_id': row['ligand_id'],
            'best_target': row['best_target'],
            'best_affinity': best,
            'second_best_affinity': second_best,
            'selectivity_score': selectivity_score,
            'classification': classification,
            'selective_for': selective_for
        })

    return pd.DataFrame(selectivity_data)


def create_multi_target_summary(
    target_results: Dict,
    comparative_df: pd.DataFrame,
    selectivity_df: pd.DataFrame
) -> Dict:
    """Create summary statistics for multi-target screening"""

    summary = {
        'timestamp': datetime.now().isoformat(),
        'num_targets': len(target_results),
        'successful_targets': sum(1 for r in target_results.values() if r.get('success', False)),
        'failed_targets': sum(1 for r in target_results.values() if not r.get('success', False)),
        'total_molecules': len(comparative_df) if not comparative_df.empty else 0,
        'targets': {}
    }

    # Per-target stats
    for pdb_id, data in target_results.items():
        if data.get('success', False):
            summary['targets'][pdb_id] = {
                'num_hits': data['num_hits'],
                'best_affinity': float(data['best_affinity']),
                'hit_rate': (data['num_hits'] / data['num_molecules'] * 100) if data['num_molecules'] > 0 else 0
            }

    # Selectivity stats
    if not selectivity_df.empty:
        summary['selectivity'] = {
            'selective_molecules': int((selectivity_df['classification'] == 'Selective').sum()),
            'promiscuous_molecules': int((selectivity_df['classification'] == 'Promiscuous').sum()),
            'mean_selectivity_score': float(selectivity_df['selectivity_score'].mean())
        }

    # Multi-target hits
    if not comparative_df.empty and 'num_targets_hit' in comparative_df.columns:
        summary['multi_target_hits'] = {
            '1_target': int((comparative_df['num_targets_hit'] == 1).sum()),
            '2_targets': int((comparative_df['num_targets_hit'] == 2).sum()),
            '3+_targets': int((comparative_df['num_targets_hit'] >= 3).sum())
        }

    return summary


def print_multi_target_summary(summary: Dict, target_proteins: List[str]):
    """Print formatted summary to console"""

    print("\n" + "="*70)
    print("MULTI-TARGET SCREENING SUMMARY")
    print("="*70)

    print(f"\nTargets Screened: {summary['num_targets']}")
    print(f"Successful: {summary['successful_targets']}")
    print(f"Total Molecules: {summary['total_molecules']}")

    print("\nPer-Target Results:")
    print("-" * 70)

    for pdb_id in target_proteins:
        if pdb_id in summary['targets']:
            data = summary['targets'][pdb_id]
            print(f"  {pdb_id}:")
            print(f"    Hits: {data['num_hits']}")
            print(f"    Best Affinity: {data['best_affinity']:.2f} kcal/mol")
            print(f"    Hit Rate: {data['hit_rate']:.1f}%")

    if 'selectivity' in summary:
        print("\nSelectivity Analysis:")
        print("-" * 70)
        print(f"  Selective molecules: {summary['selectivity']['selective_molecules']}")
        print(f"  Promiscuous molecules: {summary['selectivity']['promiscuous_molecules']}")
        print(f"  Mean selectivity score: {summary['selectivity']['mean_selectivity_score']:.2f} kcal/mol")

    if 'multi_target_hits' in summary:
        print("\nMulti-Target Binding:")
        print("-" * 70)
        print(f"  Bind to 1 target: {summary['multi_target_hits']['1_target']}")
        print(f"  Bind to 2 targets: {summary['multi_target_hits']['2_targets']}")
        print(f"  Bind to 3+ targets: {summary['multi_target_hits']['3+_targets']}")

    print("\n" + "="*70)
