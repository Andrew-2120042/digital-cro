"""
Docking result analysis module for Digital CRO.

Functions for:
    - Ranking hits by binding affinity
    - Filtering by score thresholds
    - Diversity selection
    - Chemical clustering
    - Report generation
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import json
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw, Descriptors
from rdkit import RDLogger
import logging

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def rank_docking_results(results_csv: str, top_n: int = 50) -> pd.DataFrame:
    """
    Load and rank docking results.

    Args:
        results_csv: CSV from batch docking
        top_n: Return top N hits

    Returns:
        DataFrame with top hits, sorted by affinity

    Adds columns:
        - rank (1, 2, 3, ...)
        - binding_efficiency (affinity / molecular_weight)
        - pass_threshold (True if affinity < -7.0)
    """
    df = pd.read_csv(results_csv)

    # Filter successful only
    df = df[df['success'] == True].copy()

    if df.empty:
        logger.warning("No successful docking results found")
        return df

    # Sort by binding affinity (lowest = best)
    df = df.sort_values('binding_affinity').reset_index(drop=True)

    # Add rank
    df['rank'] = range(1, len(df) + 1)

    # Add pass threshold
    df['pass_threshold'] = df['binding_affinity'] < -7.0

    # Calculate molecular properties
    if 'smiles' in df.columns:
        df['mol_weight'] = df['smiles'].apply(lambda s: _get_mol_weight(s))
        df['binding_efficiency'] = -df['binding_affinity'] / df['mol_weight'] * 1000

    # Return top N
    return df.head(top_n)


def _get_mol_weight(smiles: str) -> float:
    """Helper: get molecular weight from SMILES"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Descriptors.MolWt(mol)
    except:
        pass
    return 0.0


def filter_hits_by_score(
    df: pd.DataFrame,
    threshold: float = -7.0,
    binding_efficiency_min: float = None
) -> pd.DataFrame:
    """
    Filter hits by binding affinity threshold.

    Typical thresholds:
        -6.0: Weak binders
        -7.0: Moderate (common threshold)
        -8.0: Strong binders
        -9.0: Very strong (rare)

    Args:
        df: DataFrame with docking results
        threshold: Binding affinity cutoff (kcal/mol)
        binding_efficiency_min: Optional min binding efficiency

    Returns:
        Filtered DataFrame
    """
    df_filtered = df[df['binding_affinity'] < threshold].copy()

    if binding_efficiency_min and 'binding_efficiency' in df_filtered.columns:
        df_filtered = df_filtered[df_filtered['binding_efficiency'] >= binding_efficiency_min]

    logger.info(f"Filtered: {len(df_filtered)}/{len(df)} hits pass threshold {threshold} kcal/mol")

    return df_filtered


def select_diverse_hits(
    df: pd.DataFrame,
    n_hits: int = 20,
    diversity_threshold: float = 0.5
) -> pd.DataFrame:
    """
    Select diverse subset of top hits.

    Balance between:
        - High binding affinity (good scores)
        - Chemical diversity (different scaffolds)

    Algorithm:
        1. Sort by affinity
        2. Take top hit
        3. For remaining, add if Tanimoto < diversity_threshold
        4. Continue until n_hits reached

    Args:
        df: DataFrame with docking results (must have 'smiles' column)
        n_hits: Target number of diverse hits
        diversity_threshold: Tanimoto similarity threshold (0-1)
            0.5 = moderate diversity
            0.3 = high diversity

    Returns:
        DataFrame with diverse molecules
    """
    if 'smiles' not in df.columns:
        logger.error("DataFrame must have 'smiles' column for diversity selection")
        return df.head(n_hits)

    smiles_list = df['smiles'].tolist()

    # Convert to fingerprints
    fps = []
    valid_indices = []

    for idx, smiles in enumerate(smiles_list):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
                fps.append(fp)
                valid_indices.append(idx)
        except:
            continue

    if not fps:
        logger.warning("No valid molecules for diversity selection")
        return df.head(n_hits)

    # Select diverse subset
    selected_indices = [0]  # Start with best (index 0 in valid_indices)

    for i in range(1, len(fps)):
        if len(selected_indices) >= n_hits:
            break

        # Calculate similarity to all selected
        current_fp = fps[i]
        max_similarity = 0.0

        for selected_idx in selected_indices:
            similarity = DataStructs.TanimotoSimilarity(current_fp, fps[selected_idx])
            max_similarity = max(max_similarity, similarity)

        # Add if sufficiently different
        if max_similarity < diversity_threshold:
            selected_indices.append(i)

    # Map back to original DataFrame indices
    df_indices = [valid_indices[i] for i in selected_indices]

    logger.info(f"Selected {len(df_indices)} diverse hits from {len(df)} total")

    return df.iloc[df_indices].copy()


def cluster_hits_by_similarity(
    smiles_list: List[str],
    max_clusters: int = 10
) -> Dict[int, List[int]]:
    """
    Cluster hits by chemical similarity.

    For understanding chemical diversity in hit set.

    Args:
        smiles_list: List of SMILES from hits
        max_clusters: Target number of clusters

    Returns:
        dict: {cluster_id: [list of indices]}

    Uses Butina clustering on Morgan fingerprints.
    """
    from rdkit.ML.Cluster import Butina

    # Generate fingerprints
    fps = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
                fps.append(fp)
            else:
                fps.append(None)
        except:
            fps.append(None)

    # Calculate distance matrix
    nfps = len(fps)
    dists = []
    for i in range(nfps):
        for j in range(i):
            if fps[i] and fps[j]:
                dist = 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
            else:
                dist = 1.0
            dists.append(dist)

    # Butina clustering
    clusters = Butina.ClusterData(dists, nfps, 0.5, isDistData=True)

    # Convert to dict
    cluster_dict = {}
    for cluster_id, cluster in enumerate(clusters):
        cluster_dict[cluster_id] = list(cluster)

    logger.info(f"Formed {len(cluster_dict)} clusters from {len(smiles_list)} molecules")

    return cluster_dict


def generate_hit_report(
    df: pd.DataFrame,
    output_dir: str,
    top_n: int = 50
):
    """
    Generate comprehensive hit analysis report.

    Creates:
        - hits_summary.csv (top N with all data)
        - binding_distribution.png (histogram)
        - top_structures.png (grid of 2D structures)
        - diversity_analysis.json (clustering info)

    Args:
        df: DataFrame with docking results
        output_dir: Where to save report files
        top_n: Number of top hits to include
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Save summary CSV
    df_top = df.head(top_n)
    summary_csv = output_dir / 'hits_summary.csv'
    df_top.to_csv(summary_csv, index=False)
    logger.info(f"Saved summary: {summary_csv}")

    # 2. Binding affinity distribution
    plt.figure(figsize=(10, 6))
    plt.hist(df['binding_affinity'], bins=30, edgecolor='black')
    plt.axvline(-7.0, color='red', linestyle='--', label='Threshold (-7.0)')
    plt.xlabel('Binding Affinity (kcal/mol)')
    plt.ylabel('Count')
    plt.title(f'Binding Affinity Distribution ({len(df)} molecules)')
    plt.legend()
    plt.tight_layout()

    dist_png = output_dir / 'binding_distribution.png'
    plt.savefig(dist_png, dpi=150)
    plt.close()
    logger.info(f"Saved distribution: {dist_png}")

    # 3. Top structures (if SMILES available)
    if 'smiles' in df.columns:
        molecules = []
        legends = []

        for idx, row in df.head(20).iterrows():
            try:
                mol = Chem.MolFromSmiles(row['smiles'])
                if mol:
                    molecules.append(mol)
                    legends.append(f"{row['ligand_id']}\n{row['binding_affinity']:.2f} kcal/mol")
            except:
                continue

        if molecules:
            img = Draw.MolsToGridImage(
                molecules,
                molsPerRow=5,
                subImgSize=(200, 200),
                legends=legends
            )

            struct_png = output_dir / 'top_structures.png'
            img.save(struct_png)
            logger.info(f"Saved structures: {struct_png}")

    # 4. Diversity analysis
    if 'smiles' in df.columns and len(df) >= 10:
        smiles_list = df.head(100)['smiles'].tolist()
        clusters = cluster_hits_by_similarity(smiles_list)

        diversity_json = output_dir / 'diversity_analysis.json'
        with open(diversity_json, 'w') as f:
            json.dump({
                'num_molecules': len(smiles_list),
                'num_clusters': len(clusters),
                'cluster_sizes': {k: len(v) for k, v in clusters.items()},
                'avg_cluster_size': np.mean([len(v) for v in clusters.values()])
            }, f, indent=2)

        logger.info(f"Saved diversity: {diversity_json}")

    logger.info("Hit report generation complete")
