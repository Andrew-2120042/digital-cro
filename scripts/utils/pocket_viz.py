"""
Pocket Comparison Visualization

Create plots showing pocket similarities and properties.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict

sns.set_style("whitegrid")


def plot_similarity_matrix(
    comparison_matrix: Dict[str, float],
    target_names: List[str],
    output_path: str,
    dpi: int = 300
):
    """
    Create heatmap of pocket similarity matrix.
    """
    # Create matrix
    n = len(target_names)
    matrix = np.ones((n, n))

    for i, target1 in enumerate(target_names):
        for j, target2 in enumerate(target_names):
            if i == j:
                continue

            key = f"{target1}_vs_{target2}" if i < j else f"{target2}_vs_{target1}"

            if key in comparison_matrix:
                matrix[i, j] = comparison_matrix[key]
                matrix[j, i] = comparison_matrix[key]

    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))

    im = ax.imshow(matrix, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')

    # Labels
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(target_names, rotation=45, ha='right')
    ax.set_yticklabels(target_names)

    # Add values
    for i in range(n):
        for j in range(n):
            text = ax.text(j, i, f'{matrix[i, j]:.2f}',
                          ha='center', va='center',
                          color='black' if matrix[i, j] > 0.5 else 'white',
                          fontsize=10, fontweight='bold')

    ax.set_title('Protein Pocket Similarity Matrix',
                fontsize=14, fontweight='bold', pad=20)

    plt.colorbar(im, ax=ax, label='Similarity Score')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_pocket_properties(
    pockets: Dict[str, Dict],
    output_path: str,
    dpi: int = 300
):
    """
    Bar chart comparing pocket properties.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Pocket Property Comparison', fontsize=16, fontweight='bold')

    target_names = list(pockets.keys())

    # 1. Pocket size
    ax = axes[0, 0]
    sizes = [pockets[t]['num_residues'] for t in target_names]
    ax.bar(target_names, sizes, color='steelblue', alpha=0.7)
    ax.set_ylabel('Number of Residues', fontsize=11)
    ax.set_title('Pocket Size', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    # 2. Hydrophobicity
    ax = axes[0, 1]
    hydrophobic = [pockets[t]['properties']['hydrophobic_fraction'] for t in target_names]
    ax.bar(target_names, hydrophobic, color='orange', alpha=0.7)
    ax.set_ylabel('Fraction', fontsize=11)
    ax.set_title('Hydrophobic Character', fontsize=12, fontweight='bold')
    ax.set_ylim([0, 1])
    ax.grid(axis='y', alpha=0.3)

    # 3. Charge distribution
    ax = axes[1, 0]
    positive = [pockets[t]['properties']['positive_fraction'] for t in target_names]
    negative = [pockets[t]['properties']['negative_fraction'] for t in target_names]

    x = np.arange(len(target_names))
    width = 0.35

    ax.bar(x - width/2, positive, width, label='Positive', color='blue', alpha=0.7)
    ax.bar(x + width/2, negative, width, label='Negative', color='red', alpha=0.7)

    ax.set_ylabel('Fraction', fontsize=11)
    ax.set_title('Charge Distribution', fontsize=12, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(target_names)
    ax.legend()
    ax.set_ylim([0, 0.5])
    ax.grid(axis='y', alpha=0.3)

    # 4. Polarity
    ax = axes[1, 1]
    polar = [pockets[t]['properties']['polar_fraction'] for t in target_names]
    ax.bar(target_names, polar, color='green', alpha=0.7)
    ax.set_ylabel('Fraction', fontsize=11)
    ax.set_title('Polar Character', fontsize=12, fontweight='bold')
    ax.set_ylim([0, 1])
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_comparison_summary(
    comparisons: List[Dict],
    selectivity: Dict,
    output_path: str,
    dpi: int = 300
):
    """
    Summary visualization of pocket comparisons.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Pocket Comparison Summary', fontsize=16, fontweight='bold')

    # 1. Similarity distribution
    ax = axes[0, 0]
    similarities = [c['overall_similarity'] for c in comparisons]

    ax.hist(similarities, bins=10, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(selectivity['avg_similarity'], color='red', linestyle='--',
              linewidth=2, label=f'Average: {selectivity["avg_similarity"]:.2f}')

    ax.set_xlabel('Similarity Score', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Pocket Similarity Distribution', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    # 2. Similarity components
    ax = axes[0, 1]

    components = ['overall', 'residue', 'property', 'size']
    avg_scores = [
        np.mean([c['overall_similarity'] for c in comparisons]),
        np.mean([c['residue_similarity'] for c in comparisons]),
        np.mean([c['property_similarity'] for c in comparisons]),
        np.mean([c['size_similarity'] for c in comparisons])
    ]

    colors_comp = ['steelblue', 'orange', 'green', 'purple']

    ax.bar(components, avg_scores, color=colors_comp, alpha=0.7, edgecolor='black')
    ax.set_ylabel('Average Similarity', fontsize=11)
    ax.set_title('Similarity Components', fontsize=12, fontweight='bold')
    ax.set_ylim([0, 1])
    ax.grid(axis='y', alpha=0.3)

    # 3. Top similar pairs
    ax = axes[1, 0]

    top_pairs = sorted(comparisons, key=lambda x: x['overall_similarity'], reverse=True)[:5]

    pair_labels = [f"{c['target1'][:4]}-{c['target2'][:4]}" for c in top_pairs]
    pair_scores = [c['overall_similarity'] for c in top_pairs]

    ax.barh(pair_labels, pair_scores, color='green', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Similarity Score', fontsize=11)
    ax.set_title('Top 5 Most Similar Pairs', fontsize=12, fontweight='bold')
    ax.set_xlim([0, 1])
    ax.grid(axis='x', alpha=0.3)

    # 4. Selectivity indicator
    ax = axes[1, 1]
    ax.axis('off')

    # Selectivity gauge
    selectivity_score = 1 - selectivity['avg_similarity']  # Invert (high similarity = low selectivity)

    if selectivity_score >= 0.7:
        color = 'green'
        level = 'HIGH'
    elif selectivity_score >= 0.4:
        color = 'orange'
        level = 'MODERATE'
    else:
        color = 'red'
        level = 'LOW'

    text = f"""
    SELECTIVITY ASSESSMENT

    Level: {level}
    Score: {selectivity_score:.2f}

    {selectivity['selectivity_level']}

    Comparisons: {selectivity['num_comparisons']}
    """

    ax.text(0.5, 0.5, text,
           ha='center', va='center',
           fontsize=12, family='monospace',
           bbox=dict(boxstyle='round', facecolor=color, alpha=0.3, pad=1))

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
