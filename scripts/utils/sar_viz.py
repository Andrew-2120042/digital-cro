"""
SAR Visualization

Create plots for structure-activity relationship analysis.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict

sns.set_style("whitegrid")


def plot_scaffold_comparison(
    scaffold_analyses: List[Dict],
    output_path: str,
    dpi: int = 300
):
    """
    Bar chart comparing scaffolds by best affinity.
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    # Prepare data (top 10 scaffolds)
    top_scaffolds = scaffold_analyses[:10]

    names = [f"Scaffold {i+1}\n(n={s['n_molecules']})"
             for i, s in enumerate(top_scaffolds)]
    best_affinities = [s['best_affinity'] for s in top_scaffolds]
    mean_affinities = [s['mean_affinity'] for s in top_scaffolds]

    x = np.arange(len(names))
    width = 0.35

    # Best affinity bars
    bars1 = ax.bar(x - width/2, best_affinities, width,
                   label='Best Binder', color='green', alpha=0.7)

    # Mean affinity bars
    bars2 = ax.bar(x + width/2, mean_affinities, width,
                   label='Mean Affinity', color='steelblue', alpha=0.7)

    ax.set_xlabel('Scaffold', fontsize=12, fontweight='bold')
    ax.set_ylabel('Binding Affinity (kcal/mol)', fontsize=12, fontweight='bold')
    ax.set_title('Scaffold Comparison by Binding Affinity',
                fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.1f}',
                   ha='center', va='bottom', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_mmp_network(
    mmps: List[Dict],
    output_path: str,
    top_n: int = 15,
    dpi: int = 300
):
    """
    Visualize matched molecular pairs as network/scatter.
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    # Take top N pairs
    top_mmps = mmps[:top_n]

    # Prepare data
    improvements = []
    mol1_affinities = []
    mol2_affinities = []
    labels = []

    for mmp in top_mmps:
        mol1_affinities.append(mmp['mol1_affinity'])
        mol2_affinities.append(mmp['mol2_affinity'])
        improvements.append(mmp['delta_affinity'])
        labels.append(f"{mmp['mol1_id'][:8]}→{mmp['mol2_id'][:8]}")

    # Color by improvement
    colors = ['green' if d < 0 else 'red' for d in improvements]

    # Scatter plot
    scatter = ax.scatter(mol1_affinities, mol2_affinities,
                        c=colors, s=100, alpha=0.6, edgecolors='black')

    # Diagonal line (no change)
    min_val = min(mol1_affinities + mol2_affinities)
    max_val = max(mol1_affinities + mol2_affinities)
    ax.plot([min_val, max_val], [min_val, max_val],
           'k--', alpha=0.3, label='No Change')

    # Labels for interesting pairs
    for i in range(min(5, len(top_mmps))):
        ax.annotate(labels[i],
                   (mol1_affinities[i], mol2_affinities[i]),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, alpha=0.7)

    ax.set_xlabel('Molecule 1 Affinity (kcal/mol)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Molecule 2 Affinity (kcal/mol)', fontsize=12, fontweight='bold')
    ax.set_title('Matched Molecular Pairs\n(Green = Improved, Red = Worsened)',
                fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_property_correlations(
    correlations: pd.DataFrame,
    output_path: str,
    dpi: int = 300
):
    """
    Bar chart of property-activity correlations.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Sort by absolute correlation
    corr_sorted = correlations.sort_values('correlation', key=abs, ascending=False)

    colors = ['green' if c < 0 else 'red' for c in corr_sorted['correlation']]

    bars = ax.barh(corr_sorted['property'], corr_sorted['correlation'],
                   color=colors, alpha=0.7, edgecolor='black')

    # Mark significant correlations
    for i, (_, row) in enumerate(corr_sorted.iterrows()):
        if row['significant'] == 'Yes':
            ax.text(row['correlation'], i, ' *',
                   ha='left' if row['correlation'] > 0 else 'right',
                   va='center', fontsize=16, fontweight='bold')

    ax.set_xlabel('Correlation with Binding Affinity', fontsize=12, fontweight='bold')
    ax.set_ylabel('Molecular Property', fontsize=12, fontweight='bold')
    ax.set_title('Property-Activity Relationships\n(* = statistically significant, p<0.05)',
                fontsize=14, fontweight='bold')
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1)
    ax.grid(axis='x', alpha=0.3)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='green', alpha=0.7, label='Negative (Better binding)'),
        Patch(facecolor='red', alpha=0.7, label='Positive (Worse binding)')
    ]
    ax.legend(handles=legend_elements, loc='best')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_sar_summary(
    scaffold_analyses: List[Dict],
    mmps: List[Dict],
    correlations: pd.DataFrame,
    output_path: str,
    dpi: int = 300
):
    """
    Comprehensive SAR summary visualization.
    """
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

    # Title
    fig.suptitle('Structure-Activity Relationship (SAR) Analysis Summary',
                fontsize=16, fontweight='bold')

    # 1. Scaffold comparison (top left)
    ax1 = fig.add_subplot(gs[0, :])

    if scaffold_analyses:
        top_5 = scaffold_analyses[:5]
        names = [f"S{i+1}" for i in range(len(top_5))]
        best_aff = [s['best_affinity'] for s in top_5]

        ax1.bar(names, best_aff, color='steelblue', alpha=0.7)
        ax1.set_ylabel('Best Affinity (kcal/mol)', fontsize=10)
        ax1.set_title('Top 5 Scaffolds', fontsize=12, fontweight='bold')
        ax1.grid(axis='y', alpha=0.3)

    # 2. MMP improvements (middle left)
    ax2 = fig.add_subplot(gs[1, 0])

    if mmps:
        improvements = [mmp['delta_affinity'] for mmp in mmps[:10]]
        colors_mmp = ['green' if d < 0 else 'red' for d in improvements]

        ax2.barh(range(len(improvements)), improvements, color=colors_mmp, alpha=0.7)
        ax2.set_xlabel('Δ Affinity (kcal/mol)', fontsize=10)
        ax2.set_title('Top 10 Structural Changes', fontsize=12, fontweight='bold')
        ax2.axvline(x=0, color='black', linestyle='--')
        ax2.grid(axis='x', alpha=0.3)

    # 3. Property correlations (middle right)
    ax3 = fig.add_subplot(gs[1, 1])

    if not correlations.empty:
        colors_corr = ['green' if c < 0 else 'red' for c in correlations['correlation']]

        ax3.barh(correlations['property'], correlations['correlation'],
                color=colors_corr, alpha=0.7)
        ax3.set_xlabel('Correlation', fontsize=10)
        ax3.set_title('Property Correlations', fontsize=12, fontweight='bold')
        ax3.axvline(x=0, color='black', linestyle='--')
        ax3.grid(axis='x', alpha=0.3)

    # 4. Statistics box (bottom)
    ax4 = fig.add_subplot(gs[2, :])
    ax4.axis('off')

    n_scaffolds = len(scaffold_analyses)
    n_mmps = len(mmps)
    best_aff = scaffold_analyses[0]['best_affinity'] if scaffold_analyses else 0
    largest_cliff = abs(mmps[0]['delta_affinity']) if mmps else 0
    n_sig_corr = (correlations['significant'] == 'Yes').sum() if not correlations.empty else 0

    stats_text = f"""
    SAR Analysis Statistics:

    • Total Scaffolds Identified: {n_scaffolds}
    • Matched Molecular Pairs Found: {n_mmps}
    • Best Overall Affinity: {best_aff:.2f} kcal/mol
    • Largest Activity Cliff: {largest_cliff:.2f} kcal/mol
    • Significant Property Correlations: {n_sig_corr}
    """

    ax4.text(0.5, 0.5, stats_text,
            ha='center', va='center',
            fontsize=11, family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3))

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
