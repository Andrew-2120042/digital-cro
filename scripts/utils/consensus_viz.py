"""
Consensus Docking Visualization

Create plots showing method agreement and confidence.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict

sns.set_style("whitegrid")


def plot_method_agreement(
    df: pd.DataFrame,
    methods: List[str],
    output_path: str,
    top_n: int = 20,
    dpi: int = 300
):
    """
    Scatter plot comparing docking methods.
    """
    if len(methods) < 2:
        return

    fig, axes = plt.subplots(1, len(methods)-1, figsize=(12, 5))

    if len(methods) == 2:
        axes = [axes]

    fig.suptitle('Docking Method Agreement', fontsize=14, fontweight='bold')

    # Compare each method with Vina
    base_method = methods[0]

    for idx, compare_method in enumerate(methods[1:]):
        ax = axes[idx]

        # Get top N molecules
        df_top = df.nsmallest(top_n, 'consensus_score')

        x = df_top[f'{base_method}_score']
        y = df_top[f'{compare_method}_score']

        # Scatter plot
        ax.scatter(x, y, alpha=0.6, s=50, edgecolors='black')

        # Diagonal line (perfect agreement)
        min_val = min(x.min(), y.min())
        max_val = max(x.max(), y.max())
        ax.plot([min_val, max_val], [min_val, max_val],
               'r--', alpha=0.5, label='Perfect Agreement')

        # Correlation
        corr = np.corrcoef(x, y)[0, 1]

        ax.set_xlabel(f'{base_method.capitalize()} Score (kcal/mol)', fontsize=10)
        ax.set_ylabel(f'{compare_method.capitalize()} Score (kcal/mol)', fontsize=10)
        ax.set_title(f'{base_method.capitalize()} vs {compare_method.capitalize()}\n(r={corr:.2f})',
                    fontsize=11)
        ax.legend()
        ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_confidence_distribution(
    df: pd.DataFrame,
    output_path: str,
    dpi: int = 300
):
    """
    Visualize confidence levels across molecules.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Consensus Confidence Analysis', fontsize=14, fontweight='bold')

    # 1. Confidence pie chart
    ax = axes[0, 0]

    confidence_counts = df['confidence'].value_counts()
    colors_conf = {'high': 'green', 'medium': 'orange', 'low': 'red'}
    colors_pie = [colors_conf.get(c, 'gray') for c in confidence_counts.index]

    ax.pie(confidence_counts.values, labels=confidence_counts.index,
          autopct='%1.1f%%', colors=colors_pie, startangle=90)
    ax.set_title('Confidence Distribution', fontsize=12, fontweight='bold')

    # 2. Score standard deviation histogram
    ax = axes[0, 1]

    ax.hist(df['score_std'], bins=20, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(2.0, color='red', linestyle='--', linewidth=2,
              label='Agreement Threshold (2.0)')

    ax.set_xlabel('Score Standard Deviation (kcal/mol)', fontsize=10)
    ax.set_ylabel('Count', fontsize=10)
    ax.set_title('Method Agreement', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    # 3. Consensus vs individual scores
    ax = axes[1, 0]

    df_top = df.nsmallest(20, 'consensus_score')

    x = range(len(df_top))
    ax.plot(x, df_top['consensus_score'], 'o-', label='Consensus',
           color='black', linewidth=2, markersize=8)

    # Plot individual methods
    method_cols = [col for col in df.columns if col.endswith('_score') and col != 'consensus_score']
    colors_methods = ['blue', 'orange', 'green']

    for idx, col in enumerate(method_cols):
        method_name = col.replace('_score', '').capitalize()
        ax.plot(x, df_top[col], 'o--', label=method_name,
               alpha=0.6, color=colors_methods[idx % len(colors_methods)])

    ax.set_xlabel('Molecule Rank', fontsize=10)
    ax.set_ylabel('Binding Affinity (kcal/mol)', fontsize=10)
    ax.set_title('Top 20 Molecules: Method Comparison', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # 4. Score range (disagreement)
    ax = axes[1, 1]

    df_sorted = df.sort_values('score_range', ascending=False).head(15)

    y_pos = range(len(df_sorted))
    colors_bar = ['red' if r > 3.0 else 'orange' if r > 2.0 else 'green'
                  for r in df_sorted['score_range']]

    ax.barh(y_pos, df_sorted['score_range'], color=colors_bar, alpha=0.7, edgecolor='black')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_sorted['ligand_id'], fontsize=8)
    ax.set_xlabel('Score Range (kcal/mol)', fontsize=10)
    ax.set_title('Method Disagreement (Top 15)', fontsize=12, fontweight='bold')
    ax.axvline(2.0, color='black', linestyle='--', alpha=0.5)
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_consensus_summary(
    df: pd.DataFrame,
    analysis: Dict,
    output_path: str,
    dpi: int = 300
):
    """
    Summary visualization of consensus results.
    """
    fig = plt.figure(figsize=(14, 8))

    # Title
    fig.suptitle('Consensus Docking Summary', fontsize=16, fontweight='bold')

    # Stats box
    ax1 = plt.subplot(2, 2, 1)
    ax1.axis('off')

    stats_text = f"""
    CONSENSUS STATISTICS

    Total Molecules: {analysis['total_molecules']}
    High Confidence: {analysis['high_confidence']} ({analysis['high_confidence_rate']*100:.1f}%)
    Low Confidence: {analysis['low_confidence']}
    Discrepant: {analysis['discrepant_molecules']}

    Mean Consensus Score: {analysis['mean_consensus_score']:.2f} kcal/mol
    Mean Agreement (Std): {analysis['mean_score_std']:.2f} kcal/mol
    """

    ax1.text(0.5, 0.5, stats_text,
            ha='center', va='center',
            fontsize=11, family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    # Top hits table
    ax2 = plt.subplot(2, 2, 2)
    ax2.axis('off')

    df_top = df.nsmallest(10, 'consensus_score')

    table_data = []
    for _, row in df_top.iterrows():
        table_data.append([
            row['ligand_id'][:15],
            f"{row['consensus_score']:.2f}",
            row['confidence']
        ])

    table = ax2.table(cellText=table_data,
                     colLabels=['Molecule', 'Score', 'Confidence'],
                     cellLoc='center',
                     loc='center',
                     colWidths=[0.5, 0.25, 0.25])

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    ax2.set_title('Top 10 Consensus Hits', fontsize=12, fontweight='bold', pad=20)

    # Confidence bar chart
    ax3 = plt.subplot(2, 2, 3)

    conf_counts = df['confidence'].value_counts()
    colors_conf = {'high': 'green', 'medium': 'orange', 'low': 'red'}
    colors_bar = [colors_conf.get(c, 'gray') for c in conf_counts.index]

    ax3.bar(conf_counts.index, conf_counts.values, color=colors_bar, alpha=0.7, edgecolor='black')
    ax3.set_ylabel('Number of Molecules', fontsize=10)
    ax3.set_title('Confidence Levels', fontsize=12, fontweight='bold')
    ax3.grid(axis='y', alpha=0.3)

    # Score distribution
    ax4 = plt.subplot(2, 2, 4)

    ax4.hist(df['consensus_score'], bins=20, color='steelblue', alpha=0.7, edgecolor='black')
    ax4.axvline(df['consensus_score'].median(), color='red', linestyle='--',
               linewidth=2, label=f'Median: {df["consensus_score"].median():.2f}')

    ax4.set_xlabel('Consensus Score (kcal/mol)', fontsize=10)
    ax4.set_ylabel('Count', fontsize=10)
    ax4.set_title('Consensus Score Distribution', fontsize=12, fontweight='bold')
    ax4.legend()
    ax4.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
