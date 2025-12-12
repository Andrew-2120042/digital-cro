"""
Pharmacophore Visualization

Create plots showing feature distributions and correlations.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
import numpy as np

sns.set_style("whitegrid")


def plot_feature_distribution(
    feature_df: pd.DataFrame,
    output_path: str,
    dpi: int = 300
):
    """
    Plot distribution of pharmacophore features.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Pharmacophore Feature Distribution', fontsize=16, fontweight='bold')

    features = [
        ('num_donors', 'H-Bond Donors'),
        ('num_acceptors', 'H-Bond Acceptors'),
        ('num_aromatic', 'Aromatic Rings'),
        ('num_hydrophobic', 'Hydrophobic Regions'),
        ('num_positive', 'Positive Charges'),
        ('num_negative', 'Negative Charges')
    ]

    for idx, (col, label) in enumerate(features):
        ax = axes[idx // 3, idx % 3]

        if col in feature_df.columns:
            data = feature_df[col]

            # Histogram
            max_val = int(data.max()) + 2
            bins = range(max_val) if max_val > 0 else [0, 1]

            ax.hist(data, bins=bins,
                   alpha=0.7, color='steelblue', edgecolor='black')

            ax.set_xlabel(label, fontsize=12)
            ax.set_ylabel('Count', fontsize=12)
            ax.set_title(f'{label} (Mean: {data.mean():.1f})', fontsize=11)
            ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_feature_importance(
    importance_df: pd.DataFrame,
    output_path: str,
    dpi: int = 300
):
    """
    Plot feature importance (correlation with binding).
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Bar plot
    colors = ['green' if x < 0 else 'red' for x in importance_df['correlation']]

    ax.barh(importance_df['feature_type'], importance_df['correlation'], color=colors, alpha=0.7)

    ax.set_xlabel('Correlation with Binding Affinity', fontsize=12, fontweight='bold')
    ax.set_ylabel('Feature Type', fontsize=12, fontweight='bold')
    ax.set_title('Feature Importance\n(Negative correlation = better binding)',
                fontsize=14, fontweight='bold')

    ax.axvline(x=0, color='black', linestyle='--', linewidth=1)
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_pharmacophore_summary(
    feature_df: pd.DataFrame,
    hypothesis: dict,
    output_path: str,
    dpi: int = 300
):
    """
    Create summary visualization with hypothesis.
    """
    fig = plt.figure(figsize=(12, 8))

    # Title
    fig.suptitle('Pharmacophore Hypothesis Summary', fontsize=16, fontweight='bold')

    # Hypothesis text box
    ax_text = plt.subplot(3, 1, 1)
    ax_text.axis('off')

    hypothesis_text = hypothesis['hypothesis_text']

    # Wrap text for better display
    import textwrap
    wrapped_text = textwrap.fill(hypothesis_text, width=80)

    ax_text.text(0.5, 0.5, wrapped_text,
                ha='center', va='center',
                fontsize=13,
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    # Feature occurrence
    ax_occur = plt.subplot(3, 2, 3)

    feature_cols = ['num_donors', 'num_acceptors', 'num_aromatic', 'num_hydrophobic']
    available_cols = [col for col in feature_cols if col in feature_df.columns]
    means = [feature_df[col].mean() for col in available_cols]
    labels = [col.replace('num_', '').title() for col in available_cols]

    if means:
        ax_occur.bar(labels, means, color='steelblue', alpha=0.7)
        ax_occur.set_ylabel('Average Count', fontsize=11)
        ax_occur.set_title('Average Feature Counts', fontsize=12, fontweight='bold')
        ax_occur.tick_params(axis='x', rotation=45)

    # Feature presence (% of molecules)
    ax_presence = plt.subplot(3, 2, 4)

    presence = [(feature_df[col] > 0).sum() / len(feature_df) * 100
                for col in available_cols]

    if presence:
        ax_presence.bar(labels, presence, color='coral', alpha=0.7)
        ax_presence.set_ylabel('% of Molecules', fontsize=11)
        ax_presence.set_title('Feature Presence', fontsize=12, fontweight='bold')
        ax_presence.tick_params(axis='x', rotation=45)
        ax_presence.set_ylim([0, 100])

    # Feature distribution heatmap
    ax_heat = plt.subplot(3, 1, 3)

    # Create matrix: molecules × features
    if available_cols:
        matrix_data = feature_df[available_cols].values

        im = ax_heat.imshow(matrix_data.T, aspect='auto', cmap='YlOrRd', interpolation='nearest')

        ax_heat.set_yticks(range(len(labels)))
        ax_heat.set_yticklabels(labels)
        ax_heat.set_xlabel('Molecule Index', fontsize=11)
        ax_heat.set_title('Feature Distribution Across Molecules', fontsize=12, fontweight='bold')

        plt.colorbar(im, ax=ax_heat, label='Feature Count')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
