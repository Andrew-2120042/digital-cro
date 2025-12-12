"""
Test pharmacophore analysis
"""

from pathlib import Path
from scripts.utils.pharmacophore import analyze_pharmacophore_from_results
from scripts.utils.pharmacophore_viz import (
    plot_feature_distribution,
    plot_feature_importance,
    plot_pharmacophore_summary
)

def test_pharmacophore():
    """Test pharmacophore analysis on existing results"""

    # Use complete_workflow results
    results_csv = "data/outputs/complete_workflow/final_results.csv"
    output_dir = "data/outputs/complete_workflow/pharmacophore_test"

    print("🔬 Testing Pharmacophore Analysis")
    print("=" * 60)

    # Run analysis
    print("\n1. Analyzing features...")
    results = analyze_pharmacophore_from_results(
        results_csv=results_csv,
        top_n=10,  # Analyze top 10 hits
        min_occurrence=0.7,  # 70% occurrence threshold
        output_dir=output_dir
    )

    print(f"\n✓ Analyzed {results['hypothesis']['n_molecules']} molecules")

    # Display hypothesis
    print("\n2. Pharmacophore Hypothesis:")
    print("-" * 60)
    print(results['hypothesis']['hypothesis_text'])
    print("\nCommon features:")
    for feat, count in results['hypothesis']['common_features'].items():
        print(f"  - {feat}: {count}")

    # Feature importance
    print("\n3. Feature Importance:")
    print("-" * 60)
    print(results['importance'].to_string(index=False))

    # Generate visualizations
    print("\n4. Generating visualizations...")

    plot_feature_distribution(
        results['feature_df'],
        f"{output_dir}/feature_distribution.png"
    )
    print("  ✓ Feature distribution plot")

    plot_feature_importance(
        results['importance'],
        f"{output_dir}/feature_importance.png"
    )
    print("  ✓ Feature importance plot")

    plot_pharmacophore_summary(
        results['feature_df'],
        results['hypothesis'],
        f"{output_dir}/pharmacophore_summary.png"
    )
    print("  ✓ Summary plot")

    print(f"\n✅ Pharmacophore analysis complete!")
    print(f"📁 Results saved to: {output_dir}")

    # Verify files exist
    print("\n5. Verifying output files...")
    expected_files = [
        'pharmacophore_features.csv',
        'feature_importance.csv',
        'pharmacophore_hypothesis.txt',
        'feature_distribution.png',
        'feature_importance.png',
        'pharmacophore_summary.png'
    ]

    output_path = Path(output_dir)
    all_exist = True
    for filename in expected_files:
        file_path = output_path / filename
        if file_path.exists():
            print(f"  ✓ {filename}")
        else:
            print(f"  ✗ {filename} MISSING")
            all_exist = False

    if all_exist:
        print("\n✅ All files created successfully!")
    else:
        print("\n❌ Some files are missing!")

    return results

if __name__ == "__main__":
    test_pharmacophore()
