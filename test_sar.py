"""
Test SAR Analysis
"""

from pathlib import Path
from scripts.utils.sar_analysis import analyze_sar_from_results
from scripts.utils.sar_viz import (
    plot_scaffold_comparison,
    plot_mmp_network,
    plot_property_correlations,
    plot_sar_summary
)

def test_sar():
    """Test SAR analysis on existing results"""

    # Use complete_workflow results
    results_csv = "data/outputs/complete_workflow/final_results.csv"
    output_dir = "data/outputs/complete_workflow/sar_test"

    print("🧬 Testing SAR Analysis")
    print("=" * 60)

    # Run analysis
    print("\n1. Analyzing structure-activity relationships...")
    results = analyze_sar_from_results(
        results_csv=results_csv,
        min_molecules_per_scaffold=2,
        output_dir=output_dir
    )

    print(f"\n✓ Analysis complete!")
    print(f"  - Scaffolds identified: {results['n_scaffolds']}")
    print(f"  - Matched pairs found: {len(results['mmps'])}")

    # Display scaffolds
    print("\n2. Top Scaffolds:")
    print("-" * 60)
    for i, scaffold in enumerate(results['scaffolds'][:5], 1):
        print(f"  {i}. {scaffold['n_molecules']} molecules, best: {scaffold['best_affinity']:.2f} kcal/mol")

    # Display MMPs
    if results['mmps']:
        print("\n3. Top Matched Molecular Pairs:")
        print("-" * 60)
        for i, mmp in enumerate(results['mmps'][:5], 1):
            print(f"  {i}. {mmp['mol1_id']} → {mmp['mol2_id']}")
            print(f"     Δ Affinity: {mmp['delta_affinity']:.2f} kcal/mol ({mmp['improvement']})")

    # Display correlations
    if not results['correlations'].empty:
        print("\n4. Property Correlations:")
        print("-" * 60)
        print(results['correlations'].to_string(index=False))

    # Display recommendations
    print("\n5. Recommendations:")
    print("-" * 60)
    for i, rec in enumerate(results['recommendations'], 1):
        print(f"  {i}. {rec}")

    # Generate visualizations
    print("\n6. Generating visualizations...")

    if results['scaffolds']:
        plot_scaffold_comparison(
            results['scaffolds'],
            f"{output_dir}/scaffold_comparison.png"
        )
        print("  ✓ Scaffold comparison plot")

    if results['mmps']:
        plot_mmp_network(
            results['mmps'],
            f"{output_dir}/matched_pairs.png"
        )
        print("  ✓ Matched pairs plot")

    if not results['correlations'].empty:
        plot_property_correlations(
            results['correlations'],
            f"{output_dir}/property_correlations.png"
        )
        print("  ✓ Property correlations plot")

    # Summary plot
    plot_sar_summary(
        results['scaffolds'],
        results['mmps'],
        results['correlations'],
        f"{output_dir}/sar_summary.png"
    )
    print("  ✓ Summary plot")

    print(f"\n✅ SAR analysis complete!")
    print(f"📁 Results saved to: {output_dir}")

    # Verify files exist
    print("\n7. Verifying output files...")
    expected_files = [
        'scaffold_summary.csv',
        'matched_pairs.csv',
        'property_correlations.csv',
        'sar_recommendations.txt',
        'scaffold_comparison.png',
        'matched_pairs.png',
        'property_correlations.png',
        'sar_summary.png'
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
    test_sar()
