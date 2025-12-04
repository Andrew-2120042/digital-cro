"""
Test Pocket Comparison Analysis
"""

from pathlib import Path
from scripts.utils.pocket_comparison import compare_multi_target_pockets
from scripts.utils.pocket_viz import (
    plot_similarity_matrix,
    plot_pocket_properties,
    plot_comparison_summary
)


def test_pocket_comparison():
    """Test pocket comparison on existing multi-target results"""

    # Use multi_target_demo results
    multi_target_dir = "data/outputs/multi_target_demo"
    output_dir = "data/outputs/multi_target_demo/pocket_test"

    print("🧬 Testing Pocket Comparison Analysis")
    print("=" * 60)

    # Run analysis
    print("\n1. Analyzing protein pockets...")
    results = compare_multi_target_pockets(
        multi_target_dir=multi_target_dir,
        output_dir=output_dir
    )

    print(f"\n✓ Analysis complete!")
    print(f"  - Targets analyzed: {len(results['target_names'])}")
    print(f"  - Pairwise comparisons: {len(results['comparisons'])}")

    # Display target names
    print(f"\n2. Targets:")
    print("-" * 60)
    for target in results['target_names']:
        pocket = results['pockets'][target]
        print(f"  {target}: {pocket['num_residues']} residues")

    # Display comparisons
    if results['comparisons']:
        print(f"\n3. Pocket Comparisons:")
        print("-" * 60)
        for comp in results['comparisons']:
            print(f"  {comp['target1']} vs {comp['target2']}")
            print(f"     Overall similarity: {comp['overall_similarity']:.2f}")
            print(f"     {comp['interpretation']}")
            print()

    # Display selectivity
    print(f"4. Selectivity Assessment:")
    print("-" * 60)
    selectivity = results['selectivity']
    print(f"  Average similarity: {selectivity['avg_similarity']:.2f}")
    print(f"  Level: {selectivity['selectivity_level']}")

    # Display recommendations
    if results['recommendations']:
        print(f"\n5. Recommendations:")
        print("-" * 60)
        for i, rec in enumerate(results['recommendations'], 1):
            print(f"  {i}. {rec}")

    # Generate visualizations
    print(f"\n6. Generating visualizations...")

    if results['comparison_matrix']:
        plot_similarity_matrix(
            results['comparison_matrix'],
            results['target_names'],
            f"{output_dir}/similarity_matrix.png"
        )
        print("  ✓ Similarity matrix")

    if results['pockets']:
        plot_pocket_properties(
            results['pockets'],
            f"{output_dir}/pocket_properties.png"
        )
        print("  ✓ Pocket properties")

    if results['comparisons']:
        plot_comparison_summary(
            results['comparisons'],
            results['selectivity'],
            f"{output_dir}/comparison_summary.png"
        )
        print("  ✓ Comparison summary")

    print(f"\n✅ Pocket comparison analysis complete!")
    print(f"📁 Results saved to: {output_dir}")

    # Verify files exist
    print("\n7. Verifying output files...")
    expected_files = [
        'pocket_comparison.json',
        'repurposing_recommendations.txt',
        'similarity_matrix.png',
        'pocket_properties.png',
        'comparison_summary.png'
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
    test_pocket_comparison()
