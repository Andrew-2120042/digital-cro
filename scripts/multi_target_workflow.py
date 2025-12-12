"""
Multi-Target Screening Command-Line Interface

Usage:
    python scripts/multi_target_workflow.py \
        --targets 1HSG,3CL5,1M17 \
        --library data/molecules.csv \
        --output results/multi_target
"""

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.utils.multi_target import run_multi_target_screening


def main():
    parser = argparse.ArgumentParser(
        description='Multi-Target Molecular Screening',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Screen against 3 targets
  python scripts/multi_target_workflow.py \\
    --targets 1HSG,3CL5,1M17 \\
    --library data/molecules.csv \\
    --output results/multi_target

  # Run targets in parallel (faster, more RAM)
  python scripts/multi_target_workflow.py \\
    --targets 1HSG,3CL5 \\
    --library data/molecules.csv \\
    --parallel-targets \\
    --workers 2
        """
    )

    parser.add_argument(
        '--targets',
        type=str,
        required=True,
        help='Comma-separated PDB IDs (e.g., 1HSG,3CL5,1M17)'
    )

    parser.add_argument(
        '--library',
        type=str,
        required=True,
        help='Path to molecule library CSV or SMI file'
    )

    parser.add_argument(
        '--output',
        type=str,
        default='data/outputs/multi_target',
        help='Output directory (default: data/outputs/multi_target)'
    )

    parser.add_argument(
        '--project',
        type=str,
        default='Multi-Target Screening',
        help='Project name for reports'
    )

    parser.add_argument(
        '--client',
        type=str,
        default='Client',
        help='Client name for reports'
    )

    parser.add_argument(
        '--threshold',
        type=float,
        default=-7.0,
        help='Binding affinity threshold in kcal/mol (default: -7.0)'
    )

    parser.add_argument(
        '--workers',
        type=int,
        default=4,
        help='Workers per target for parallel docking (default: 4)'
    )

    parser.add_argument(
        '--parallel-targets',
        action='store_true',
        help='Run targets in parallel (faster but uses more RAM)'
    )

    args = parser.parse_args()

    # Parse target list
    targets = [t.strip().upper() for t in args.targets.split(',')]

    # Validate library path
    library_path = Path(args.library)
    if not library_path.exists():
        print(f"Error: Library file not found: {args.library}")
        sys.exit(1)

    print("\n" + "="*70)
    print("MULTI-TARGET SCREENING CONFIGURATION")
    print("="*70)
    print(f"  Targets: {', '.join(targets)} ({len(targets)} total)")
    print(f"  Library: {args.library}")
    print(f"  Output: {args.output}")
    print(f"  Threshold: {args.threshold} kcal/mol")
    print(f"  Workers per target: {args.workers}")
    print(f"  Parallel targets: {'Yes' if args.parallel_targets else 'No (sequential)'}")
    print("="*70)
    print()

    # Run multi-target screening
    results = run_multi_target_screening(
        target_proteins=targets,
        ligand_library_path=str(library_path),
        output_dir=args.output,
        project_name=args.project,
        client_name=args.client,
        affinity_threshold=args.threshold,
        max_workers_per_target=args.workers,
        run_targets_parallel=args.parallel_targets
    )

    print(f"\n✅ Multi-target screening complete!")
    print(f"\nResults saved to: {results['output_dir']}")
    print(f"\nGenerated files:")
    print(f"  - multi_target_results.csv     (Comparative table)")
    print(f"  - selectivity_analysis.csv     (Selectivity classification)")
    print(f"  - multi_target_summary.json    (Summary statistics)")
    print(f"  - target_*/                    (Per-target results)")
    print()


if __name__ == '__main__':
    main()
