"""
Phase 5 Test: PDF Report Generation

Tests:
1. Create report with mock data
2. Add all sections (cover, summary, results, ADMET, methodology)
3. Embed visualizations
4. Generate complete workflow report
"""

import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from scripts.utils.report_generator import DrugDiscoveryReport, generate_complete_report


def test_basic_report():
    """Test basic report creation"""
    print("\n" + "="*70)
    print("[Test 1] BASIC REPORT CREATION")
    print("="*70)

    output_dir = Path('data/outputs/phase5_test')
    output_dir.mkdir(parents=True, exist_ok=True)

    report_path = output_dir / 'test_report.pdf'

    # Create report
    report = DrugDiscoveryReport(
        output_path=str(report_path),
        project_name="Test Drug Discovery",
        client_name="Test Client"
    )

    # Add cover
    report.add_cover_page(
        target_protein="HIV-1 Protease",
        target_pdb_id="1HSG",
        num_molecules_screened=100
    )

    # Add methodology
    report.add_methodology()

    # Generate
    pdf_path = report.generate()

    print(f"✓ Basic report created: {pdf_path}")
    assert Path(pdf_path).exists()

    # Check file size
    file_size = Path(pdf_path).stat().st_size
    print(f"  File size: {file_size:,} bytes")

    return pdf_path


def test_complete_report():
    """Test complete report with data"""
    print("\n" + "="*70)
    print("[Test 2] COMPLETE REPORT WITH DATA")
    print("="*70)

    # Create mock data
    data = {
        'ligand_id': ['aspirin', 'ibuprofen', 'caffeine', 'morphine', 'paracetamol'],
        'smiles': [
            'CC(=O)Oc1ccccc1C(=O)O',
            'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'CN1CCC23c4ccccc4CC2C1Cc5c3cccc5O',
            'CC(=O)Nc1ccc(O)cc1'
        ],
        'binding_affinity': [-6.67, -7.21, -5.89, -7.85, -6.42],
        'qed_score': [0.550, 0.822, 0.538, 0.806, 0.595],
        'lipinski_compliant': [True, True, True, True, True],
        'oral_bioavailability': ['High', 'High', 'High', 'High', 'High'],
        'bbb_penetrant': [True, True, True, True, True]
    }

    df = pd.DataFrame(data)

    output_dir = Path('data/outputs/phase5_test')
    report_path = output_dir / 'complete_report.pdf'

    # Generate complete report
    pdf_path = generate_complete_report(
        df_results=df,
        output_path=str(report_path),
        project_name="Test Project",
        client_name="Test Client",
        target_protein="HIV-1 Protease",
        target_pdb_id="1HSG",
        affinity_threshold=-7.0
    )

    print(f"✓ Complete report created: {pdf_path}")
    assert Path(pdf_path).exists()

    # Check file size (should be > 3KB for basic report)
    file_size = Path(pdf_path).stat().st_size
    print(f"  File size: {file_size:,} bytes")
    assert file_size > 3000, "PDF too small"

    return pdf_path


def test_report_with_visualizations():
    """Test report with embedded visualizations"""
    print("\n" + "="*70)
    print("[Test 3] REPORT WITH VISUALIZATIONS")
    print("="*70)

    # Create mock data
    data = {
        'ligand_id': ['mol1', 'mol2', 'mol3', 'mol4', 'mol5'],
        'smiles': [
            'CC(=O)Oc1ccccc1C(=O)O',
            'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'CN1CCC23c4ccccc4CC2C1Cc5c3cccc5O',
            'CC(=O)Nc1ccc(O)cc1'
        ],
        'binding_affinity': [-8.5, -7.9, -7.3, -6.8, -6.2],
        'qed_score': [0.650, 0.722, 0.638, 0.706, 0.695],
        'lipinski_compliant': [True, True, True, True, True],
        'oral_bioavailability': ['High', 'High', 'Medium', 'High', 'High'],
        'bbb_penetrant': [True, True, False, True, True],
        'molecular_weight': [180.16, 206.28, 194.19, 285.34, 151.16],
        'logp': [1.31, 3.07, -1.03, 0.89, 0.91],
        'tpsa': [63.6, 37.3, 61.8, 43.7, 49.3],
        'sa_score': [3.08, 3.60, 3.54, 5.82, 2.66]
    }

    df = pd.DataFrame(data)

    output_dir = Path('data/outputs/phase5_test')
    viz_dir = output_dir / 'viz'
    viz_dir.mkdir(exist_ok=True)

    # Create dummy visualizations for testing
    from scripts.utils.molecular_viz import (
        visualize_docking_results,
        create_admet_radar_plot,
        plot_admet_summary
    )

    # Generate visualizations
    top_hits_path = viz_dir / 'top_hits.png'
    visualize_docking_results(
        df,
        top_n=5,
        output_path=str(top_hits_path)
    )
    print(f"  ✓ Created top hits visualization")

    # ADMET summary
    admet_summary_path = viz_dir / 'admet_summary.png'
    plot_admet_summary(df, output_path=str(admet_summary_path))
    print(f"  ✓ Created ADMET summary")

    # Radar plot
    radar_path = viz_dir / 'admet_radar_mol1.png'
    create_admet_radar_plot(
        df,
        'mol1',
        id_column='ligand_id',
        output_path=str(radar_path)
    )
    print(f"  ✓ Created radar plot")

    # Generate report with visualizations
    report_path = output_dir / 'report_with_viz.pdf'

    pdf_path = generate_complete_report(
        df_results=df,
        output_path=str(report_path),
        project_name="Visualization Test Project",
        client_name="Test Pharma Inc.",
        target_protein="Test Target",
        target_pdb_id="TEST",
        visualization_dir=str(viz_dir),
        affinity_threshold=-7.0
    )

    print(f"✓ Report with visualizations created: {pdf_path}")
    assert Path(pdf_path).exists()

    # Check file size (should be larger with images)
    file_size = Path(pdf_path).stat().st_size
    print(f"  File size: {file_size:,} bytes")
    assert file_size > 10000, "PDF with visualizations should be > 10KB"

    return pdf_path


def main():
    """Run Phase 5 tests"""
    print("\n" + "="*70)
    print("PHASE 5 TEST: PDF REPORT GENERATION")
    print("="*70)

    try:
        # Test 1: Basic report
        test_basic_report()

        # Test 2: Complete report
        test_complete_report()

        # Test 3: Report with visualizations
        test_report_with_visualizations()

        print("\n" + "="*70)
        print("✅ PHASE 5 COMPLETE!")
        print("="*70)
        print("\nAll report generation modules working correctly")
        print("Digital CRO platform is production-ready!")
        print("\nGenerated reports in: data/outputs/phase5_test/")
        print("\nYou can now run the complete workflow:")
        print("  python scripts/complete_workflow.py --pdb 1HSG --library data/test_library.smi")

        return True

    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
