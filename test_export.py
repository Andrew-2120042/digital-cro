"""
Test export formats functionality
"""

from pathlib import Path
from scripts.utils.export_formats import export_all_formats

def test_export():
    """Test exporting docking results to all formats"""

    # Test with complete_workflow results
    result_path = Path("data/outputs/complete_workflow")

    # Find protein
    protein_pdb = list(result_path.glob("proteins/*_cleaned.pdb"))[0]
    print(f"✓ Found protein: {protein_pdb}")

    # Find ligands (top 3 for quick test)
    ligand_pdbqts = sorted((result_path / "docking" / "results").glob("*_docked.pdbqt"))[:3]
    print(f"✓ Found {len(ligand_pdbqts)} ligands")

    # Get ligand names
    ligand_names = [f.stem.replace('_docked', '') for f in ligand_pdbqts]
    print(f"✓ Ligand names: {ligand_names}")

    # Export to all formats
    print("\n🔄 Exporting to all formats...")

    export_dir = result_path / "exports_test"

    exports = export_all_formats(
        protein_pdb=str(protein_pdb),
        ligand_pdbqts=[str(p) for p in ligand_pdbqts],
        output_dir=str(export_dir),
        ligand_names=ligand_names
    )

    print("\n✅ Export complete!")
    print(f"\nExported files:")
    for format_name, files in exports.items():
        print(f"\n{format_name.upper()}:")
        for f in files:
            print(f"  - {Path(f).name}")

    # Verify files exist
    print("\n🔍 Verifying files...")
    all_good = True
    for format_name, files in exports.items():
        for f in files:
            if not Path(f).exists():
                print(f"❌ Missing: {f}")
                all_good = False

    if all_good:
        print("✅ All files created successfully!")
    else:
        print("❌ Some files are missing!")

    return exports

if __name__ == "__main__":
    test_export()
