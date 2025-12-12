"""
Quick smoke test for docking functionality.
Uses existing Phase 3.2 test data to verify batch docking works.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from scripts.utils.docking_batch import run_batch_docking
from scripts.utils.vina_wrapper import get_vina_box_from_pocket

# Use existing Phase 3.2 test data
output_dir = Path('data/outputs/phase3_test')
receptor = output_dir / 'proteins/1HSG_receptor.pdbqt'
ligand_dir = output_dir / 'ligands'

# Pocket from Phase 3.2
pocket = {
    'center': [13.072668075561523, 22.467355728149414, 5.55748987197876],
    'volume': 1965.136962890625
}

# Build ligand list (just aspirin for smoke test)
ligands = [{
    'id': 'aspirin',
    'smiles': 'CC(=O)Oc1ccccc1C(=O)O',
    'pdbqt_path': str(ligand_dir / 'aspirin.pdbqt')
}]

print("Smoke Test: Single Ligand Docking")
print(f"Receptor: {receptor}")
print(f"Ligand: {ligands[0]['pdbqt_path']}")

# Check files exist
if not receptor.exists():
    print(f"ERROR: Receptor not found: {receptor}")
    sys.exit(1)

if not Path(ligands[0]['pdbqt_path']).exists():
    print(f"ERROR: Ligand not found: {ligands[0]['pdbqt_path']}")
    sys.exit(1)

# Get box parameters
box = get_vina_box_from_pocket(pocket)
print(f"Box center: ({box['center_x']:.2f}, {box['center_y']:.2f}, {box['center_z']:.2f})")
print(f"Box size: ({box['size_x']:.1f}, {box['size_y']:.1f}, {box['size_z']:.1f})")

# Run docking
smoke_test_dir = output_dir / 'smoke_test'
result = run_batch_docking(
    receptor_pdbqt=str(receptor),
    ligands=ligands,
    output_dir=str(smoke_test_dir),
    docking_box=box,
    max_workers=1,
    exhaustiveness=8
)

print("\nResults:")
print(f"Success: {result['successful']}/1")
print(f"Failed: {result['failed']}/1")

if result['successful'] > 0:
    for r in result['results']:
        if r['success']:
            print(f"\n✓ Aspirin docked successfully!")
            print(f"  Binding affinity: {r['binding_affinity']:.2f} kcal/mol")
            print(f"  Runtime: {r['runtime']:.1f}s")
    print("\n✅ SMOKE TEST PASSED")
    sys.exit(0)
else:
    print("\n✗ SMOKE TEST FAILED")
    for r in result['results']:
        print(f"  Error: {r['error']}")
    sys.exit(1)
