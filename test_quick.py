#!/usr/bin/env python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

print("Testing batch docking...")

from scripts.utils.docking_batch import run_batch_docking
from scripts.utils.vina_wrapper import get_vina_box_from_pocket

output_dir = Path('data/outputs/phase3_test')
receptor = output_dir / 'proteins/1HSG_receptor.pdbqt'
ligand_pdbqt = output_dir / 'ligands/aspirin.pdbqt'

print(f"Receptor exists: {receptor.exists()}")
print(f"Ligand exists: {ligand_pdbqt.exists()}")

if not receptor.exists() or not ligand_pdbqt.exists():
    print("ERROR: Files not found")
    sys.exit(1)

pocket = {
    'center': [13.072668075561523, 22.467355728149414, 5.55748987197876],
    'volume': 1965.136962890625
}

ligands = [{
    'id': 'aspirin',
    'smiles': 'CC(=O)Oc1ccccc1C(=O)O',
    'pdbqt_path': str(ligand_pdbqt)
}]

box = get_vina_box_from_pocket(pocket)
print(f"Box: center=({box['center_x']:.1f},{box['center_y']:.1f},{box['center_z']:.1f}), size=({box['size_x']:.1f},{box['size_y']:.1f},{box['size_z']:.1f})")

smoke_dir = output_dir / 'smoke_test'
print(f"Running docking to: {smoke_dir}")

result = run_batch_docking(
    receptor_pdbqt=str(receptor),
    ligands=ligands,
    output_dir=str(smoke_dir),
    docking_box=box,
    max_workers=1,
    exhaustiveness=8
)

print(f"\nResults: {result['successful']}/{result['total_ligands']} successful")

if result['successful'] > 0:
    for r in result['results']:
        if r['success']:
            print(f"✓ Affinity: {r['binding_affinity']:.2f} kcal/mol")
            print("✅ SMOKE TEST PASSED")
else:
    print("✗ FAILED")
    for r in result['results']:
        print(f"  Error: {r['error']}")
    sys.exit(1)
