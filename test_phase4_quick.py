#!/usr/bin/env python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

print("Testing Phase 4 modules...")

# Test imports
try:
    print("\n1. Importing admet_predictions...")
    from scripts.utils.admet_predictions import (
        calculate_lipinski_properties,
        predict_admet_properties
    )
    print("   ✓ admet_predictions imported")
except Exception as e:
    print(f"   ✗ Failed: {e}")
    sys.exit(1)

try:
    print("\n2. Importing molecular_viz...")
    from scripts.utils.molecular_viz import (
        draw_molecules_grid,
        create_admet_radar_plot
    )
    print("   ✓ molecular_viz imported")
except Exception as e:
    print(f"   ✗ Failed: {e}")
    sys.exit(1)

# Quick functional test
try:
    print("\n3. Testing Lipinski calculation...")
    smiles = 'CC(=O)Oc1ccccc1C(=O)O'  # Aspirin
    props = calculate_lipinski_properties(smiles)
    print(f"   MW: {props['molecular_weight']:.2f} Da")
    print(f"   LogP: {props['logp']:.2f}")
    print(f"   Lipinski compliant: {props['lipinski_compliant']}")
    print("   ✓ Lipinski calculation works")
except Exception as e:
    print(f"   ✗ Failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

try:
    print("\n4. Testing ADMET prediction...")
    admet = predict_admet_properties(smiles)
    print(f"   QED: {admet['qed_score']:.3f}")
    print(f"   SA: {admet['sa_score']:.2f}")
    print(f"   BBB: {admet['bbb_penetrant']}")
    print(f"   Bioavailability: {admet['oral_bioavailability']}")
    print("   ✓ ADMET prediction works")
except Exception as e:
    print(f"   ✗ Failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n✅ Phase 4 quick test PASSED!")
print("All modules working correctly")
