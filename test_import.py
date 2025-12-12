#!/usr/bin/env python
"""Minimal import test for Phase 4"""
import sys
from pathlib import Path

print("=" * 70)
print("PHASE 4 IMPORT TEST")
print("=" * 70)

# Test 1: Import admet_predictions
print("\n[1] Testing admet_predictions import...")
try:
    from scripts.utils.admet_predictions import (
        calculate_lipinski_properties,
        predict_admet_properties,
        calculate_qed,
        predict_bbb_penetration
    )
    print("✓ admet_predictions imports successfully")
except ImportError as e:
    print(f"✗ Import failed: {e}")
    sys.exit(1)

# Test 2: Import molecular_viz
print("\n[2] Testing molecular_viz import...")
try:
    from scripts.utils.molecular_viz import (
        draw_molecules_grid,
        create_admet_radar_plot,
        plot_admet_summary
    )
    print("✓ molecular_viz imports successfully")
except ImportError as e:
    print(f"✗ Import failed: {e}")
    sys.exit(1)

# Test 3: Quick functionality test
print("\n[3] Testing basic functionality...")
try:
    smiles = 'CC(=O)Oc1ccccc1C(=O)O'  # Aspirin
    props = calculate_lipinski_properties(smiles)
    print(f"✓ Aspirin MW: {props['molecular_weight']:.2f} Da")
    print(f"✓ Aspirin LogP: {props['logp']:.2f}")
    print(f"✓ Lipinski compliant: {props['lipinski_compliant']}")
except Exception as e:
    print(f"✗ Function test failed: {e}")
    sys.exit(1)

# Test 4: ADMET prediction
print("\n[4] Testing ADMET prediction...")
try:
    admet = predict_admet_properties(smiles)
    print(f"✓ QED Score: {admet['qed_score']:.3f}")
    print(f"✓ SA Score: {admet['sa_score']:.2f}")
    print(f"✓ BBB Penetrant: {admet['bbb_penetrant']}")
    print(f"✓ Oral Bioavailability: {admet['oral_bioavailability']}")
except Exception as e:
    print(f"✗ ADMET test failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 70)
print("✅ ALL TESTS PASSED!")
print("=" * 70)
print("\nPhase 4 modules are working correctly.")
print("Ready to run full test suite: python scripts/tests/test_phase4.py")
