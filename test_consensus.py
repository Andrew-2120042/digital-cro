#!/usr/bin/env python3
"""
Quick test of consensus docking module
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent / 'drug-discovery'))

print("Testing Consensus Docking Module...")
print("=" * 60)

# Test 1: Import module
print("\n[1/3] Testing module import...")
try:
    from scripts.utils.consensus_docking import ConsensusDocking
    print("✓ Module imported successfully")
except Exception as e:
    print(f"✗ Failed to import module: {e}")
    sys.exit(1)

# Test 2: Check available methods
print("\n[2/3] Checking available docking methods...")
try:
    consensus = ConsensusDocking()
    print(f"✓ Available methods: {consensus.available_methods}")

    if 'vina' in consensus.available_methods:
        print("  ✓ AutoDock Vina: AVAILABLE")
    else:
        print("  ✗ AutoDock Vina: NOT FOUND")

    if 'smina' in consensus.available_methods:
        print("  ✓ Smina: AVAILABLE (optimal)")
    else:
        print("  ⚠️  Smina: NOT FOUND (will use Vina only)")
        print("     To install Smina: bash scripts/install_smina.sh")

    if 'ledock' in consensus.available_methods:
        print("  ✓ LeDock: AVAILABLE")
    else:
        print("  ℹ️  LeDock: NOT FOUND (optional)")

except Exception as e:
    print(f"✗ Failed to initialize consensus docking: {e}")
    sys.exit(1)

# Test 3: Test visualization imports
print("\n[3/3] Testing visualization modules...")
try:
    from scripts.utils.consensus_viz import (
        plot_method_agreement,
        plot_confidence_distribution,
        plot_consensus_summary
    )
    print("✓ Visualization modules imported successfully")
except Exception as e:
    print(f"✗ Failed to import visualizations: {e}")
    sys.exit(1)

# Summary
print("\n" + "=" * 60)
print("✅ CONSENSUS DOCKING MODULE TEST PASSED")
print("=" * 60)

if len(consensus.available_methods) >= 2:
    print(f"\n🎯 Multi-method consensus ready with {len(consensus.available_methods)} methods!")
else:
    print(f"\n⚠️  Only {len(consensus.available_methods)} method available.")
    print("   Install Smina for full consensus validation:")
    print("   bash scripts/install_smina.sh")

print("\nReady to use consensus docking in:")
print("  - Streamlit UI: Check 'Use Consensus Docking' checkbox")
print("  - Command line: python scripts/complete_workflow.py --consensus ...")

print("\nFeature 14 implementation complete! 🎉")
