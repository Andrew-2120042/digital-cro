# Phase 5: PDF Report Generation - Setup Guide

## Installation

Phase 5 requires the ReportLab library for PDF generation.

```bash
pip install reportlab
```

## Testing Phase 5

Run the Phase 5 test suite:

```bash
cd "/Users/nareshnallabothula/digital cro/drug-discovery"
python scripts/tests/test_phase5.py
```

Expected output:
```
======================================================================
PHASE 5 TEST: PDF REPORT GENERATION
======================================================================

======================================================================
[Test 1] BASIC REPORT CREATION
======================================================================
✓ Basic report created: data/outputs/phase5_test/test_report.pdf
  File size: 12,345 bytes

======================================================================
[Test 2] COMPLETE REPORT WITH DATA
======================================================================
✓ Complete report created: data/outputs/phase5_test/complete_report.pdf
  File size: 25,678 bytes

======================================================================
[Test 3] REPORT WITH VISUALIZATIONS
======================================================================
  ✓ Created top hits visualization
  ✓ Created ADMET summary
  ✓ Created radar plot
✓ Report with visualizations created: data/outputs/phase5_test/report_with_viz.pdf
  File size: 156,789 bytes

======================================================================
✅ PHASE 5 COMPLETE!
======================================================================

All report generation modules working correctly
Digital CRO platform is production-ready!
```

## Running the Complete Workflow

Once Phase 5 tests pass, you can run the complete end-to-end workflow:

```bash
python scripts/complete_workflow.py \
  --pdb 1HSG \
  --library data/outputs/phase3_test/test_library.smi \
  --output data/outputs/complete_workflow \
  --project "HIV-1 Protease Inhibitor Discovery" \
  --client "Biotech Research Corp" \
  --threshold -7.0 \
  --workers 4
```

This will:
1. Download and prepare the 1HSG protein structure
2. Detect the binding pocket
3. Prepare the receptor for docking
4. Prepare ligands from the SMILES library
5. Run batch molecular docking
6. Rank and filter hits
7. Run ADMET predictions
8. Generate professional visualizations
9. Create a client-ready PDF report

## Output Files

The workflow creates:
- `final_results.csv` - Complete results with docking + ADMET data
- `visualizations/` - All generated plots and figures
- `1HSG_report.pdf` - Professional PDF report ready for clients

## Features

### PDF Report Includes:
- ✅ Professional cover page with project details
- ✅ Executive summary with key findings
- ✅ Top hits table with rankings
- ✅ Molecular docking results section
- ✅ Embedded visualizations (molecular structures)
- ✅ ADMET analysis with compliance rates
- ✅ ADMET summary dashboard
- ✅ Radar plots for top molecules
- ✅ Methodology appendix

### Complete Workflow Features:
- ✅ Command-line interface
- ✅ Parallel docking (configurable workers)
- ✅ Comprehensive error handling
- ✅ Progress logging
- ✅ Modular and extensible
- ✅ Production-ready

## Troubleshooting

### Issue: reportlab not found
```bash
pip install reportlab
```

### Issue: Missing visualizations in PDF
Make sure Phase 4 visualizations are generated correctly. Check that matplotlib and seaborn are installed:
```bash
pip install matplotlib seaborn
```

### Issue: No successful docking results
- Verify AutoDock Vina is installed: `vina --version`
- Check ligand preparation succeeded
- Review log output for specific errors

## Next Steps

With Phase 5 complete, the Digital CRO platform is **production-ready**!

You can now:
1. Run complete drug discovery campaigns
2. Generate professional PDF reports
3. Deliver client-ready results
4. Scale to larger molecule libraries
5. Customize branding and styling

**This is what you charge $20k+ for!**
