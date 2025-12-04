# Digital CRO - Computational Drug Discovery Platform

A professional computational drug discovery service designed to help biotech companies identify promising drug candidates through systematic molecular screening and analysis.

## Overview

Digital CRO provides end-to-end computational tools for early-stage drug discovery, focusing on:
- Molecular property calculation and analysis
- Drug-likeness evaluation (Lipinski's Rule of Five)
- Protein-ligand docking simulations
- Virtual screening of compound libraries
- Lead optimization and structure-activity relationship analysis

This repository implements **Phase 1** of the platform: foundational infrastructure for molecule operations and property calculations.

## Project Structure

```
drug-discovery/
├── README.md                  # Project documentation
├── requirements.txt           # Python dependencies
├── data/                      # Data storage
│   ├── molecules/            # Input molecule files (SMILES, SDF)
│   ├── proteins/             # Protein structure files (PDB)
│   └── outputs/              # Generated results and visualizations
├── scripts/
│   ├── utils/                # Core utility modules
│   │   ├── __init__.py      # Package initialization
│   │   ├── molecule_ops.py  # Molecule loading, conversion, visualization
│   │   └── properties.py    # Property calculations and filtering
│   └── tests/               # Validation and testing scripts
│       └── test_phase1.py   # Phase 1 validation suite
└── notebooks/               # Jupyter notebooks for exploration
```

## Features

### Phase 1 Capabilities

#### Molecule Operations (`molecule_ops.py`)
- **SMILES Loading**: Parse and load molecular structures from SMILES files
- **Format Conversion**: Convert SMILES strings to RDKit molecule objects
- **2D Visualization**: Generate publication-quality 2D molecular structures
- **Grid Display**: Create organized grids of multiple molecules

#### Property Calculations (`properties.py`)
- **Lipinski Properties**: MW, LogP, HBD, HBA, rotatable bonds
- **Drug-likeness Filtering**: Automated Lipinski Rule of Five evaluation
- **Comprehensive Descriptors**: TPSA, molecular formula, ring counts, stereochemistry
- **Batch Processing**: Calculate properties for entire compound libraries

## Installation

### Prerequisites

- Python 3.11+
- macOS, Linux, or Windows
- 2GB free disk space
- Internet connection (for initial package installation)

### Step 1: Clone or Download

```bash
cd /path/to/your/workspace
# Repository already exists at: drug-discovery/
```

### Step 2: Create Virtual Environment (Recommended)

```bash
cd drug-discovery
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### Step 3: Install Dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

**Note**: RDKit installation may take 5-10 minutes as it's a large package.

### Step 4: Verify Installation

```bash
python -c "import rdkit; print(f'RDKit version: {rdkit.__version__}')"
python -c "import pandas; print(f'Pandas version: {pandas.__version__}')"
```

## Usage

### Running Phase 1 Validation Tests

The validation script tests all core functionality and generates sample outputs:

```bash
cd drug-discovery
python scripts/tests/test_phase1.py
```

**Expected Output:**
```
================================================================================
PHASE 1 VALIDATION TEST - Digital CRO
================================================================================

[1] Loading test molecules...
    ✓ Loaded 10 molecules
    ✓ 10 valid molecules, 0 failed

[2] Calculating molecular properties...
    ✓ Calculated properties for 10 molecules

[3] Filtering by Lipinski's Rule of Five...
    ✓ 7 molecules pass Lipinski rules
    ✓ 3 molecules fail Lipinski rules

[4] Generating individual molecule images...
    ✓ Generated 10 individual molecule images

[5] Generating molecule grid...
    ✓ Generated grid with 10 molecules

[6] Saving properties to CSV...
    ✓ Saved properties for 10 molecules

[7] Generating property distribution plots...
    ✓ Generated property distribution plots

[8] Summary Statistics:
...
```

### Generated Outputs

After running the test script, check `data/outputs/phase1_test/`:

1. **individual_molecules/** - PNG images of each molecule
   - `Aspirin.png`
   - `Caffeine.png`
   - `Ibuprofen.png`
   - ... (10 total)

2. **molecule_grid.png** - Grid visualization of all 10 test molecules

3. **properties.csv** - Complete property data for all molecules
   ```csv
   name,molecular_weight,logp,hbd,hba,rotatable_bonds,passes_lipinski,...
   Aspirin,180.16,1.19,1,4,3,True,...
   ```

4. **property_distributions.png** - Statistical distribution plots
   - Molecular weight histogram
   - LogP distribution
   - H-bond donor/acceptor counts
   - TPSA distribution
   - Lipinski pass/fail summary

### Using the Utility Modules

#### Example 1: Load and Analyze Molecules

```python
from scripts.utils.molecule_ops import smiles_to_mol, draw_molecule
from scripts.utils.properties import calculate_all_properties

# Convert SMILES to molecule
smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
mol = smiles_to_mol(smiles)

# Calculate properties
props = calculate_all_properties(mol)
print(f"Molecular Weight: {props['molecular_weight']:.2f}")
print(f"LogP: {props['logp']:.2f}")
print(f"Passes Lipinski: {props['passes_lipinski']}")

# Visualize
draw_molecule(mol, filename="data/outputs/aspirin.png")
```

#### Example 2: Batch Processing

```python
import pandas as pd
from scripts.utils.molecule_ops import smiles_to_mol
from scripts.utils.properties import calculate_all_properties

# Load SMILES file
smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]
names = ["Ethanol", "Acetic Acid", "Benzene"]

# Process batch
results = []
for name, smi in zip(names, smiles_list):
    mol = smiles_to_mol(smi)
    if mol:
        props = calculate_all_properties(mol)
        props['name'] = name
        results.append(props)

# Create DataFrame
df = pd.DataFrame(results)
df.to_csv("data/outputs/batch_results.csv", index=False)
```

#### Example 3: Filter Drug-like Molecules

```python
from scripts.utils.molecule_ops import smiles_to_mol
from scripts.utils.properties import passes_lipinski

compounds = {
    "Drug A": "CC(=O)Oc1ccccc1C(=O)O",
    "Drug B": "CCCCCCCCCCCCCCCCCC",  # Too lipophilic
}

for name, smiles in compounds.items():
    mol = smiles_to_mol(smiles)
    if mol and passes_lipinski(mol):
        print(f"{name}: PASS - Drug-like")
    else:
        print(f"{name}: FAIL - Not drug-like")
```

## Dependencies

Core packages and their purposes:

| Package | Version | Purpose |
|---------|---------|---------|
| rdkit | ≥2023.9.0 | Cheminformatics toolkit (molecule handling, properties) |
| pandas | ≥2.1.0 | Data manipulation and analysis |
| numpy | ≥1.24.0 | Numerical computing |
| matplotlib | ≥3.8.0 | Data visualization and plotting |
| biopython | ≥1.83 | Biological computation (protein handling) |
| pillow | ≥10.0.0 | Image processing |
| jupyter | ≥1.0.0 | Interactive notebooks |

## API Reference

### molecule_ops.py

#### `load_smiles_file(filepath: str) -> pd.DataFrame`
Load SMILES file into DataFrame. Supports comma, tab, and space-separated formats.

**Parameters:**
- `filepath`: Path to SMILES file

**Returns:** DataFrame with columns `['smiles', 'mol_id', ...]`

---

#### `smiles_to_mol(smiles: str) -> Optional[Chem.Mol]`
Convert SMILES string to RDKit molecule object.

**Parameters:**
- `smiles`: SMILES string

**Returns:** RDKit Mol object or None if invalid

---

#### `draw_molecule(mol, size=(300,300), filename=None)`
Draw 2D structure of molecule.

**Parameters:**
- `mol`: RDKit Mol object
- `size`: Image dimensions (width, height)
- `filename`: Optional save path (PNG)

**Returns:** PIL Image if no filename, else None

---

#### `draw_molecule_grid(mols, labels=None, mols_per_row=4, filename=None)`
Draw grid of multiple molecules.

**Parameters:**
- `mols`: List of RDKit Mol objects
- `labels`: Optional list of labels
- `mols_per_row`: Molecules per row
- `filename`: Optional save path (PNG)

**Returns:** PIL Image if no filename, else None

### properties.py

#### `calculate_lipinski_properties(mol) -> dict`
Calculate Lipinski Rule of Five properties.

**Returns:** Dictionary with `molecular_weight`, `logp`, `hbd`, `hba`, `rotatable_bonds`

---

#### `passes_lipinski(mol) -> bool`
Check if molecule passes Lipinski criteria.

**Criteria:**
- MW: 150-500 Da
- LogP: -0.4 to 5.6
- HBD: ≤5
- HBA: ≤10

**Returns:** True if passes, False otherwise

---

#### `calculate_all_properties(mol) -> dict`
Calculate comprehensive molecular properties.

**Returns:** Dictionary with 15+ properties including Lipinski parameters, TPSA, ring counts, formula, stereochemistry

## Testing

### Test Molecules

The Phase 1 validation uses 10 known pharmaceutical compounds:

1. **Aspirin** - Analgesic/anti-inflammatory
2. **Caffeine** - Stimulant
3. **Ibuprofen** - NSAID
4. **Paracetamol** (Acetaminophen) - Analgesic
5. **Penicillin** - Antibiotic
6. **Atorvastatin** - Statin (cholesterol)
7. **Metformin** - Antidiabetic
8. **Warfarin** - Anticoagulant
9. **Sildenafil** (Viagra) - PDE5 inhibitor
10. **Lisinopril** - ACE inhibitor

### Validation Criteria

Phase 1 is successfully validated if:
- ✓ All 10 molecules load without errors
- ✓ Properties calculated for all molecules
- ✓ Individual PNG images generated
- ✓ Grid visualization created
- ✓ Properties CSV file saved
- ✓ Distribution plots generated
- ✓ Summary statistics printed
- ✓ Lipinski filtering works correctly

## Troubleshooting

### RDKit Installation Issues

**Issue:** RDKit fails to install via pip

**Solution:**
```bash
# Try conda instead
conda create -n drug-discovery python=3.11
conda activate drug-discovery
conda install -c conda-forge rdkit
pip install -r requirements.txt
```

### Import Errors

**Issue:** `ModuleNotFoundError: No module named 'scripts'`

**Solution:** Run scripts from project root:
```bash
cd drug-discovery
python scripts/tests/test_phase1.py  # Correct
# NOT: cd scripts/tests && python test_phase1.py
```

### Permission Errors

**Issue:** Cannot write to `data/outputs/`

**Solution:**
```bash
chmod -R 755 data/
```

## Roadmap

### Phase 2: Protein-Ligand Docking (Planned)
- PDB file parsing and protein preparation
- Docking grid generation
- AutoDock Vina integration
- Binding affinity prediction

### Phase 3: Virtual Screening (Planned)
- High-throughput docking
- Parallel processing
- Hit identification and ranking
- Interactive result visualization

### Phase 4: Lead Optimization (Planned)
- Structure-activity relationship (SAR) analysis
- Similarity searching
- Scaffold hopping
- ADMET prediction

### Phase 5: Reporting & Deployment (Planned)
- Automated report generation
- RESTful API
- Web interface
- Cloud deployment

## Contributing

This is a production research platform. Code standards:
- Full docstrings for all functions
- Type hints for function signatures
- Comprehensive error handling
- Unit tests for new features
- PEP 8 style compliance

## License

Proprietary - Digital CRO Platform
All rights reserved.

## Contact

For questions or support regarding Digital CRO:
- Technical issues: Check troubleshooting section
- Feature requests: Document in project roadmap
- Collaboration: Contact project maintainers

## Acknowledgments

Built with:
- **RDKit** - Open-source cheminformatics toolkit
- **Python** - Scientific computing ecosystem
- **Lipinski's Rule of Five** - Drug-likeness criteria

---

**Version:** 1.0.0 (Phase 1)
**Last Updated:** 2025
**Status:** Production-ready foundation
