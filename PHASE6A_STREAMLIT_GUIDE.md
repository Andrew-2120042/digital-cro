# Phase 6A: Streamlit Web Interface - Setup Guide

## Overview

Phase 6A adds a professional web-based user interface to the Digital CRO platform using Streamlit. This allows clients to run virtual screening campaigns through their browser without any coding knowledge.

## Features

### 🏠 Home Page
- Platform overview and capabilities
- Quick start guide
- Platform statistics
- Call-to-action buttons

### 🔬 New Screening
- Interactive form for project configuration
- Protein input (PDB ID or file upload)
- Molecule library upload (CSV/SMI format)
- Parameter configuration
- Real-time progress tracking
- Results summary with visualizations
- Download buttons (PDF, CSV, ZIP)

### 📊 Results Browser
- Browse all past screening campaigns
- View detailed results and metrics
- Filter and sort molecules
- Display visualizations
- Download previous results

### 📖 Documentation
- Comprehensive user guide
- File format specifications
- Best practices
- FAQ section
- Troubleshooting tips

## Installation

Streamlit has been added to requirements.txt:

```bash
pip install streamlit>=1.28.0
```

All dependencies should already be installed if you followed Phase 1-5 setup.

## Running the Application

### Method 1: Command Line

```bash
cd "/Users/nareshnallabothula/digital cro/drug-discovery"
streamlit run app.py
```

The application will open in your default browser at: http://localhost:8501

### Method 2: Background Mode

```bash
streamlit run app.py --server.headless true --server.port 8501
```

### Method 3: Custom Port

```bash
streamlit run app.py --server.port 8080
```

## Usage Guide

### Running a New Screening

1. **Navigate to "🔬 New Screening"**

2. **Fill in Project Configuration:**
   - Project Name: e.g., "HIV-1 Protease Screening"
   - Client Name: e.g., "Biotech Research Corp"
   - Affinity Threshold: Default -7.0 kcal/mol
   - Parallel Workers: 4-8 recommended

3. **Select Protein Input:**
   - **Option A:** Enter PDB ID (e.g., `1HSG`)
   - **Option B:** Upload PDB file

4. **Upload Molecule Library:**
   - CSV format with `smiles` column
   - Optional `id` column for molecule names
   - Preview shows first 10 molecules

5. **Click "Run Screening"**
   - Progress bar shows real-time status
   - Detailed logs available in expander
   - Results appear when complete

6. **Download Results:**
   - 📄 PDF Report
   - 📊 CSV Data
   - 📦 Complete ZIP archive

### Viewing Past Results

1. **Navigate to "📊 Results Browser"**

2. **Select a screening campaign** from dropdown

3. **View results in tabs:**
   - 📋 All Results - Full dataset with filters
   - 🏆 Top Hits - Filtered and ranked
   - 📊 Visualizations - All generated figures
   - 📥 Downloads - PDF, CSV, ZIP files

## File Structure

```
drug-discovery/
├── app.py                          # Main Streamlit application
├── pages/
│   ├── __init__.py                # Package init
│   ├── home.py                    # Home page
│   ├── screening.py               # Screening workflow
│   ├── results.py                 # Results browser
│   └── documentation.py           # Documentation page
├── scripts/
│   ├── complete_workflow.py       # Backend workflow
│   └── utils/                     # Utility modules
├── data/
│   └── outputs/                   # Screening results
│       └── streamlit_YYYYMMDD_HHMMSS/  # Timestamped results
└── requirements.txt               # Python dependencies
```

## Output Directory Structure

Each screening creates a timestamped directory:

```
data/outputs/streamlit_20241130_143022/
├── 1HSG_report.pdf                # Professional PDF report
├── final_results.csv              # Complete results
├── proteins/
│   ├── 1HSG.pdb                  # Downloaded protein
│   ├── 1HSG_cleaned.pdb          # Cleaned protein
│   └── 1HSG_receptor.pdbqt       # Docking-ready receptor
├── ligands/
│   └── *.pdbqt                   # Prepared ligands
├── docking/
│   ├── docking_results.csv       # Summary results
│   ├── docking_results_detailed.json
│   └── results/
│       └── *_docked.pdbqt        # Docked poses
└── visualizations/
    ├── top_hits.png              # Molecule grid
    ├── admet_summary.png         # ADMET bar charts
    └── admet_radar_*.png         # Radar plots
```

## Customization

### Branding

Edit the CSS in `app.py` to customize colors and styling:

```python
.main-header {
    color: #1f4788;  # Change header color
}

.stButton>button {
    background-color: #1f4788;  # Change button color
}
```

### Adding Pages

1. Create new file in `pages/` directory
2. Implement `show()` function
3. Import and add to navigation in `app.py`

### Modifying Workflow

The screening page calls `scripts/complete_workflow.py` as a subprocess. Modify that script to change the backend behavior.

## Deployment Options

### Option 1: Streamlit Cloud (Free)

1. Push code to GitHub
2. Connect to Streamlit Cloud: https://streamlit.io/cloud
3. Deploy with one click
4. Share public URL with clients

### Option 2: Docker

```dockerfile
FROM python:3.10-slim

WORKDIR /app
COPY . /app

RUN pip install -r requirements.txt

EXPOSE 8501

CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
```

### Option 3: Local Server

```bash
streamlit run app.py --server.port 8501 --server.address 0.0.0.0
```

Access from network: http://your-ip:8501

## Troubleshooting

### Issue: Port already in use

**Solution:**
```bash
streamlit run app.py --server.port 8502
```

### Issue: Module import errors

**Solution:**
```bash
pip install -r requirements.txt
```

### Issue: Streamlit not found

**Solution:**
```bash
pip install streamlit>=1.28.0
```

### Issue: Screening fails

**Check:**
- AutoDock Vina is installed (`vina --version`)
- Valid SMILES in library file
- Valid PDB file format
- Sufficient disk space

## Performance Tips

1. **Parallel Workers:**
   - Use 50-75% of available CPU cores
   - 4-8 workers optimal for most systems
   - More workers ≠ always faster (diminishing returns)

2. **Library Size:**
   - Start small (20-100 molecules) for testing
   - Scale up after validation
   - Consider chunking large libraries (>1000 molecules)

3. **Progress Tracking:**
   - Real-time logs help monitor progress
   - Expand "Detailed Logs" to see full output
   - Check for errors during execution

## Next Steps

After Phase 6A, you can:

1. **Test the interface:**
   ```bash
   streamlit run app.py
   ```

2. **Run a test screening:**
   - Use PDB: 1HSG
   - Upload test library from `data/outputs/phase3_test/test_library.smi`

3. **Customize branding:**
   - Change colors in app.py CSS
   - Add your logo
   - Update contact information

4. **Deploy to production:**
   - Streamlit Cloud for easy hosting
   - Docker for self-hosting
   - Internal server for enterprise

## Support

For issues or questions:
- Check the 📖 Documentation page in the app
- Review this guide
- Check Phase 1-5 setup if backend issues occur

## Version History

- **v1.0.0 (Phase 6A):** Initial Streamlit interface
  - Home page
  - Screening workflow
  - Results browser
  - Documentation

## License

© 2024 Digital CRO Platform. All rights reserved.
