# Digital CRO Platform - Streamlit Web Interface

## Quick Start

### Run the Application

```bash
cd "/Users/nareshnallabothula/digital cro/drug-discovery"
streamlit run app.py
```

The application will automatically open in your browser at: **http://localhost:8501**

---

## What You'll See

### 🏠 Home Page
- Platform overview and capabilities
- Quick start guide
- Platform statistics

### 🔬 New Screening
**Run a virtual screening campaign:**
1. Enter project details (name, client, parameters)
2. Input protein:
   - Option A: PDB ID (e.g., `1HSG`)
   - Option B: Upload PDB file
3. Upload molecule library (CSV with SMILES)
4. Click "Run Screening"
5. Watch real-time progress
6. Download results (PDF, CSV, ZIP)

### 📊 Results Browser
**View past screenings:**
- Browse all completed campaigns
- View detailed results and metrics
- Download previous reports
- See visualizations

### 📖 Documentation
**Complete user guide:**
- How to use the platform
- Interpret results
- File formats
- FAQ and troubleshooting

---

## Test It Now

### Quick Test Run

1. **Start Streamlit:**
   ```bash
   streamlit run app.py
   ```

2. **Navigate to "🔬 New Screening"**

3. **Fill in the form:**
   - Project Name: `HIV-1 Test Screening`
   - Client Name: `Test Lab`
   - PDB ID: `1HSG`
   - Upload library: `data/outputs/phase3_test/test_library.smi`

4. **Click "Run Screening"**

5. **Watch the magic happen!**
   - Real-time progress bar
   - Live log output
   - Results appear when complete

6. **Download your results:**
   - PDF report
   - CSV data
   - Complete ZIP archive

---

## Features

✅ **No Coding Required** - Browser-based interface
✅ **Real-Time Progress** - Live updates during screening
✅ **Professional Reports** - PDF deliverables
✅ **Data Export** - CSV and ZIP downloads
✅ **Results Browser** - View all past campaigns
✅ **Comprehensive Docs** - Built-in help system

---

## File Structure

```
app.py                      # Main Streamlit application
pages/
├── home.py                # Home page
├── screening.py           # Screening workflow
├── results.py             # Results browser
└── documentation.py       # Documentation
```

---

## Platform Capabilities

### Molecular Docking
- AutoDock Vina integration
- Automated protein preparation
- Pocket detection
- Parallel processing

### ADMET Predictions
- Lipinski's Rule of Five
- Drug-likeness (QED)
- BBB penetration
- Oral bioavailability

### Deliverables
- Professional PDF reports
- Complete CSV data
- Visualization figures
- Docked structures

---

## Deployment Options

### Local (Current)
```bash
streamlit run app.py
```

### Network Access
```bash
streamlit run app.py --server.address 0.0.0.0
```
Access from any device: `http://your-ip:8501`

### Streamlit Cloud (Free Hosting)
1. Push to GitHub
2. Connect to Streamlit Cloud
3. Deploy with one click
4. Share public URL

### Docker
```dockerfile
FROM python:3.10-slim
WORKDIR /app
COPY . .
RUN pip install -r requirements.txt
EXPOSE 8501
CMD ["streamlit", "run", "app.py", "--server.port=8501"]
```

---

## Troubleshooting

### Port Already in Use
```bash
streamlit run app.py --server.port 8502
```

### Import Errors
```bash
pip install -r requirements.txt
```

### Screening Fails
- Check AutoDock Vina is installed: `vina --version`
- Verify valid SMILES in library
- Ensure sufficient disk space

---

## Next Steps

1. ✅ **Test the interface** (you are here!)
2. 🎨 **Customize branding** in `app.py`
3. 🚀 **Deploy to production**
4. 💰 **Start offering to clients**

---

## Support

For detailed documentation, see:
- `PHASE6A_STREAMLIT_GUIDE.md` - Complete setup guide
- In-app documentation page (📖 Documentation)
- Phase 1-5 documentation for backend details

---

## License

© 2024 Digital CRO Platform. All rights reserved.

**Version:** 1.0.0 (Phase 6A)
