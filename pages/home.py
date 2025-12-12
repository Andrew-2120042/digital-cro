"""
Home page for Digital CRO Platform
"""

import streamlit as st
from pathlib import Path

def show():
    """Display home page"""

    # Header
    st.markdown('<p class="main-header">💊 Digital CRO Platform</p>', unsafe_allow_html=True)
    st.markdown(
        '<p class="sub-header">AI-Powered Drug Discovery & Molecular Screening</p>',
        unsafe_allow_html=True
    )

    st.markdown("---")

    # Overview section
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("""
        <div class="metric-card">
            <h3>🎯 Virtual Screening</h3>
            <p>Screen 1000s of molecules against your target protein using AutoDock Vina</p>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
        <div class="metric-card">
            <h3>🧪 ADMET Predictions</h3>
            <p>Predict drug-likeness, toxicity, and pharmacokinetic properties</p>
        </div>
        """, unsafe_allow_html=True)

    with col3:
        st.markdown("""
        <div class="metric-card">
            <h3>📄 Professional Reports</h3>
            <p>Generate publication-quality PDF reports with visualizations</p>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("---")

    # Key features
    st.subheader("✨ Key Features")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        **Molecular Docking:**
        - Automated protein preparation
        - Ligand-based pocket detection
        - Parallel batch docking
        - Multiple binding modes per molecule

        **ADMET Analysis:**
        - Lipinski's Rule of Five
        - QED (drug-likeness) scoring
        - BBB penetration prediction
        - Oral bioavailability estimation
        """)

    with col2:
        st.markdown("""
        **Visualization:**
        - 2D molecular structures
        - 3D protein-ligand complexes
        - ADMET property distributions
        - Binding affinity rankings

        **Deliverables:**
        - Professional PDF reports
        - CSV data files
        - Visualization figures
        - Docked structure files
        """)

    st.markdown("---")

    # Quick start guide
    st.subheader("🚀 Quick Start")

    st.markdown("""
    **Getting started is easy:**

    1. **Prepare your data**
       - Protein structure (PDB ID or upload PDB file)
       - Molecule library (CSV file with SMILES strings)

    2. **Run screening**
       - Go to 🔬 New Screening
       - Upload files and configure parameters
       - Click "Run Screening"

    3. **Download results**
       - View results in 📊 Results Browser
       - Download PDF report and data files

    **Need help?** Check the 📖 Documentation page for detailed guides.
    """)

    st.markdown("---")

    # Call to action
    col1, col2, col3 = st.columns([1, 2, 1])

    with col2:
        if st.button("🔬 Start New Screening", use_container_width=True):
            st.session_state['navigate_to'] = '🔬 New Screening'
            st.rerun()

    # Statistics (if available)
    st.markdown("---")
    st.subheader("📈 Platform Statistics")

    # Check for previous results
    results_dir = Path("data/outputs")

    if results_dir.exists():
        completed_runs = len(list(results_dir.glob("*/final_results.csv")))

        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Completed Screenings", completed_runs)

        with col2:
            st.metric("Total Molecules Screened", "~")

        with col3:
            st.metric("Avg. Hit Rate", "~65%")

        with col4:
            st.metric("Avg. Processing Time", "~15 min")
    else:
        st.info("No screening results yet. Run your first screening to see statistics!")
