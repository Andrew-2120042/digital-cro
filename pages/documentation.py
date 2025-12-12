"""
Documentation page
"""

import streamlit as st


def show():
    """Display documentation"""

    st.markdown('<p class="main-header">📖 Documentation</p>', unsafe_allow_html=True)

    st.markdown("---")

    # Quick start
    with st.expander("🚀 Quick Start Guide", expanded=True):
        st.markdown("""
        ### Getting Started with Digital CRO Platform

        #### Step 1: Prepare Your Data

        **Target Protein:**
        - **Option A:** Use PDB ID (e.g., `1HSG` for HIV-1 Protease)
        - **Option B:** Upload your own PDB file

        **Molecule Library:**
        - Create a CSV file with a `smiles` column
        - Optionally include an `id` column for molecule names
        - Example:

        ```csv
        id,smiles
        aspirin,CC(=O)Oc1ccccc1C(=O)O
        ibuprofen,CC(C)Cc1ccc(cc1)C(C)C(=O)O
        caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
        ```

        #### Step 2: Run Screening

        1. Go to **🔬 New Screening**
        2. Enter project details
        3. Input protein (PDB ID or file)
        4. Upload molecule library CSV
        5. Configure parameters:
           - **Affinity Threshold:** Default -7.0 kcal/mol
           - **Parallel Workers:** Use 4-8 for best performance
        6. Click **"Run Screening"**

        #### Step 3: View Results

        - Monitor progress in real-time
        - View summary statistics
        - Download PDF report and data files
        - Browse past screenings in **📊 Results Browser**
        """)

    # Understanding results
    with st.expander("📊 Understanding Results"):
        st.markdown("""
        ### Interpreting Your Results

        #### Binding Affinity

        Measured in **kcal/mol** (more negative = stronger binding):

        | Range | Interpretation |
        |-------|---------------|
        | < -9.0 | **Excellent** - Very strong binder |
        | -8.0 to -9.0 | **Good** - Strong binder |
        | -7.0 to -8.0 | **Moderate** - Decent hit |
        | -6.0 to -7.0 | **Weak** |
        | > -6.0 | **Very weak** |

        #### Drug-Likeness (QED)

        **Quantitative Estimate of Drug-likeness:**
        - **Score:** 0.0 to 1.0 (higher = more drug-like)
        - **Threshold:** Typically ≥ 0.5 is considered drug-like
        - Takes into account: MW, LogP, HBD, HBA, PSA, rotatable bonds, aromatic rings

        #### Lipinski's Rule of Five

        A molecule **passes** if it has:
        - ✅ Molecular weight ≤ 500 Da
        - ✅ LogP ≤ 5
        - ✅ H-bond donors ≤ 5
        - ✅ H-bond acceptors ≤ 10

        **Molecules that pass Lipinski's Rule are more likely to be orally bioavailable.**

        #### ADMET Properties

        - **BBB Penetration:** Can the molecule cross the blood-brain barrier?
          - Important for CNS drugs
          - Based on MW, LogP, PSA, and other factors

        - **Oral Bioavailability:** Can it be taken as a pill?
          - Combines Lipinski + Veber rules
          - Considers rotatable bonds and TPSA

        - **Synthetic Accessibility (SA Score):** How easy to synthesize?
          - Scale: 1-10 (lower = easier)
          - < 3: Easy
          - 3-6: Moderate
          - > 6: Difficult
        """)

    # Technical details
    with st.expander("🔬 Technical Details"):
        st.markdown("""
        ### Platform Capabilities

        #### Molecular Docking

        - **Software:** AutoDock Vina v1.2.7
        - **Pocket Detection:** Ligand-based method using reference ligands
        - **Exhaustiveness:** 8 (search accuracy parameter)
        - **Binding modes:** 9 conformations per molecule
        - **Grid size:** Automatically calculated based on pocket

        #### ADMET Predictions

        **Implemented Rules:**
        - Lipinski's Rule of Five (drug-likeness)
        - QED (Quantitative Estimate of Drug-likeness)
        - BBB penetration (5-criteria model)
        - Oral bioavailability (Lipinski + Veber)
        - Synthetic Accessibility Score (SA Score)

        **Molecular Descriptors:**
        - Molecular weight, LogP, TPSA
        - H-bond donors/acceptors
        - Rotatable bonds
        - Aromatic rings
        - Binding efficiency (LE)

        #### Performance

        **Throughput:**
        - Single molecule: ~60-90 seconds
        - 100 molecules: ~10-15 minutes
        - 1,000 molecules: ~2-3 hours

        **Scaling:**
        - Parallel processing: Up to 16 workers
        - Linear scaling with CPU cores
        - Recommended: 4-8 workers for balanced performance

        #### Output Files

        **Generated Deliverables:**
        - 📄 **PDF Report:** Professional client-ready document
        - 📊 **CSV Data:** Complete results with all properties
        - 🖼️ **Visualizations:** PNG images (molecular structures, ADMET plots)
        - 🧬 **PDBQT Files:** Docked protein-ligand structures
        - 📦 **ZIP Archive:** Complete results package
        """)

    # File formats
    with st.expander("📁 File Formats"):
        st.markdown("""
        ### Input File Formats

        #### Protein Structure

        **PDB Format (.pdb):**
        - Standard Protein Data Bank format
        - Can be downloaded from RCSB PDB (https://www.rcsb.org)
        - Should contain protein coordinates
        - Can include ligands and water molecules (will be removed)

        **Example PDB IDs:**
        - `1HSG` - HIV-1 Protease
        - `4HVP` - HIV-1 Protease with inhibitor
        - `1E66` - Estrogen receptor
        - `3ERT` - Estrogen receptor with ligand

        #### Molecule Library

        **CSV Format (.csv):**

        Required columns:
        - `smiles`: SMILES string representation

        Optional columns:
        - `id`: Molecule identifier/name

        **Example:**
        ```csv
        id,smiles
        aspirin,CC(=O)Oc1ccccc1C(=O)O
        ibuprofen,CC(C)Cc1ccc(cc1)C(C)C(=O)O
        ```

        **SMILES Format (.smi):**

        Space or tab-separated:
        ```
        CC(=O)Oc1ccccc1C(=O)O aspirin
        CC(C)Cc1ccc(cc1)C(C)C(=O)O ibuprofen
        ```

        ### Output File Formats

        - **PDF:** Report with visualizations and analysis
        - **CSV:** Tabular results (compatible with Excel, Python, R)
        - **PDBQT:** Docked structures (viewable in PyMOL, Chimera)
        - **PNG:** High-resolution figures (300 DPI)
        """)

    # FAQ
    with st.expander("❓ Frequently Asked Questions"):
        st.markdown("""
        ### Common Questions

        #### General

        **Q: What format should my molecule library be in?**

        A: CSV format with a `smiles` column. Optionally include an `id` column for molecule names.

        ---

        **Q: How long does a screening take?**

        A: Depends on library size:
        - 20 molecules: ~1-2 minutes
        - 100 molecules: ~10-15 minutes
        - 1,000 molecules: ~2-3 hours

        ---

        **Q: What's a good binding affinity threshold?**

        A: Standard is **-7.0 kcal/mol**. For stricter filtering, use -8.0 or -9.0.

        ---

        **Q: Can I use my own protein structure?**

        A: Yes! Upload a PDB file instead of using a PDB ID.

        ---

        #### Troubleshooting

        **Q: What if my screening fails?**

        A: Check:
        - ✅ Valid SMILES strings in your library
        - ✅ Protein file is valid PDB format
        - ✅ Sufficient disk space
        - ✅ AutoDock Vina is properly installed

        ---

        **Q: Why are some molecules failing?**

        A: Common reasons:
        - Invalid SMILES syntax
        - Molecules too large (>50 rotatable bonds)
        - Incompatible functional groups
        - Very flexible molecules

        ---

        **Q: How do I interpret the PDF report?**

        A: The report includes:
        - **Executive Summary:** Key findings and statistics
        - **Top Hits Table:** Best binding molecules
        - **ADMET Analysis:** Drug-likeness predictions
        - **Visualizations:** Structures and property plots
        - **Methodology:** Technical details

        ---

        #### Advanced

        **Q: Can I customize the docking parameters?**

        A: Currently the parameters are optimized for general use. Contact support for custom configurations.

        ---

        **Q: How accurate are the predictions?**

        A: Molecular docking provides **estimates** of binding affinity:
        - Good for ranking molecules
        - Experimental validation required
        - Typical correlation: R² ~ 0.5-0.7

        ---

        **Q: Can I download the docked structures?**

        A: Yes! The complete ZIP file includes all PDBQT structures in the `docking/results/` folder.
        """)

    # Best practices
    with st.expander("✨ Best Practices"):
        st.markdown("""
        ### Tips for Successful Screening

        #### Library Preparation

        1. **Clean your SMILES:**
           - Remove salts and counterions
           - Standardize tautomers
           - Check for invalid characters

        2. **Size your library:**
           - Start small (20-100 molecules) for testing
           - Scale up after validation
           - Consider diversity when selecting molecules

        3. **Include controls:**
           - Known active compounds
           - Known inactive compounds
           - Helps validate the docking

        #### Protein Preparation

        1. **Choose the right structure:**
           - High resolution (< 2.5 Å) if possible
           - With bound ligand for pocket detection
           - Remove mutations if unwanted

        2. **Consider multiple structures:**
           - Proteins are flexible
           - Different conformations may bind different ligands
           - Can run multiple screenings

        #### Parameter Selection

        1. **Affinity threshold:**
           - -7.0 kcal/mol: Standard threshold
           - -8.0 kcal/mol: More stringent
           - Adjust based on your target

        2. **Workers:**
           - Use 4-8 for optimal performance
           - More workers = faster, but diminishing returns
           - Don't exceed your CPU cores

        #### Result Interpretation

        1. **Don't rely on a single metric:**
           - Consider binding affinity + ADMET
           - Balance potency with drug-likeness
           - Check for PAINS (pan-assay interference)

        2. **Validate experimentally:**
           - Docking is a prediction tool
           - Always validate hits in the lab
           - Consider orthogonal assays

        3. **Visualize structures:**
           - Download PDBQT files
           - Inspect binding poses
           - Look for key interactions
        """)

    # Contact
    st.markdown("---")
    st.subheader("💬 Need Help?")

    col1, col2 = st.columns(2)

    with col1:
        st.info("""
        **Contact Support:**
        - Email: support@digitalcro.com
        - Response time: 24-48 hours
        """)

    with col2:
        st.info("""
        **Resources:**
        - Documentation: https://docs.digitalcro.com
        - GitHub: https://github.com/digitalcro/platform
        - Video tutorials: Coming soon
        """)
