"""
Results browser page
"""

import streamlit as st
from pathlib import Path
import pandas as pd
from datetime import datetime
import zipfile
import tempfile


def show():
    """Display results browser"""

    st.markdown('<p class="main-header">📊 Results Browser</p>', unsafe_allow_html=True)

    st.markdown("---")

    # Find all result directories
    results_dir = Path("data/outputs")

    if not results_dir.exists():
        st.info("No screening results found. Run your first screening to see results here!")
        return

    # Get all directories with results
    result_folders = []

    for folder in results_dir.iterdir():
        if folder.is_dir():
            csv_path = folder / "final_results.csv"
            pdf_path = list(folder.glob("*_report.pdf"))

            if csv_path.exists():
                result_folders.append({
                    'path': folder,
                    'name': folder.name,
                    'csv': csv_path,
                    'pdf': pdf_path[0] if pdf_path else None,
                    'modified': datetime.fromtimestamp(csv_path.stat().st_mtime)
                })

    if not result_folders:
        st.info("No screening results found.")
        return

    # Sort by date (newest first)
    result_folders.sort(key=lambda x: x['modified'], reverse=True)

    st.write(f"**Found {len(result_folders)} screening result(s)**")

    # Folder selector
    selected_folder = st.selectbox(
        "Select Screening Campaign",
        options=result_folders,
        format_func=lambda x: f"{x['name']} - {x['modified'].strftime('%Y-%m-%d %H:%M')}"
    )

    if selected_folder:
        display_result_details(selected_folder)


def display_result_details(result_info):
    """Display detailed results for selected campaign"""

    st.markdown("---")

    # Load data
    df = pd.read_csv(result_info['csv'])

    # Summary metrics
    st.subheader("📈 Campaign Summary")

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Total Molecules", len(df))

    with col2:
        hits = df[df['binding_affinity'] <= -7.0]
        hit_rate = (len(hits) / len(df)) * 100 if len(df) > 0 else 0
        st.metric("Hit Rate", f"{hit_rate:.1f}%", f"{len(hits)} hits")

    with col3:
        best_affinity = df['binding_affinity'].min()
        best_mol = df.loc[df['binding_affinity'].idxmin(), 'ligand_id']
        st.metric("Best Affinity", f"{best_affinity:.2f} kcal/mol", f"{best_mol}")

    with col4:
        if 'qed_score' in df.columns:
            mean_qed = df['qed_score'].mean()
            st.metric("Mean QED", f"{mean_qed:.3f}")
        else:
            st.metric("Mean QED", "N/A")

    # Tabs for different views
    tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
        "📋 All Results",
        "🏆 Top Hits",
        "🔬 3D Viewer",
        "📊 Visualizations",
        "💊 Pharmacophore",
        "🧬 SAR Analysis",
        "📥 Downloads"
    ])

    with tab1:
        display_all_results(df)

    with tab2:
        display_top_hits(df)

    with tab3:
        display_3d_viewer_tab(result_info['path'], df)

    with tab4:
        display_visualizations(result_info['path'])

    with tab5:
        display_pharmacophore_tab(result_info['path'])

    with tab6:
        display_sar_tab(result_info['path'])

    with tab7:
        display_downloads(result_info)


def display_all_results(df):
    """Display all results with filtering"""

    st.subheader("📋 All Results")

    # Add filters
    col1, col2, col3 = st.columns(3)

    with col1:
        if 'lipinski_compliant' in df.columns:
            filter_lipinski = st.checkbox("Only Lipinski Compliant", value=False)
        else:
            filter_lipinski = False

    with col2:
        min_affinity = st.number_input(
            "Min Affinity (kcal/mol)",
            value=-12.0,
            max_value=0.0,
            step=0.5
        )

    with col3:
        if 'qed_score' in df.columns:
            min_qed = st.slider("Min QED Score", 0.0, 1.0, 0.0, 0.05)
        else:
            min_qed = 0.0

    # Apply filters
    filtered_df = df.copy()

    if filter_lipinski and 'lipinski_compliant' in df.columns:
        filtered_df = filtered_df[filtered_df['lipinski_compliant'] == True]

    filtered_df = filtered_df[filtered_df['binding_affinity'] >= min_affinity]

    if min_qed > 0 and 'qed_score' in df.columns:
        filtered_df = filtered_df[filtered_df['qed_score'] >= min_qed]

    st.write(f"**Showing {len(filtered_df)} of {len(df)} molecules**")

    # Sort by binding affinity
    filtered_df = filtered_df.sort_values('binding_affinity')

    # Display table
    st.dataframe(filtered_df, use_container_width=True, height=400)


def display_top_hits(df):
    """Display top hits with filtering options"""

    st.subheader("🏆 Top Hits")

    col1, col2 = st.columns(2)

    with col1:
        threshold = st.slider(
            "Affinity Threshold (kcal/mol)",
            min_value=-12.0,
            max_value=-4.0,
            value=-7.0,
            step=0.1
        )

    with col2:
        top_n = st.number_input(
            "Number of hits to show",
            min_value=5,
            max_value=100,
            value=20,
            step=5
        )

    # Filter and sort
    hits = df[df['binding_affinity'] <= threshold]
    hits = hits.nsmallest(top_n, 'binding_affinity')

    st.write(f"**Showing {len(hits)} molecules with affinity ≤ {threshold} kcal/mol**")

    # Create display dataframe
    display_cols = ['ligand_id', 'binding_affinity']

    # Add optional columns if they exist
    optional_cols = ['qed_score', 'molecular_weight', 'logp', 'lipinski_compliant',
                    'oral_bioavailability', 'bbb_penetrant']

    for col in optional_cols:
        if col in hits.columns:
            display_cols.append(col)

    # Rename columns for display
    display_df = hits[display_cols].copy()

    column_names = {
        'ligand_id': 'Molecule',
        'binding_affinity': 'Affinity (kcal/mol)',
        'qed_score': 'QED',
        'molecular_weight': 'MW (Da)',
        'logp': 'LogP',
        'lipinski_compliant': 'Lipinski Pass',
        'oral_bioavailability': 'Oral Bioavailability',
        'bbb_penetrant': 'BBB Penetrant'
    }

    display_df = display_df.rename(columns=column_names)

    # Display with styling
    st.dataframe(
        display_df,
        use_container_width=True,
        hide_index=True,
        height=400
    )

    # Distribution charts
    st.markdown("### 📊 Property Distributions")

    col1, col2 = st.columns(2)

    with col1:
        # Affinity distribution
        st.bar_chart(
            df['binding_affinity'].value_counts().sort_index(),
            use_container_width=True
        )
        st.caption("Binding Affinity Distribution")

    with col2:
        if 'qed_score' in df.columns:
            st.bar_chart(
                df['qed_score'].value_counts().sort_index(),
                use_container_width=True
            )
            st.caption("QED Score Distribution")


def display_visualizations(result_path):
    """Display visualization images"""

    st.subheader("📊 Visualizations")

    viz_dir = result_path / "visualizations"

    if not viz_dir.exists():
        st.info("No visualizations found for this campaign.")
        return

    # Find all PNG files
    viz_files = list(viz_dir.glob("*.png"))

    if not viz_files:
        st.info("No visualization images found.")
        return

    # Display images
    for viz_file in sorted(viz_files):
        st.markdown(f"**{viz_file.stem.replace('_', ' ').title()}**")
        st.image(str(viz_file), use_column_width=True)
        st.markdown("---")


def display_downloads(result_info):
    """Display download options"""

    st.subheader("📥 Download Files")

    col1, col2 = st.columns(2)

    with col1:
        # PDF report
        if result_info['pdf'] and result_info['pdf'].exists():
            with open(result_info['pdf'], 'rb') as f:
                st.download_button(
                    "📄 Download PDF Report",
                    data=f,
                    file_name=result_info['pdf'].name,
                    mime="application/pdf",
                    use_container_width=True
                )
        else:
            st.button("📄 PDF Not Available", disabled=True, use_container_width=True)

        # CSV results
        if result_info['csv'].exists():
            with open(result_info['csv'], 'rb') as f:
                st.download_button(
                    "📊 Download CSV Results",
                    data=f,
                    file_name=result_info['csv'].name,
                    mime="text/csv",
                    use_container_width=True
                )
        else:
            st.button("📊 CSV Not Available", disabled=True, use_container_width=True)

    with col2:
        # ZIP all files
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix='.zip') as tmp:
                with zipfile.ZipFile(tmp.name, 'w', zipfile.ZIP_DEFLATED) as zipf:
                    for file in result_info['path'].rglob('*'):
                        if file.is_file():
                            zipf.write(
                                file,
                                file.relative_to(result_info['path'])
                            )

                with open(tmp.name, 'rb') as f:
                    st.download_button(
                        "📦 Download Complete Results (ZIP)",
                        data=f,
                        file_name=f"{result_info['name']}_complete.zip",
                        mime="application/zip",
                        use_container_width=True
                    )
        except Exception as e:
            st.button(f"📦 ZIP Failed: {str(e)}", disabled=True, use_container_width=True)

        # View log files
        log_files = list(result_info['path'].glob("*.log"))
        if log_files:
            if st.button("📝 View Logs", use_container_width=True, key="view_logs_btn"):
                for log_file in log_files:
                    with open(log_file, 'r') as f:
                        st.text_area(
                            f"Log: {log_file.name}",
                            value=f.read(),
                            height=200
                        )

    # Professional visualization export
    st.markdown("---")
    st.markdown("### 🔬 Export for Professional Software")

    st.info("""
Export docked structures for visualization in:
- **PyMOL** - Publication-quality rendering
- **Chimera** - Structural analysis
- **VMD** - Molecular dynamics visualization
    """)

    col1, col2, col3 = st.columns(3)

    with col1:
        if st.button("📦 Export to PyMOL", use_container_width=True, key="export_pymol_btn"):
            export_to_pymol(result_info['path'])

    with col2:
        if st.button("📦 Export to Chimera", use_container_width=True, key="export_chimera_btn"):
            export_to_chimera(result_info['path'])

    with col3:
        if st.button("📦 Export to VMD", use_container_width=True, key="export_vmd_btn"):
            export_to_vmd(result_info['path'])


def display_3d_viewer_tab(result_path, df):
    """Display 3D protein-ligand complex viewer"""

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from utils.mol_viewer import create_complex_viewer, create_multi_pose_viewer
    import streamlit.components.v1 as components
    import json

    st.subheader("🔬 3D Protein-Ligand Complex")

    # Find protein and ligand files
    protein_dir = result_path / 'proteins'
    docking_dir = result_path / 'docking' / 'results'

    # Get protein file
    pdb_files = list(protein_dir.glob('*_cleaned.pdb'))
    if not pdb_files:
        pdb_files = list(protein_dir.glob('*.pdb'))

    if not pdb_files:
        st.warning("Protein structure not found in results.")
        return

    protein_pdb = pdb_files[0]

    # Get docked results
    if not docking_dir.exists():
        st.warning("No docking results found.")
        return

    ligand_files = sorted(docking_dir.glob('*_docked.pdbqt'))

    if not ligand_files:
        st.warning("No docked ligand structures found.")
        return

    # Molecule selector
    col1, col2 = st.columns([2, 1])

    with col1:
        selected_molecule = st.selectbox(
            "Select Molecule to View",
            options=df['ligand_id'].tolist(),
            format_func=lambda x: f"{x} ({df[df['ligand_id']==x]['binding_affinity'].values[0]:.2f} kcal/mol)"
        )

    with col2:
        view_mode = st.radio(
            "View Mode",
            ["Single Pose", "All Poses"],
            horizontal=True
        )

    # Find ligand file(s) for selected molecule
    molecule_files = [f for f in ligand_files if selected_molecule in f.stem]

    if not molecule_files:
        st.warning(f"Structure file not found for {selected_molecule}")
        return

    # Get molecule info
    mol_info = df[df['ligand_id'] == selected_molecule].iloc[0]

    # Display molecule info
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Binding Affinity", f"{mol_info['binding_affinity']:.2f} kcal/mol")

    with col2:
        if 'qed_score' in mol_info:
            st.metric("QED Score", f"{mol_info['qed_score']:.2f}")

    with col3:
        if 'molecular_weight' in mol_info:
            st.metric("Molecular Weight", f"{mol_info['molecular_weight']:.1f} Da")

    with col4:
        if 'lipinski_compliant' in mol_info:
            status = "✅ Pass" if mol_info['lipinski_compliant'] else "❌ Fail"
            st.metric("Lipinski", status)

    # Read pocket info
    pocket_center = None
    pocket_size = None

    docking_json = result_path / 'docking' / 'docking_results_detailed.json'
    if docking_json.exists():
        try:
            with open(docking_json, 'r') as f:
                docking_data = json.load(f)
                if 'pocket_center' in docking_data:
                    pocket_center = tuple(docking_data['pocket_center'])
                if 'pocket_size' in docking_data:
                    pocket_size = tuple(docking_data['pocket_size'])
        except:
            pass

    # Create viewer
    try:
        if view_mode == "Single Pose":
            # Show best pose
            viewer = create_complex_viewer(
                protein_pdb_path=str(protein_pdb),
                ligand_pdbqt_path=str(molecule_files[0]),
                pocket_center=pocket_center,
                pocket_size=pocket_size,
                width=900,
                height=650,
                show_pocket=True
            )
        else:
            # Show all poses
            viewer = create_multi_pose_viewer(
                protein_pdb_path=str(protein_pdb),
                ligand_pdbqt_paths=[str(f) for f in molecule_files[:9]],
                width=900,
                height=650
            )

        # Render viewer
        components.html(viewer._make_html(), height=700, scrolling=False)

        # Instructions
        st.caption("""
        **💡 Interactive Controls:**
        - 🖱️ **Left-click + drag** to rotate
        - 🖱️ **Right-click + drag** or **scroll** to zoom
        - **Protein** shown in gray (cartoon representation)
        - **Ligand** shown in colored sticks
        - **Binding pocket residues** highlighted in orange
        """)

        if view_mode == "All Poses":
            st.info(f"Showing {len(molecule_files[:9])} binding poses with different colors. Each pose represents a potential binding configuration.")

    except Exception as e:
        st.error(f"Error creating 3D viewer: {e}")
        st.exception(e)


def export_to_pymol(result_path):
    """Export results to PyMOL format"""
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from scripts.utils.export_formats import create_pymol_session_files

    result_path = Path(result_path)

    # Find protein and ligands
    protein_pdb_files = list(result_path.glob("proteins/*_cleaned.pdb"))
    if not protein_pdb_files:
        st.error("Protein structure not found")
        return

    protein_pdb = protein_pdb_files[0]

    docking_results_dir = result_path / "docking" / "results"
    if not docking_results_dir.exists():
        st.error("Docking results not found")
        return

    ligand_pdbqts = sorted(docking_results_dir.glob("*_docked.pdbqt"))[:10]  # Top 10

    if not ligand_pdbqts:
        st.error("No docked ligands found")
        return

    # Get ligand names from results
    results_csv = result_path / "final_results.csv"
    if results_csv.exists():
        df = pd.read_csv(results_csv)
        ligand_names = df.head(10)['ligand_id'].tolist()
    else:
        ligand_names = [f.stem.replace('_docked', '') for f in ligand_pdbqts]

    # Create export directory
    export_dir = result_path / "exports" / "pymol"
    export_dir.mkdir(parents=True, exist_ok=True)

    # Generate scripts
    with st.spinner("Generating PyMOL scripts..."):
        try:
            scripts = create_pymol_session_files(
                protein_pdb=str(protein_pdb),
                ligand_pdbqts=[str(p) for p in ligand_pdbqts],
                output_dir=str(export_dir),
                ligand_names=ligand_names
            )

            # Create ZIP
            zip_path = result_path / "exports" / "pymol_export.zip"

            with zipfile.ZipFile(zip_path, 'w') as zipf:
                # Add scripts
                for script in scripts:
                    zipf.write(script, Path(script).name)

                # Add protein
                zipf.write(protein_pdb, protein_pdb.name)

                # Add ligands
                for ligand in ligand_pdbqts:
                    zipf.write(ligand, ligand.name)

            # Download button
            with open(zip_path, 'rb') as f:
                st.download_button(
                    "📥 Download PyMOL Package",
                    data=f,
                    file_name="pymol_export.zip",
                    mime="application/zip",
                    use_container_width=True,
                    key="download_pymol_zip"
                )

            st.success(f"✓ Created {len(scripts)} PyMOL scripts for top {len(ligand_pdbqts)} hits!")

            with st.expander("📖 How to use in PyMOL"):
                st.code("""
# 1. Unzip the downloaded package
# 2. Open PyMOL
# 3. File → Run Script → Select a .pml file
# Or use PyMOL console:
@all_ligands_pymol.pml

# Individual ligands:
@ligand_name_pymol.pml
                """, language="bash")

        except Exception as e:
            st.error(f"Export failed: {e}")


def export_to_chimera(result_path):
    """Export to Chimera format"""
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from scripts.utils.export_formats import create_chimera_script

    result_path = Path(result_path)

    protein_pdb_files = list(result_path.glob("proteins/*_cleaned.pdb"))
    if not protein_pdb_files:
        st.error("Protein structure not found")
        return

    protein_pdb = protein_pdb_files[0]

    docking_results_dir = result_path / "docking" / "results"
    if not docking_results_dir.exists():
        st.error("Docking results not found")
        return

    ligand_pdbqts = sorted(docking_results_dir.glob("*_docked.pdbqt"))[:10]

    if not ligand_pdbqts:
        st.error("No docked ligands found")
        return

    results_csv = result_path / "final_results.csv"
    if results_csv.exists():
        df = pd.read_csv(results_csv)
        ligand_names = df.head(10)['ligand_id'].tolist()
    else:
        ligand_names = [f.stem.replace('_docked', '') for f in ligand_pdbqts]

    export_dir = result_path / "exports" / "chimera"
    export_dir.mkdir(parents=True, exist_ok=True)

    scripts = []

    with st.spinner("Generating Chimera scripts..."):
        try:
            for ligand_path, name in zip(ligand_pdbqts, ligand_names):
                script = create_chimera_script(
                    protein_pdb=str(protein_pdb),
                    ligand_pdbqt=str(ligand_path),
                    output_path=str(export_dir / f"{name}_chimera.py"),
                    ligand_name=name
                )
                scripts.append(script)

            zip_path = result_path / "exports" / "chimera_export.zip"

            with zipfile.ZipFile(zip_path, 'w') as zipf:
                for script in scripts:
                    zipf.write(script, Path(script).name)
                zipf.write(protein_pdb, protein_pdb.name)
                for ligand in ligand_pdbqts:
                    zipf.write(ligand, ligand.name)

            with open(zip_path, 'rb') as f:
                st.download_button(
                    "📥 Download Chimera Package",
                    data=f,
                    file_name="chimera_export.zip",
                    mime="application/zip",
                    use_container_width=True,
                    key="download_chimera_zip"
                )

            st.success(f"✓ Created {len(scripts)} Chimera scripts for top {len(ligand_pdbqts)} hits!")

            with st.expander("📖 How to use in Chimera"):
                st.code("""
# 1. Unzip the downloaded package
# 2. Open UCSF Chimera
# 3. File → Open → Select a .py script
# Or use Chimera command line:
chimera --script ligand_name_chimera.py
                """, language="bash")

        except Exception as e:
            st.error(f"Export failed: {e}")


def export_to_vmd(result_path):
    """Export to VMD format"""
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from scripts.utils.export_formats import create_vmd_script

    result_path = Path(result_path)

    protein_pdb_files = list(result_path.glob("proteins/*_cleaned.pdb"))
    if not protein_pdb_files:
        st.error("Protein structure not found")
        return

    protein_pdb = protein_pdb_files[0]

    docking_results_dir = result_path / "docking" / "results"
    if not docking_results_dir.exists():
        st.error("Docking results not found")
        return

    ligand_pdbqts = sorted(docking_results_dir.glob("*_docked.pdbqt"))[:10]

    if not ligand_pdbqts:
        st.error("No docked ligands found")
        return

    results_csv = result_path / "final_results.csv"
    if results_csv.exists():
        df = pd.read_csv(results_csv)
        ligand_names = df.head(10)['ligand_id'].tolist()
    else:
        ligand_names = [f.stem.replace('_docked', '') for f in ligand_pdbqts]

    export_dir = result_path / "exports" / "vmd"
    export_dir.mkdir(parents=True, exist_ok=True)

    scripts = []

    with st.spinner("Generating VMD scripts..."):
        try:
            for ligand_path, name in zip(ligand_pdbqts, ligand_names):
                script = create_vmd_script(
                    protein_pdb=str(protein_pdb),
                    ligand_pdbqt=str(ligand_path),
                    output_path=str(export_dir / f"{name}_vmd.tcl"),
                    ligand_name=name
                )
                scripts.append(script)

            zip_path = result_path / "exports" / "vmd_export.zip"

            with zipfile.ZipFile(zip_path, 'w') as zipf:
                for script in scripts:
                    zipf.write(script, Path(script).name)
                zipf.write(protein_pdb, protein_pdb.name)
                for ligand in ligand_pdbqts:
                    zipf.write(ligand, ligand.name)

            with open(zip_path, 'rb') as f:
                st.download_button(
                    "📥 Download VMD Package",
                    data=f,
                    file_name="vmd_export.zip",
                    mime="application/zip",
                    use_container_width=True,
                    key="download_vmd_zip"
                )

            st.success(f"✓ Created {len(scripts)} VMD scripts for top {len(ligand_pdbqts)} hits!")

            with st.expander("📖 How to use in VMD"):
                st.code("""
# 1. Unzip the downloaded package
# 2. Open VMD
# 3. Extensions → Tk Console
# 4. In console, type:
source ligand_name_vmd.tcl

# Or from command line:
vmd -e ligand_name_vmd.tcl
                """, language="bash")

        except Exception as e:
            st.error(f"Export failed: {e}")


def display_pharmacophore_tab(result_path):
    """Display pharmacophore analysis"""

    st.subheader("💊 Pharmacophore Analysis")

    st.info("""
    **Pharmacophore hypothesis** identifies common 3D features across top-binding molecules.
    This helps understand WHY molecules bind and guides future drug design.

    Features analyzed:
    - Hydrogen bond donors/acceptors
    - Aromatic rings
    - Hydrophobic regions
    - Ionizable groups
    """)

    # Check if analysis exists
    result_path = Path(result_path)
    pharmaco_dir = result_path / "pharmacophore"

    if not pharmaco_dir.exists():
        # Generate analysis
        st.markdown("### 🔬 Generate Analysis")

        col1, col2 = st.columns(2)

        with col1:
            top_n = st.number_input(
                "Number of top hits to analyze",
                min_value=5,
                max_value=50,
                value=20,
                step=5
            )

        with col2:
            min_occurrence = st.slider(
                "Minimum occurrence (%)",
                min_value=50,
                max_value=100,
                value=70,
                step=10
            ) / 100

        if st.button("🔬 Generate Pharmacophore Analysis", type="primary", use_container_width=True, key="btn_gen_pharmaco"):
            generate_pharmacophore_analysis(result_path, top_n, min_occurrence)
            st.rerun()
    else:
        # Display existing analysis
        display_pharmacophore_results(pharmaco_dir, result_path)


def generate_pharmacophore_analysis(result_path, top_n=20, min_occurrence=0.7):
    """Generate pharmacophore analysis for results"""
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from scripts.utils.pharmacophore import analyze_pharmacophore_from_results
    from scripts.utils.pharmacophore_viz import (
        plot_feature_distribution,
        plot_feature_importance,
        plot_pharmacophore_summary
    )

    result_path = Path(result_path)
    results_csv = result_path / "final_results.csv"
    pharmaco_dir = result_path / "pharmacophore"
    pharmaco_dir.mkdir(exist_ok=True)

    with st.spinner("Analyzing pharmacophore features..."):
        try:
            # Run analysis
            results = analyze_pharmacophore_from_results(
                results_csv=str(results_csv),
                top_n=top_n,
                min_occurrence=min_occurrence,
                output_dir=str(pharmaco_dir)
            )

            # Generate visualizations
            plot_feature_distribution(
                results['feature_df'],
                str(pharmaco_dir / "feature_distribution.png")
            )

            plot_feature_importance(
                results['importance'],
                str(pharmaco_dir / "feature_importance.png")
            )

            plot_pharmacophore_summary(
                results['feature_df'],
                results['hypothesis'],
                str(pharmaco_dir / "pharmacophore_summary.png")
            )

            st.success("✓ Pharmacophore analysis complete!")

        except Exception as e:
            st.error(f"Analysis failed: {e}")
            import traceback
            st.exception(e)


def display_pharmacophore_results(pharmaco_dir, result_path):
    """Display pharmacophore analysis results"""

    # Hypothesis
    hypothesis_file = pharmaco_dir / "pharmacophore_hypothesis.txt"

    if hypothesis_file.exists():
        with open(hypothesis_file, 'r') as f:
            hypothesis_text = f.read()

        st.markdown("### 🎯 Pharmacophore Hypothesis")
        st.success(hypothesis_text)

    st.markdown("---")

    # Visualizations
    st.markdown("### 📊 Feature Analysis")

    # Summary
    summary_img = pharmaco_dir / "pharmacophore_summary.png"
    if summary_img.exists():
        st.image(str(summary_img), use_column_width=True)
        st.markdown("---")

    # Tabs for detailed views
    detail_tab1, detail_tab2 = st.tabs(["Feature Distribution", "Feature Importance"])

    with detail_tab1:
        dist_img = pharmaco_dir / "feature_distribution.png"
        if dist_img.exists():
            st.image(str(dist_img), use_column_width=True)
            st.caption("Distribution of pharmacophore features across analyzed molecules")

    with detail_tab2:
        imp_img = pharmaco_dir / "feature_importance.png"
        if imp_img.exists():
            st.image(str(imp_img), use_column_width=True)
            st.caption("Correlation between features and binding affinity (negative = better binding)")

    # Download data
    st.markdown("---")
    st.markdown("### 📥 Download Analysis Data")

    col1, col2, col3 = st.columns(3)

    with col1:
        features_csv = pharmaco_dir / "pharmacophore_features.csv"
        if features_csv.exists():
            with open(features_csv, 'rb') as f:
                st.download_button(
                    "📊 Feature Data (CSV)",
                    data=f,
                    file_name="pharmacophore_features.csv",
                    mime="text/csv",
                    use_container_width=True,
                    key="download_pharma_features"
                )

    with col2:
        importance_csv = pharmaco_dir / "feature_importance.csv"
        if importance_csv.exists():
            with open(importance_csv, 'rb') as f:
                st.download_button(
                    "📈 Feature Importance (CSV)",
                    data=f,
                    file_name="feature_importance.csv",
                    mime="text/csv",
                    use_container_width=True,
                    key="download_pharma_importance"
                )

    with col3:
        if hypothesis_file.exists():
            with open(hypothesis_file, 'rb') as f:
                st.download_button(
                    "📄 Hypothesis (TXT)",
                    data=f,
                    file_name="pharmacophore_hypothesis.txt",
                    mime="text/plain",
                    use_container_width=True,
                    key="download_pharma_hypothesis"
                )

    # Regenerate option
    st.markdown("---")
    if st.button("🔄 Regenerate Analysis", use_container_width=False):
        import shutil
        shutil.rmtree(pharmaco_dir)
        st.success("Analysis cleared! Click 'Generate' to create new analysis.")
        st.rerun()


def display_sar_tab(result_path):
    """Display SAR analysis"""

    st.subheader("🧬 Structure-Activity Relationship (SAR) Analysis")

    st.info("""
    **SAR Analysis** identifies molecular scaffolds and analyzes how structural changes affect binding.
    This guides which analogs to synthesize next for lead optimization.

    **Features:**
    - Scaffold identification and analysis
    - Matched molecular pair (MMP) detection
    - Property-activity correlations
    - Synthesis recommendations
    """)

    # Check if analysis exists
    result_path = Path(result_path)
    sar_dir = result_path / "sar_analysis"

    if not sar_dir.exists():
        # Configuration
        st.markdown("### 🔬 Generate Analysis")

        col1, col2 = st.columns(2)

        with col1:
            min_scaffold_size = st.number_input(
                "Minimum molecules per scaffold",
                min_value=2,
                max_value=10,
                value=2,
                help="Ignore scaffolds with fewer molecules"
            )

        # Generate button
        if st.button("🔬 Generate SAR Analysis", type="primary", use_container_width=True, key="btn_gen_sar"):
            generate_sar_analysis(result_path, min_scaffold_size)
            st.rerun()
    else:
        # Display existing analysis
        display_sar_results(sar_dir)

        # Regenerate option
        st.markdown("---")
        if st.button("🔄 Regenerate Analysis", use_container_width=False, key="btn_regen_sar"):
            import shutil
            shutil.rmtree(sar_dir)
            st.success("Analysis cleared! Generate new analysis above.")
            st.rerun()


def generate_sar_analysis(result_path, min_scaffold_size):
    """Generate SAR analysis"""
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from scripts.utils.sar_analysis import analyze_sar_from_results
    from scripts.utils.sar_viz import (
        plot_scaffold_comparison,
        plot_mmp_network,
        plot_property_correlations,
        plot_sar_summary
    )

    result_path = Path(result_path)
    results_csv = result_path / "final_results.csv"
    sar_dir = result_path / "sar_analysis"
    sar_dir.mkdir(exist_ok=True)

    with st.spinner("Analyzing structure-activity relationships..."):
        try:
            # Run analysis
            results = analyze_sar_from_results(
                results_csv=str(results_csv),
                min_molecules_per_scaffold=min_scaffold_size,
                output_dir=str(sar_dir)
            )

            # Generate visualizations
            if results['scaffolds']:
                plot_scaffold_comparison(
                    results['scaffolds'],
                    str(sar_dir / "scaffold_comparison.png")
                )

            if results['mmps']:
                plot_mmp_network(
                    results['mmps'],
                    str(sar_dir / "matched_pairs.png")
                )

            if not results['correlations'].empty:
                plot_property_correlations(
                    results['correlations'],
                    str(sar_dir / "property_correlations.png")
                )

            # Summary plot
            plot_sar_summary(
                results['scaffolds'],
                results['mmps'],
                results['correlations'],
                str(sar_dir / "sar_summary.png")
            )

            st.success("✓ SAR analysis complete!")

        except Exception as e:
            st.error(f"Analysis failed: {e}")
            import traceback
            st.exception(e)


def display_sar_results(sar_dir):
    """Display SAR analysis results"""

    # Summary visualization
    st.markdown("### 📊 SAR Summary")

    summary_img = sar_dir / "sar_summary.png"
    if summary_img.exists():
        st.image(str(summary_img), use_column_width=True)

    st.markdown("---")

    # Recommendations
    rec_file = sar_dir / "sar_recommendations.txt"
    if rec_file.exists():
        with open(rec_file, 'r') as f:
            recommendations = f.read()

        st.markdown("### 🎯 Optimization Recommendations")
        st.success(recommendations)

    st.markdown("---")

    # Detailed tabs
    detail_tab1, detail_tab2, detail_tab3 = st.tabs([
        "Scaffold Analysis",
        "Matched Pairs",
        "Property Correlations"
    ])

    with detail_tab1:
        st.markdown("#### Molecular Scaffolds")

        st.info("Core structures (Murcko scaffolds) shared by multiple molecules")

        scaffold_img = sar_dir / "scaffold_comparison.png"
        if scaffold_img.exists():
            st.image(str(scaffold_img), use_column_width=True)

        scaffold_csv = sar_dir / "scaffold_summary.csv"
        if scaffold_csv.exists():
            df = pd.read_csv(scaffold_csv)
            st.dataframe(df, use_container_width=True, hide_index=True)

    with detail_tab2:
        st.markdown("#### Matched Molecular Pairs")

        st.info("Pairs of molecules differing by single structural change")

        mmp_img = sar_dir / "matched_pairs.png"
        if mmp_img.exists():
            st.image(str(mmp_img), use_column_width=True)

        mmp_csv = sar_dir / "matched_pairs.csv"
        if mmp_csv.exists():
            df = pd.read_csv(mmp_csv)
            st.dataframe(df.head(20), use_container_width=True, hide_index=True)
            st.caption(f"Showing top 20 of {len(df)} matched pairs")

    with detail_tab3:
        st.markdown("#### Property-Activity Correlations")

        st.info("Correlation between molecular properties and binding affinity")

        corr_img = sar_dir / "property_correlations.png"
        if corr_img.exists():
            st.image(str(corr_img), use_column_width=True)

        corr_csv = sar_dir / "property_correlations.csv"
        if corr_csv.exists():
            df = pd.read_csv(corr_csv)
            st.dataframe(df, use_container_width=True, hide_index=True)

    # Downloads
    st.markdown("---")
    st.markdown("### 📥 Download SAR Data")

    col1, col2, col3 = st.columns(3)

    files = [
        ('scaffold_summary.csv', 'Scaffold Data'),
        ('matched_pairs.csv', 'Matched Pairs'),
        ('property_correlations.csv', 'Correlations')
    ]

    for idx, (filename, label) in enumerate(files):
        file_path = sar_dir / filename

        if file_path.exists():
            with [col1, col2, col3][idx]:
                with open(file_path, 'rb') as f:
                    st.download_button(
                        f"📊 {label}",
                        data=f,
                        file_name=filename,
                        mime="text/csv",
                        use_container_width=True,
                        key=f"download_sar_{filename}"
                    )
