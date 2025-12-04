"""
Screening workflow page
"""

import streamlit as st
from pathlib import Path
import pandas as pd
import tempfile
import shutil
from datetime import datetime
import sys
import subprocess
import json

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def show():
    """Display screening workflow page"""

    st.markdown('<p class="main-header">🔬 New Virtual Screening</p>', unsafe_allow_html=True)

    st.markdown("---")

    # Import library manager
    from scripts.utils.library_manager import LibraryManager

    # ========== SCREENING CONFIGURATION FORM ==========
    with st.form("screening_config"):
        st.subheader("📋 Project Configuration")

        col1, col2 = st.columns(2)

        with col1:
            project_name = st.text_input(
                "Project Name",
                value="Drug Discovery Project",
                help="Name for this screening campaign"
            )

            client_name = st.text_input(
                "Client/Organization Name",
                value="Research Lab",
                help="Your organization name for the report"
            )

        with col2:
            affinity_threshold = st.number_input(
                "Binding Affinity Threshold (kcal/mol)",
                value=-7.0,
                step=0.1,
                help="Molecules with affinity ≤ this value are considered hits"
            )

            max_workers = st.number_input(
                "Parallel Workers",
                value=4,
                min_value=1,
                max_value=16,
                help="Number of CPU cores to use for docking"
            )

        st.markdown("---")

        # Consensus docking option
        st.subheader("🎯 Docking Method")

        use_consensus = st.checkbox(
            "Use Consensus Docking (Multi-Method Validation)",
            value=False,
            help="Run multiple docking programs and combine results for higher confidence. Takes 2-3x longer but provides validated results."
        )

        if use_consensus:
            st.info("""
            **Consensus Mode:** Your molecules will be docked using multiple programs:
            - AutoDock Vina (default)
            - Smina (improved scoring)

            Only molecules with high agreement between methods will be flagged as high-confidence hits.

            ⏱️ Runtime: ~2-3x longer than standard docking
            💰 Value: Industry-standard validation, justifies premium pricing
            """)

        st.markdown("---")

        # Protein input section
        st.subheader("🧬 Target Protein")

        protein_input_method = st.radio(
            "Protein Input Method",
            ["PDB ID", "Upload PDB File"],
            horizontal=True
        )

        if protein_input_method == "PDB ID":
            pdb_id = st.text_input(
                "PDB ID",
                value="1HSG",
                max_chars=4,
                help="4-letter PDB code (e.g., 1HSG for HIV-1 Protease)"
            ).upper()

            protein_file = None

        else:
            pdb_id = None
            protein_file = st.file_uploader(
                "Upload PDB File",
                type=['pdb'],
                help="Upload protein structure in PDB format"
            )

        st.markdown("---")

        # ========== LIBRARY SECTION - STAYS INSIDE FORM ==========
        st.subheader("💊 Molecule Library")

        # Get saved libraries
        manager = LibraryManager()
        saved_libraries = manager.list_libraries()

        # Initialize variables
        library_file = None
        selected_library_id = None

        # Show radio button only if there are saved libraries
        if saved_libraries:
            library_source = st.radio(
                "Library Source",
                ["📚 Use Saved Library", "📤 Upload New File"],
                key="lib_source_radio",
                horizontal=True
            )
        else:
            library_source = "📤 Upload New File"
            st.info("💡 No saved libraries found. You can save libraries in the '📚 My Libraries' page for quick reuse!")

        st.write("")  # spacing

        # Show appropriate widget based on selection
        # These render EVERY time the form loads
        if library_source == "📚 Use Saved Library" and saved_libraries:
            # ========== SAVED LIBRARY OPTION ==========

            # Create display names
            library_options = {}
            for lib in saved_libraries:
                mol_count = lib.get('num_molecules', lib.get('molecule_count', 0))
                display_name = f"{lib['name']} ({mol_count} molecules)"
                library_options[display_name] = lib

            # Check if library was pre-selected from Libraries page
            default_index = 0
            if 'selected_library' in st.session_state:
                for i, lib in enumerate(saved_libraries):
                    if lib.get('library_id') == st.session_state['selected_library']:
                        default_index = i
                        break
                st.session_state.pop('selected_library')  # Clear after use

            selected_display = st.selectbox(
                "Select Library",
                options=list(library_options.keys()),
                index=default_index,
                key="library_dropdown"
            )

            selected_library = library_options[selected_display]
            selected_library_id = selected_library.get('library_id')

            # Show info about selected library
            mol_count = selected_library.get('num_molecules', selected_library.get('molecule_count', 0))
            last_used = selected_library.get('last_used', 'Never')
            if last_used != 'Never' and len(last_used) > 10:
                last_used = last_used[:10]

            st.success(f"""
            ✅ **Using:** {selected_library['name']}
            - 📊 {mol_count} molecules
            - 📝 {selected_library.get('description', 'No description')}
            - 📅 Last used: {last_used}
            """)

            library_file = None

        else:
            # ========== UPLOAD NEW FILE OPTION ==========
            # This shows when "Upload New File" is selected OR when no libraries exist

            st.markdown("**Upload Molecule Library**")

            library_file = st.file_uploader(
                "Drag and drop file here",
                type=['csv', 'smi', 'txt'],
                help="Supported formats: CSV (id,smiles), SMI (SMILES<tab>ID)",
                key="library_file_uploader"
            )

            if library_file:
                try:
                    # Preview file
                    if library_file.name.endswith('.csv'):
                        import io

                        # Read file
                        file_contents = library_file.read()
                        df = pd.read_csv(io.BytesIO(file_contents), dtype=str)

                        st.success(f"✅ **Loaded:** {len(df)} molecules from {library_file.name}")

                        # Show preview
                        with st.expander("📋 Preview (first 5 rows)"):
                            st.dataframe(df.head(), use_container_width=True)

                        # Reset file pointer for later use
                        library_file.seek(0)

                    else:
                        st.success(f"✅ **Uploaded:** {library_file.name}")

                except Exception as e:
                    st.error(f"Error reading file: {e}")

            selected_library_id = None

        st.markdown("---")

        # Job queue option
        submit_to_queue = st.checkbox(
            "📋 Submit to job queue (run in background)",
            value=False,
            help="Queue this job for background processing instead of running now"
        )

        # Submit button
        submit_label = "📋 Add to Queue" if submit_to_queue else "🚀 Run Screening"
        submitted = st.form_submit_button(
            submit_label,
            use_container_width=True,
            type="primary"
        )

    # Process screening if form submitted
    if submitted:
        # Validation
        errors = []

        if protein_input_method == "PDB ID" and not pdb_id:
            errors.append("❌ Please enter a PDB ID")

        if protein_input_method == "Upload PDB File" and not protein_file:
            errors.append("❌ Please upload a PDB file")

        # Check library - either file upload or saved library
        if library_source == "📚 Use Saved Library" and not selected_library_id:
            errors.append("❌ Please select a saved library")
        elif library_source == "📤 Upload New File" and not library_file:
            errors.append("❌ Please upload a molecule library")

        if errors:
            for error in errors:
                st.error(error)
        else:
            # Get library path
            manager = LibraryManager()

            if library_source == "📚 Use Saved Library":
                # Use saved library
                library_path_to_use = manager.get_library_path(selected_library_id)
                manager.update_usage(selected_library_id)
                library_file_to_use = None  # Not an uploaded file
            else:
                # Use uploaded library
                library_path_to_use = None
                library_file_to_use = library_file

            if submit_to_queue:
                # Submit to queue instead of running immediately
                submit_to_job_queue(
                    pdb_id=pdb_id,
                    protein_file=protein_file,
                    library_file=library_file_to_use,
                    library_path=library_path_to_use,
                    project_name=project_name,
                    client_name=client_name,
                    affinity_threshold=affinity_threshold,
                    max_workers=max_workers,
                    use_consensus=use_consensus
                )
            else:
                # Run screening immediately
                run_screening(
                    pdb_id=pdb_id,
                    protein_file=protein_file,
                    library_file=library_file_to_use,
                    library_path=library_path_to_use,
                    project_name=project_name,
                    client_name=client_name,
                    affinity_threshold=affinity_threshold,
                    max_workers=max_workers,
                    use_consensus=use_consensus
                )


def submit_to_job_queue(
    pdb_id,
    protein_file,
    library_file,
    library_path,
    project_name,
    client_name,
    affinity_threshold,
    max_workers,
    use_consensus=False
):
    """Submit screening job to queue for background processing"""

    from scripts.utils.job_queue import JobQueue

    try:
        queue = JobQueue()

        # Handle library - either uploaded file or saved library path
        if library_path:
            # Using saved library
            library_save_path = library_path
        else:
            # Save uploaded library file permanently to job queue directory
            library_dir = Path("data/job_queue/libraries")
            library_dir.mkdir(parents=True, exist_ok=True)

            # Reset file pointer
            library_file.seek(0)

            # Parse and clean the library file
            if library_file.name.endswith('.csv'):
                import io
                from rdkit import Chem

                # Read CSV with pandas - IGNORE comment lines
                file_contents = library_file.read()
                df = pd.read_csv(io.BytesIO(file_contents), dtype=str, comment='#')

                # Normalize column names
                df.columns = df.columns.str.lower().str.strip()

                # Find id and smiles columns
                id_col = None
                smiles_col = None

                for col in df.columns:
                    if col in ['id', 'name', 'mol_id', 'compound_id', 'molecule_id']:
                        id_col = col
                    elif col in ['smiles', 'smile', 'smi', 'structure']:
                        smiles_col = col

                if not id_col or not smiles_col:
                    st.error(f"❌ CSV must have 'id' and 'smiles' columns. Found: {df.columns.tolist()}")
                    return

                # Keep only id and smiles columns
                df_clean = pd.DataFrame({
                    'id': df[id_col],
                    'smiles': df[smiles_col]
                })

                # Drop rows where SMILES column is NaN
                df_clean = df_clean.dropna(subset=['smiles'])

                # Minimal validation: only remove truly empty SMILES
                # Don't validate with RDKit - let the robust 6-method pipeline handle it
                # Filter out empty entries
                df_clean = df_clean[df_clean['smiles'].notna()]
                df_clean = df_clean[df_clean['smiles'].astype(str).str.strip() != '']
                df_clean = df_clean[~df_clean['smiles'].astype(str).str.lower().isin(['nan', 'none', 'null', ''])]
                df_clean = df_clean.reset_index(drop=True)

                if len(df_clean) == 0:
                    st.error("❌ No valid molecules in library!")
                    return

                # Save as tab-separated file
                smi_path = library_dir / f"{Path(library_file.name).stem}_cleaned.smi"
                with open(smi_path, 'w') as f:
                    for _, row in df_clean.iterrows():
                        f.write(f"{row['smiles']}\t{row['id']}\n")

                library_save_path = smi_path

            elif library_file.name.endswith('.smi'):
                # SMI format - parse and re-save as tab-separated
                import io
                file_contents = library_file.read()
                df = pd.read_csv(io.BytesIO(file_contents), sep=r'\s+', header=None, names=['smiles', 'id'], dtype=str)

                smi_path = library_dir / f"{Path(library_file.name).stem}_cleaned.smi"
                with open(smi_path, 'w') as f:
                    for _, row in df.iterrows():
                        f.write(f"{row['smiles']}\t{row['id']}\n")

                library_save_path = smi_path

            else:
                # Unknown format - save as-is
                library_save_path = library_dir / library_file.name
                with open(library_save_path, 'wb') as f:
                    library_file.seek(0)
                    f.write(library_file.read())

        # Handle protein file upload (if applicable)
        if protein_file:
            # For uploaded PDB files, save them and use the file path
            protein_save_path = Path("data/job_queue/proteins") / protein_file.name
            protein_save_path.parent.mkdir(parents=True, exist_ok=True)

            with open(protein_save_path, 'wb') as f:
                f.write(protein_file.read())

            pdb_id = str(protein_save_path)

        # Submit job to queue
        job_id = queue.submit_single_target_job(
            pdb_id=pdb_id,
            ligand_library_path=str(library_save_path),
            project_name=project_name,
            client_name=client_name,
            affinity_threshold=affinity_threshold,
            max_workers=max_workers,
            use_consensus=use_consensus
        )

        st.success(f"""
        ✅ **Job Submitted to Queue!**

        **Job ID:** `{job_id}`

        Your job has been queued for processing. It will be processed in the order it was received (FIFO).

        **Next Steps:**
        1. Go to the **📋 My Jobs** page to monitor progress
        2. You can submit more jobs while this one processes
        3. Download results when the job completes

        **To process jobs:** Run the job processor in a terminal:
        ```
        python scripts/job_processor.py
        ```

        Or run in continuous daemon mode:
        ```
        python scripts/job_processor.py --daemon
        ```
        """)

        # Add button to navigate to jobs page
        if st.button("📋 Go to My Jobs", type="primary"):
            st.session_state['navigate_to'] = "📋 My Jobs"
            st.rerun()

    except Exception as e:
        st.error(f"❌ Failed to submit job to queue: {str(e)}")
        st.exception(e)


def run_screening(
    pdb_id,
    protein_file,
    library_file,
    library_path,
    project_name,
    client_name,
    affinity_threshold,
    max_workers,
    use_consensus=False
):
    """Run complete screening workflow with progress tracking"""

    # Create temporary directory for this run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(f"data/outputs/streamlit_{timestamp}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Progress tracking
    progress_container = st.container()

    with progress_container:
        progress_bar = st.progress(0, text="Starting screening workflow...")
        status_text = st.empty()
        log_expander = st.expander("📋 Detailed Logs", expanded=False)

    # Save uploaded files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Handle protein file
        if pdb_id:
            target_protein = f"PDB {pdb_id}"
            target_pdb_id = pdb_id
        else:
            # Save uploaded PDB file
            protein_path = tmpdir / protein_file.name
            with open(protein_path, 'wb') as f:
                f.write(protein_file.read())

            target_protein = protein_file.name
            target_pdb_id = protein_file.name.replace('.pdb', '')

        # Handle library - either uploaded file or saved library path
        if library_path:
            # Using saved library - just use the path directly
            final_library_path = library_path
        else:
            # Reset file pointer before reading
            library_file.seek(0)

            # CRITICAL: Parse CSV properly and keep only id,smiles columns
            if library_file.name.endswith('.csv'):
                import io
                from rdkit import Chem

                # Read CSV with pandas - IGNORE comment lines
                file_contents = library_file.read()
                df = pd.read_csv(io.BytesIO(file_contents), dtype=str, comment='#')

                # Normalize column names
                df.columns = df.columns.str.lower().str.strip()

                # Find id and smiles columns
                id_col = None
                smiles_col = None

                for col in df.columns:
                    if col in ['id', 'name', 'mol_id', 'compound_id', 'molecule_id']:
                        id_col = col
                    elif col in ['smiles', 'smile', 'smi', 'structure']:
                        smiles_col = col

                if not id_col or not smiles_col:
                    st.error(f"❌ CSV must have 'id' and 'smiles' columns. Found: {df.columns.tolist()}")
                    return

                # Keep only id and smiles columns, rename to standard
                df_clean = pd.DataFrame({
                    'id': df[id_col],
                    'smiles': df[smiles_col]
                })

                # Drop rows where SMILES column is NaN
                df_clean = df_clean.dropna(subset=['smiles'])

                # Minimal validation: only remove truly empty SMILES
                # Don't validate with RDKit - let the robust 6-method pipeline handle it
                molecule_count = len(df_clean)

                st.info(f"📋 Processing {molecule_count} molecules from library...")

                # Remove only empty/invalid string entries
                empty_count = 0
                for idx, row in df_clean.iterrows():
                    smiles_str = str(row['smiles']).strip()
                    if not smiles_str or smiles_str.lower() in ['nan', 'none', '', 'null']:
                        empty_count += 1

                # Filter out empty entries
                df_clean = df_clean[df_clean['smiles'].notna()]
                df_clean = df_clean[df_clean['smiles'].astype(str).str.strip() != '']
                df_clean = df_clean[~df_clean['smiles'].astype(str).str.lower().isin(['nan', 'none', 'null'])]
                df_clean = df_clean.reset_index(drop=True)

                actual_count = len(df_clean)

                if empty_count > 0:
                    st.warning(f"⚠️ Removed {empty_count} empty SMILES entries")

                if actual_count == 0:
                    st.error("❌ No valid molecules in library!")
                    return

                st.success(f"✅ Loaded {actual_count} molecules for screening")
                st.info("ℹ️ Complex molecules will be processed using advanced 6-method 3D generation pipeline")

                # Save as tab-separated file (required by batch_prepare_ligands)
                smi_path = tmpdir / "library_cleaned.smi"
                with open(smi_path, 'w') as f:
                    for _, row in df_clean.iterrows():
                        f.write(f"{row['smiles']}\t{row['id']}\n")

                final_library_path = str(smi_path)

            elif library_file.name.endswith('.smi'):
                # SMI format - parse as tab-separated
                import io
                file_contents = library_file.read()
                df = pd.read_csv(io.BytesIO(file_contents), sep=r'\s+', header=None, names=['smiles', 'id'], dtype=str)

                # Save as tab-separated for consistency
                smi_path = tmpdir / "library_cleaned.smi"
                with open(smi_path, 'w') as f:
                    for _, row in df.iterrows():
                        f.write(f"{row['smiles']}\t{row['id']}\n")

                final_library_path = str(smi_path)

            else:
                # Unknown format - save as-is and let complete_workflow handle it
                library_file_path = tmpdir / library_file.name
                with open(library_file_path, 'wb') as f:
                    library_file.seek(0)
                    f.write(library_file.read())
                final_library_path = str(library_file_path)

        try:
            # Run complete workflow using subprocess for better control
            status_text.markdown("**Status:** 🔄 Starting workflow...")
            progress_bar.progress(5, text="Initializing...")

            # Build command
            cmd = [
                sys.executable,
                "scripts/complete_workflow.py",
                "--pdb", pdb_id if pdb_id else str(protein_path),
                "--library", str(final_library_path),
                "--output", str(output_dir),
                "--project", project_name,
                "--client", client_name,
                "--threshold", str(affinity_threshold),
                "--workers", str(max_workers)
            ]

            if use_consensus:
                cmd.append("--consensus")

            # Run workflow
            with log_expander:
                log_area = st.empty()
                logs = []

            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True
            )

            # Track progress based on log output
            stage_progress = {
                "Downloading": 10,
                "Detecting": 20,
                "Preparing receptor": 30,
                "Preparing ligand": 40,
                "Running batch": 50,
                "Ranking": 70,
                "ADMET": 80,
                "Generating visualizations": 90,
                "Generating PDF": 95
            }

            current_progress = 5

            for line in process.stdout:
                logs.append(line.strip())

                # Update log display
                with log_expander:
                    log_area.code('\n'.join(logs[-20:]))  # Show last 20 lines

                # Update progress based on log content
                for keyword, progress in stage_progress.items():
                    if keyword.lower() in line.lower():
                        current_progress = progress
                        progress_bar.progress(current_progress, text=f"{keyword}...")
                        break

                # Update status text for key events
                if "✓" in line or "complete" in line.lower():
                    status_text.markdown(f"**Status:** {line.strip()}")

            process.wait()

            if process.returncode != 0:
                progress_bar.progress(0, text="Failed!")
                st.error(f"❌ Screening failed with exit code {process.returncode}")
                with st.expander("Error Details", expanded=True):
                    st.code('\n'.join(logs))
                return

            progress_bar.progress(100, text="Complete!")
            status_text.markdown("**Status:** ✅ Screening complete!")

            # Success message
            st.success(f"""
            🎉 **Screening Complete!**

            Results saved to: `{output_dir}`
            """)

            # Display results summary
            display_results_summary(output_dir, project_name)

            # Display 3D visualization
            display_3d_visualization(output_dir)

        except Exception as e:
            progress_bar.progress(0, text="Failed!")
            status_text.markdown("")
            st.error(f"❌ Screening failed: {str(e)}")
            st.exception(e)


def display_results_summary(output_dir, project_name):
    """Display summary of screening results"""

    st.markdown("---")
    st.subheader("📊 Results Summary")

    # Load results
    csv_path = output_dir / "final_results.csv"

    if csv_path.exists():
        df = pd.read_csv(csv_path)

        # Metrics
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Total Molecules", len(df))

        with col2:
            hits = df[df['binding_affinity'] <= -7.0]
            hit_rate = (len(hits) / len(df)) * 100 if len(df) > 0 else 0
            st.metric("Hits Identified", len(hits), f"{hit_rate:.1f}%")

        with col3:
            best_affinity = df['binding_affinity'].min()
            st.metric("Best Affinity", f"{best_affinity:.2f} kcal/mol")

        with col4:
            if 'lipinski_compliant' in df.columns:
                lipinski_pass = df['lipinski_compliant'].sum()
                pct = (lipinski_pass / len(df)) * 100
                st.metric("Lipinski Compliant", f"{lipinski_pass}/{len(df)}", f"{pct:.0f}%")

        # Top hits table
        st.markdown("### 🏆 Top 10 Hits")

        top_hits = df.nsmallest(10, 'binding_affinity')

        display_cols = ['ligand_id', 'binding_affinity', 'qed_score', 'lipinski_compliant']
        if 'oral_bioavailability' in top_hits.columns:
            display_cols.append('oral_bioavailability')

        display_df = top_hits[display_cols].copy()

        display_df.columns = [
            col.replace('_', ' ').title() for col in display_df.columns
        ]

        st.dataframe(display_df, use_container_width=True, hide_index=True)

        # Visualizations
        st.markdown("### 📊 Visualizations")

        viz_dir = output_dir / "visualizations"
        if viz_dir.exists():
            viz_files = list(viz_dir.glob("*.png"))

            if viz_files:
                # Show visualizations in columns
                for viz_file in sorted(viz_files)[:3]:  # Show first 3
                    st.image(str(viz_file), caption=viz_file.stem.replace('_', ' ').title())

        # Download buttons
        st.markdown("---")
        st.subheader("📥 Download Results")

        col1, col2, col3 = st.columns(3)

        # Find PDF report
        pdf_files = list(output_dir.glob("*_report.pdf"))
        pdf_path = pdf_files[0] if pdf_files else None

        with col1:
            if pdf_path and pdf_path.exists():
                with open(pdf_path, 'rb') as f:
                    st.download_button(
                        "📄 Download PDF Report",
                        data=f,
                        file_name=f"{project_name.replace(' ', '_')}_report.pdf",
                        mime="application/pdf",
                        use_container_width=True
                    )
            else:
                st.button("📄 PDF Not Available", disabled=True, use_container_width=True)

        with col2:
            if csv_path.exists():
                with open(csv_path, 'rb') as f:
                    st.download_button(
                        "📊 Download CSV Data",
                        data=f,
                        file_name=f"{project_name.replace(' ', '_')}_results.csv",
                        mime="text/csv",
                        use_container_width=True
                    )
            else:
                st.button("📊 CSV Not Available", disabled=True, use_container_width=True)

        with col3:
            # Create ZIP of all results
            import zipfile

            zip_path = output_dir / "complete_results.zip"

            try:
                with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
                    for file in output_dir.rglob('*'):
                        if file.is_file() and file != zip_path:
                            zipf.write(file, file.relative_to(output_dir))

                with open(zip_path, 'rb') as f:
                    st.download_button(
                        "📦 Download All Files (ZIP)",
                        data=f,
                        file_name=f"{project_name.replace(' ', '_')}_complete.zip",
                        mime="application/zip",
                        use_container_width=True
                    )
            except Exception as e:
                st.button(f"📦 ZIP Failed: {str(e)}", disabled=True, use_container_width=True)

    else:
        st.warning("Results file not found")


def display_3d_visualization(output_dir):
    """Display 3D protein-ligand complex viewer"""

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from utils.mol_viewer import create_complex_viewer
    import streamlit.components.v1 as components
    import json

    st.markdown("---")
    st.subheader("🔬 3D Protein-Ligand Complex Viewer")

    # Find protein and ligand files
    protein_dir = output_dir / 'proteins'
    ligand_dir = output_dir / 'ligands'
    docking_dir = output_dir / 'docking' / 'results'

    # Get protein file
    pdb_files = list(protein_dir.glob('*_cleaned.pdb'))
    if not pdb_files:
        pdb_files = list(protein_dir.glob('*.pdb'))

    if not pdb_files:
        st.warning("Protein structure file not found.")
        return

    protein_pdb = pdb_files[0]

    # Get docked ligand files
    if docking_dir.exists():
        ligand_files = list(docking_dir.glob('*_docked.pdbqt'))
    else:
        ligand_files = []

    if not ligand_files:
        st.warning("No docked ligand structures found.")
        return

    # Load results to get top hit
    csv_path = output_dir / "final_results.csv"

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        top_hit = df.nsmallest(1, 'binding_affinity').iloc[0]
        top_ligand_id = top_hit['ligand_id']

        # Find corresponding PDBQT file
        top_ligand_file = None
        for lfile in ligand_files:
            if top_ligand_id in lfile.stem:
                top_ligand_file = lfile
                break

        if not top_ligand_file:
            top_ligand_file = ligand_files[0]

        # Display info
        col1, col2, col3 = st.columns(3)

        with col1:
            st.info(f"**Molecule:** {top_hit['ligand_id']}")

        with col2:
            st.info(f"**Affinity:** {top_hit['binding_affinity']:.2f} kcal/mol")

        with col3:
            qed_val = top_hit.get('qed_score', 'N/A')
            if qed_val != 'N/A':
                st.info(f"**QED Score:** {qed_val:.2f}")
            else:
                st.info(f"**QED Score:** {qed_val}")

        # Create 3D viewer
        try:
            # Read pocket info from docking results if available
            pocket_center = None
            pocket_size = None

            # Try to read from docking config or results
            docking_json = output_dir / 'docking' / 'docking_results_detailed.json'
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
            viewer = create_complex_viewer(
                protein_pdb_path=str(protein_pdb),
                ligand_pdbqt_path=str(top_ligand_file),
                pocket_center=pocket_center,
                pocket_size=pocket_size,
                width=800,
                height=600,
                show_pocket=True
            )

            # Render in Streamlit
            components.html(viewer._make_html(), height=650, scrolling=False)

            st.caption("""
            **💡 Interactive Controls:**
            - 🖱️ **Left-click + drag** to rotate
            - 🖱️ **Right-click + drag** to zoom
            - 🖱️ **Scroll wheel** to zoom in/out
            - **Protein** shown in gray cartoon
            - **Ligand** shown in green sticks
            - **Binding pocket** residues highlighted in orange
            """)

        except Exception as e:
            st.error(f"Error creating 3D viewer: {e}")
            st.exception(e)
