"""
Multi-Target Screening Page

Allows users to screen molecule library against multiple protein targets.
"""

import streamlit as st
from pathlib import Path
import pandas as pd
import tempfile
import sys
from datetime import datetime
import zipfile
import json

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.utils.multi_target import run_multi_target_screening
import plotly.express as px
import plotly.graph_objects as go


def show():
    """Display multi-target screening page"""

    st.markdown('<p class="main-header">🎯 Multi-Target Screening</p>', unsafe_allow_html=True)

    st.info("""
    **Screen your molecule library against multiple protein targets at once!**

    Perfect for:
    - Drug repurposing studies
    - Selectivity analysis
    - Multi-target drug discovery
    - Off-target prediction
    """)

    st.markdown("---")

    # Configuration form
    with st.form("multi_target_config"):
        st.subheader("📋 Project Configuration")

        col1, col2 = st.columns(2)

        with col1:
            project_name = st.text_input(
                "Project Name",
                value="Multi-Target Screening",
                help="Name for this screening campaign"
            )

            client_name = st.text_input(
                "Client/Organization Name",
                value="Research Lab",
                help="Your organization name"
            )

        with col2:
            affinity_threshold = st.number_input(
                "Binding Affinity Threshold (kcal/mol)",
                value=-7.0,
                step=0.1,
                help="Molecules with affinity ≤ this are hits"
            )

            max_workers = st.number_input(
                "Parallel Workers (per target)",
                value=4,
                min_value=1,
                max_value=16,
                help="CPU cores for docking"
            )

        st.markdown("---")

        # Target proteins section
        st.subheader("🧬 Target Proteins")

        st.write("**Select targets to screen against:**")

        # Method 1: Enter PDB IDs
        pdb_ids_input = st.text_area(
            "PDB IDs (one per line or comma-separated)",
            value="1HSG\n3CL5",
            height=100,
            help="Enter PDB IDs separated by newlines or commas"
        )

        # Parse PDB IDs
        if pdb_ids_input:
            # Split by newlines or commas, clean up
            pdb_ids = [
                p.strip().upper()
                for p in pdb_ids_input.replace(',', '\n').split('\n')
                if p.strip()
            ]

            if pdb_ids:
                st.success(f"✓ {len(pdb_ids)} targets selected: {', '.join(pdb_ids)}")
            else:
                pdb_ids = []
        else:
            pdb_ids = []

        st.markdown("---")

        # Molecule library section
        st.subheader("💊 Molecule Library")

        library_file = st.file_uploader(
            "Upload Molecule Library",
            type=['csv', 'smi', 'txt'],
            help="CSV with 'smiles' column or SMI format (SMILES<tab>ID)"
        )

        if library_file:
            # Detect file type and preview
            file_ext = Path(library_file.name).suffix.lower()

            if file_ext == '.csv':
                df_preview = pd.read_csv(library_file)
                library_file.seek(0)
            elif file_ext in ['.smi', '.txt']:
                # Read SMI file
                lines = library_file.read().decode('utf-8').strip().split('\n')
                library_file.seek(0)

                smiles_list = []
                id_list = []

                for line in lines:
                    if '\t' in line:
                        parts = line.split('\t')
                        smiles_list.append(parts[0])
                        id_list.append(parts[1] if len(parts) > 1 else f'mol_{len(id_list)}')
                    else:
                        smiles_list.append(line.strip())
                        id_list.append(f'mol_{len(id_list)}')

                df_preview = pd.DataFrame({
                    'id': id_list,
                    'smiles': smiles_list
                })

            st.write(f"**Preview:** {len(df_preview)} molecules loaded")
            st.dataframe(df_preview.head(10), use_container_width=True)

        st.markdown("---")

        # Advanced options
        with st.expander("⚙️ Advanced Options"):
            run_parallel = st.checkbox(
                "Run targets in parallel",
                value=False,
                help="Faster but uses more RAM (recommended only for <5 targets)"
            )

        st.markdown("---")

        # Job queue option
        submit_to_queue = st.checkbox(
            "📋 Submit to job queue (run in background)",
            value=False,
            help="Queue this job for background processing instead of running now"
        )

        # Submit button
        submit_label = "📋 Add to Queue" if submit_to_queue else "🚀 Run Multi-Target Screening"
        submitted = st.form_submit_button(
            submit_label,
            use_container_width=True,
            type="primary"
        )

    # Process if submitted
    if submitted:
        # Clear any previous results when starting new screening
        if 'multi_target_results' in st.session_state:
            del st.session_state['multi_target_results']
            del st.session_state['multi_target_output_dir']
            del st.session_state['multi_target_pdb_ids']

        # Validation
        errors = []

        if not pdb_ids:
            errors.append("Please enter at least 2 PDB IDs")
        elif len(pdb_ids) < 2:
            errors.append("Multi-target screening requires at least 2 targets")

        if not library_file:
            errors.append("Please upload a molecule library")

        if errors:
            for error in errors:
                st.error(error)
        else:
            if submit_to_queue:
                # Submit to queue instead of running immediately
                submit_to_multi_target_queue(
                    pdb_ids=pdb_ids,
                    library_file=library_file,
                    project_name=project_name,
                    client_name=client_name,
                    affinity_threshold=affinity_threshold,
                    max_workers=max_workers
                )
            else:
                # Run multi-target screening immediately
                run_multi_target(
                    pdb_ids=pdb_ids,
                    library_file=library_file,
                    project_name=project_name,
                    client_name=client_name,
                    affinity_threshold=affinity_threshold,
                    max_workers=max_workers,
                    run_parallel=run_parallel
                )

    # CRITICAL FIX: Display cached results if they exist (for button clicks after initial form submission)
    elif 'multi_target_results' in st.session_state:
        display_multi_target_results(
            st.session_state['multi_target_results'],
            st.session_state['multi_target_output_dir'],
            st.session_state['multi_target_pdb_ids']
        )


def submit_to_multi_target_queue(
    pdb_ids,
    library_file,
    project_name,
    client_name,
    affinity_threshold,
    max_workers
):
    """Submit multi-target screening job to queue for background processing"""

    from scripts.utils.job_queue import JobQueue

    try:
        queue = JobQueue()

        # Save library file permanently to job queue directory
        library_save_path = Path("data/job_queue/libraries") / library_file.name
        library_save_path.parent.mkdir(parents=True, exist_ok=True)

        with open(library_save_path, 'wb') as f:
            f.write(library_file.read())

        # Convert .smi or .txt to CSV if needed
        file_ext = Path(library_file.name).suffix.lower()

        if file_ext in ['.smi', '.txt']:
            # Parse SMI file
            df = pd.read_csv(library_save_path, sep=r'\s+', header=None, names=['smiles', 'id'])
            csv_path = library_save_path.parent / f"{library_save_path.stem}.csv"
            df.to_csv(csv_path, index=False)
            library_save_path = csv_path

        # Submit job to queue
        job_id = queue.submit_multi_target_job(
            target_proteins=pdb_ids,
            ligand_library_path=str(library_save_path),
            project_name=project_name,
            client_name=client_name,
            affinity_threshold=affinity_threshold,
            max_workers=max_workers
        )

        st.success(f"""
        ✅ **Multi-Target Job Submitted to Queue!**

        **Job ID:** `{job_id}`
        **Targets:** {', '.join(pdb_ids)} ({len(pdb_ids)} targets)

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


def run_multi_target(
    pdb_ids,
    library_file,
    project_name,
    client_name,
    affinity_threshold,
    max_workers,
    run_parallel
):
    """Execute multi-target screening with progress tracking"""

    # Create output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(f"data/outputs/multi_target_{timestamp}")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save library file
    library_path = output_dir / library_file.name
    with open(library_path, 'wb') as f:
        f.write(library_file.getvalue())

    # Progress tracking
    progress_container = st.container()

    with progress_container:
        st.markdown("### 🔄 Screening Progress")

        # Overall progress
        progress_bar = st.progress(0)
        status_text = st.empty()

        # Per-target status
        target_status = {}
        target_cols = st.columns(min(len(pdb_ids), 4))

        for idx, pdb_id in enumerate(pdb_ids):
            col_idx = idx % len(target_cols)
            with target_cols[col_idx]:
                target_status[pdb_id] = st.empty()
                target_status[pdb_id].info(f"**{pdb_id}**\n⏳ Queued")

        st.markdown("---")

    try:
        # Update status
        status_text.text(f"Starting multi-target screening for {len(pdb_ids)} targets...")
        progress_bar.progress(5)

        # Run screening
        results = run_multi_target_screening(
            target_proteins=pdb_ids,
            ligand_library_path=str(library_path),
            output_dir=str(output_dir),
            project_name=project_name,
            client_name=client_name,
            affinity_threshold=affinity_threshold,
            max_workers_per_target=max_workers,
            run_targets_parallel=run_parallel
        )

        progress_bar.progress(100)
        status_text.text("✅ Multi-target screening complete!")

        # Update target status
        for pdb_id in pdb_ids:
            if pdb_id in results['target_results']:
                data = results['target_results'][pdb_id]

                if data.get('success', False):
                    target_status[pdb_id].success(
                        f"**{pdb_id}**\n"
                        f"✅ {data['num_hits']} hits\n"
                        f"Best: {data['best_affinity']:.2f}"
                    )
                else:
                    target_status[pdb_id].error(f"**{pdb_id}**\n❌ Failed")

        # Success message
        st.success(f"""
        🎉 **Multi-Target Screening Complete!**

        - Targets screened: {len(pdb_ids)}
        - Total molecules: {results['summary']['total_molecules']}
        - Results saved to: `{output_dir.name}`
        """)

        # CRITICAL FIX: Store results in session state so they persist across reruns
        st.session_state['multi_target_results'] = results
        st.session_state['multi_target_output_dir'] = str(output_dir)
        st.session_state['multi_target_pdb_ids'] = pdb_ids

        # Display results
        display_multi_target_results(results, output_dir, pdb_ids)

    except Exception as e:
        progress_bar.progress(0)
        status_text.text("")
        st.error(f"❌ Multi-target screening failed: {str(e)}")
        st.exception(e)


def display_multi_target_results(results, output_dir, pdb_ids):
    """Display multi-target screening results"""

    st.markdown("---")
    st.markdown("## 📊 Results Summary")

    # Summary metrics
    summary = results['summary']

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Targets Screened", summary['successful_targets'])

    with col2:
        st.metric("Total Molecules", summary['total_molecules'])

    with col3:
        if 'selectivity' in summary:
            st.metric("Selective Molecules", summary['selectivity']['selective_molecules'])

    with col4:
        if 'multi_target_hits' in summary:
            multi_hits = summary['multi_target_hits'].get('2_targets', 0) + \
                        summary['multi_target_hits'].get('3+_targets', 0)
            st.metric("Multi-Target Hits", multi_hits)

    # Tabs for different views
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "🔥 Comparative Heatmap",
        "🎯 Selectivity Analysis",
        "🧬 Pocket Comparison",
        "📋 All Results",
        "📥 Downloads"
    ])

    with tab1:
        display_comparative_heatmap(results['comparative_df'], pdb_ids)

    with tab2:
        display_selectivity_analysis(results['selectivity_df'])

    with tab3:
        display_pocket_comparison_tab(output_dir, pdb_ids)

    with tab4:
        st.subheader("Complete Results Table")
        st.dataframe(results['comparative_df'], use_container_width=True)

        # Download comparative CSV
        csv = results['comparative_df'].to_csv(index=False)
        st.download_button(
            "📊 Download Comparative Results (CSV)",
            data=csv,
            file_name="multi_target_comparative.csv",
            mime="text/csv"
        )

    with tab5:
        display_downloads(results, output_dir, pdb_ids)


def display_comparative_heatmap(df, pdb_ids):
    """Display interactive heatmap of binding affinities"""

    st.subheader("🔥 Binding Affinity Heatmap")

    # Check if dataframe is empty
    if df.empty or len(df) == 0:
        st.warning("No molecules were successfully docked. Please check your input molecules and try again.")
        return

    st.write("""
    **How to read:**
    - Darker colors = stronger binding (more negative affinity)
    - Each row is a molecule, each column is a target
    - Click cells to see exact values
    """)

    # Prepare data for heatmap
    affinity_cols = [f'{pdb_id}_affinity' for pdb_id in pdb_ids]

    # Get top N molecules
    max_molecules = min(50, len(df))
    default_molecules = min(20, len(df))

    if max_molecules < 5:
        # If less than 5 molecules, just show all
        top_n = len(df)
    else:
        top_n = st.slider("Number of molecules to show", 5, max_molecules, default_molecules)

    df_top = df.nsmallest(top_n, 'best_affinity')

    # Create heatmap data
    heatmap_data = df_top[affinity_cols].values
    molecules = df_top['ligand_id'].values

    # Create Plotly heatmap
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data,
        x=pdb_ids,
        y=molecules,
        colorscale='RdYlGn_r',  # Red = weak, Green = strong
        reversescale=False,
        colorbar=dict(title="Binding Affinity<br>(kcal/mol)"),
        hovertemplate='Molecule: %{y}<br>Target: %{x}<br>Affinity: %{z:.2f} kcal/mol<extra></extra>'
    ))

    fig.update_layout(
        title=f"Top {top_n} Molecules Across {len(pdb_ids)} Targets",
        xaxis_title="Target Protein",
        yaxis_title="Molecule",
        height=max(400, top_n * 20),
        width=800
    )

    st.plotly_chart(fig, use_container_width=True)


def display_selectivity_analysis(selectivity_df):
    """Display selectivity analysis results"""

    st.subheader("🎯 Selectivity Analysis")

    if selectivity_df.empty:
        st.warning("No selectivity data available")
        return

    # Summary
    col1, col2 = st.columns(2)

    with col1:
        # Pie chart: Selective vs Promiscuous
        counts = selectivity_df['classification'].value_counts()

        fig = px.pie(
            values=counts.values,
            names=counts.index,
            title="Selectivity Classification",
            color_discrete_map={
                'Selective': '#28a745',
                'Promiscuous': '#ffc107'
            }
        )

        st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Selectivity score distribution
        fig = px.histogram(
            selectivity_df,
            x='selectivity_score',
            nbins=20,
            title="Selectivity Score Distribution",
            labels={'selectivity_score': 'Selectivity Score (kcal/mol)'}
        )

        fig.add_vline(
            x=2.0,
            line_dash="dash",
            line_color="red",
            annotation_text="Threshold (2.0)"
        )

        st.plotly_chart(fig, use_container_width=True)

    # Selective molecules table
    st.markdown("### 🏆 Selective Molecules")

    selective = selectivity_df[selectivity_df['classification'] == 'Selective']

    if not selective.empty:
        display_cols = [
            'ligand_id', 'selective_for', 'selectivity_score',
            'best_affinity', 'second_best_affinity'
        ]

        st.dataframe(
            selective[display_cols].sort_values('selectivity_score', ascending=False),
            use_container_width=True,
            hide_index=True
        )
    else:
        st.info("No highly selective molecules found (all are promiscuous binders)")


def display_downloads(results, output_dir, pdb_ids):
    """Display download options for all reports"""

    st.subheader("📥 Download Results")

    output_dir = Path(output_dir)

    # Individual target reports
    st.markdown("### 📄 Individual Target Reports")

    cols = st.columns(min(len(pdb_ids), 3))

    for idx, pdb_id in enumerate(pdb_ids):
        col_idx = idx % len(cols)

        with cols[col_idx]:
            target_dir = output_dir / f"target_{pdb_id}"
            pdf_files = list(target_dir.glob("*_report.pdf"))

            if pdf_files:
                pdf_path = pdf_files[0]

                with open(pdf_path, 'rb') as f:
                    st.download_button(
                        f"📄 {pdb_id} Report",
                        data=f,
                        file_name=f"{pdb_id}_report.pdf",
                        mime="application/pdf",
                        use_container_width=True,
                        key=f"pdf_{pdb_id}"
                    )

    st.markdown("---")

    # All results as ZIP
    st.markdown("### 📦 Complete Results Package")

    col1, col2 = st.columns(2)

    with col1:
        # Comparative results CSV
        csv_path = output_dir / "multi_target_results.csv"

        if csv_path.exists():
            with open(csv_path, 'rb') as f:
                st.download_button(
                    "📊 Comparative Results (CSV)",
                    data=f,
                    file_name="multi_target_results.csv",
                    mime="text/csv",
                    use_container_width=True
                )

    with col2:
        # Selectivity analysis CSV
        sel_path = output_dir / "selectivity_analysis.csv"

        if sel_path.exists():
            with open(sel_path, 'rb') as f:
                st.download_button(
                    "🎯 Selectivity Analysis (CSV)",
                    data=f,
                    file_name="selectivity_analysis.csv",
                    mime="text/csv",
                    use_container_width=True
                )

    st.markdown("---")

    # Everything as ZIP
    zip_path = output_dir / "complete_results.zip"

    if not zip_path.exists():
        with zipfile.ZipFile(zip_path, 'w') as zipf:
            for file in output_dir.rglob('*'):
                if file.is_file() and file != zip_path:
                    zipf.write(file, file.relative_to(output_dir))

    with open(zip_path, 'rb') as f:
        st.download_button(
            "📦 Download Everything (ZIP)",
            data=f,
            file_name="multi_target_complete.zip",
            mime="application/zip",
            use_container_width=True,
            type="primary"
        )


def display_pocket_comparison_tab(output_dir, target_proteins):
    """Display protein pocket comparison analysis"""

    st.subheader("🧬 Protein Pocket Comparison")

    st.info("""
    **Pocket Comparison** analyzes structural similarity between binding sites
    to predict cross-reactivity, selectivity, and drug repurposing opportunities.
    """)

    output_path = Path(output_dir)
    pocket_dir = output_path / "pocket_comparison"

    # Check if analysis already exists
    def check_analysis_exists():
        """Helper to check if analysis is complete"""
        return pocket_dir.exists() and (pocket_dir / "pocket_comparison.json").exists()

    analysis_exists = check_analysis_exists()  # Initial check

    if not analysis_exists:
        # Analysis not yet generated
        st.warning("⚠️ Pocket comparison analysis has not been generated yet.")

        col1, col2, col3 = st.columns([1, 2, 1])
        with col2:
            analyze_clicked = st.button(
                "🔬 Analyze Pocket Similarity",
                key="analyze_pocket_button",
                type="primary",
                use_container_width=True
            )

        if analyze_clicked:
            with st.spinner("🔬 Analyzing protein pocket similarity..."):
                try:
                    from scripts.utils.pocket_comparison import compare_multi_target_pockets

                    # Run comparison
                    pocket_dir.mkdir(exist_ok=True)

                    results = compare_multi_target_pockets(
                        str(output_path),
                        output_dir=str(pocket_dir)
                    )

                    # Generate visualizations
                    from scripts.utils.pocket_viz import (
                        plot_similarity_matrix,
                        plot_pocket_properties,
                        plot_comparison_summary
                    )

                    if results.get('comparison_matrix'):
                        plot_similarity_matrix(
                            results['comparison_matrix'],
                            results['target_names'],
                            str(pocket_dir / "similarity_matrix.png")
                        )

                    if results.get('pockets'):
                        plot_pocket_properties(
                            results['pockets'],
                            str(pocket_dir / "pocket_properties.png")
                        )

                    if results.get('comparisons'):
                        plot_comparison_summary(
                            results['comparisons'],
                            results['selectivity'],
                            str(pocket_dir / "comparison_summary.png")
                        )

                    st.success("✅ Pocket analysis complete!")

                    # RE-CHECK after generation (CRITICAL FIX!)
                    analysis_exists = check_analysis_exists()

                except Exception as e:
                    st.error(f"❌ Analysis failed: {str(e)}")
                    import traceback
                    with st.expander("📋 Error Details"):
                        st.code(traceback.format_exc())

    # Display results if analysis exists
    # This check now happens AFTER potential generation above
    if analysis_exists:
        display_pocket_comparison_results(pocket_dir)

        # Regenerate option
        st.write("")
        col1, col2, col3 = st.columns([1, 2, 1])
        with col2:
            if st.button("🔄 Regenerate Analysis", key="regen_pocket"):
                import shutil
                shutil.rmtree(pocket_dir)
                st.rerun()  # Only rerun for regenerate!


def display_pocket_comparison_results(pocket_dir):
    """Display pocket comparison results"""

    # Load results
    results_file = pocket_dir / "pocket_comparison.json"

    if results_file.exists():
        with open(results_file, 'r') as f:
            results = json.load(f)

        # Selectivity assessment
        st.markdown("### 🎯 Selectivity Assessment")

        selectivity = results['selectivity']

        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric("Average Similarity", f"{selectivity['avg_similarity']:.2f}")

        with col2:
            st.metric("Comparisons", selectivity['num_comparisons'])

        with col3:
            # Selectivity score (inverse of similarity)
            selectivity_score = 1 - selectivity['avg_similarity']
            st.metric("Selectivity Score", f"{selectivity_score:.2f}")

        # Level indicator
        if selectivity_score >= 0.7:
            st.success(f"✅ {selectivity['selectivity_level']}")
        elif selectivity_score >= 0.4:
            st.warning(f"⚠️ {selectivity['selectivity_level']}")
        else:
            st.error(f"❌ {selectivity['selectivity_level']}")

    # Recommendations
    rec_file = pocket_dir / "repurposing_recommendations.txt"

    if rec_file.exists():
        with open(rec_file, 'r') as f:
            recommendations = f.read()

        st.markdown("### 💡 Repurposing Recommendations")
        st.info(recommendations)

    # Visualizations
    st.markdown("### 📊 Analysis Visualizations")

    viz_tab1, viz_tab2, viz_tab3 = st.tabs([
        "Similarity Matrix",
        "Pocket Properties",
        "Summary"
    ])

    with viz_tab1:
        matrix_img = pocket_dir / "similarity_matrix.png"
        if matrix_img.exists():
            st.image(str(matrix_img), use_container_width=True)
            st.caption("Heatmap showing pairwise pocket similarity scores")

    with viz_tab2:
        props_img = pocket_dir / "pocket_properties.png"
        if props_img.exists():
            st.image(str(props_img), use_container_width=True)
            st.caption("Chemical and physical properties of binding pockets")

    with viz_tab3:
        summary_img = pocket_dir / "comparison_summary.png"
        if summary_img.exists():
            st.image(str(summary_img), use_container_width=True)
            st.caption("Comprehensive pocket comparison summary")

    # Download
    st.markdown("### 📥 Download Analysis")

    col1, col2 = st.columns(2)

    with col1:
        if results_file.exists():
            with open(results_file, 'rb') as f:
                st.download_button(
                    "📊 Download Results (JSON)",
                    data=f,
                    file_name="pocket_comparison.json",
                    mime="application/json",
                    use_container_width=True,
                    key="download_pocket_json"
                )

    with col2:
        if rec_file.exists():
            with open(rec_file, 'rb') as f:
                st.download_button(
                    "📄 Download Recommendations (TXT)",
                    data=f,
                    file_name="repurposing_recommendations.txt",
                    mime="text/plain",
                    use_container_width=True,
                    key="download_pocket_recs"
                )
