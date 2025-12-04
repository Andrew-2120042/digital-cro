"""
Job Queue Management Page

View and manage queued, running, and completed jobs.
"""

import streamlit as st
from pathlib import Path
import sys
from datetime import datetime
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.utils.job_queue import JobQueue, JobStatus


def show():
    """Display job management page"""

    st.markdown('<p class="main-header">📋 My Jobs</p>', unsafe_allow_html=True)

    st.write("""
    View and manage your screening jobs. Jobs are processed in the order they're submitted (FIFO).
    """)

    queue = JobQueue()

    # Queue statistics
    stats = queue.get_queue_stats()

    col1, col2, col3, col4, col5 = st.columns(5)

    with col1:
        st.metric("Total Jobs", stats['total_jobs'])

    with col2:
        st.metric("⏳ Queued", stats['queued'])

    with col3:
        st.metric("🔄 Running", stats['running'])

    with col4:
        st.metric("✅ Complete", stats['complete'])

    with col5:
        st.metric("❌ Failed", stats['failed'])

    st.markdown("---")

    # Job list tabs
    tab1, tab2, tab3, tab4 = st.tabs([
        "All Jobs",
        "⏳ Queued",
        "✅ Complete",
        "❌ Failed"
    ])

    with tab1:
        display_job_list(queue, status=None, tab_id="all_jobs_tab")

    with tab2:
        display_job_list(queue, status=JobStatus.QUEUED.value, tab_id="queued_tab")

    with tab3:
        display_job_list(queue, status=JobStatus.COMPLETE.value, tab_id="complete_tab")

    with tab4:
        display_job_list(queue, status=JobStatus.FAILED.value, tab_id="failed_tab")


def display_job_list(queue: JobQueue, status: str = None, tab_id: str = "default"):
    """Display list of jobs"""

    jobs = queue.list_jobs(status=status, limit=100)

    if not jobs:
        st.info("No jobs found")
        return

    st.write(f"**{len(jobs)} job(s)**")

    # Convert to DataFrame for display
    job_data = []

    for job in jobs:
        # Parse timestamps
        created = datetime.fromisoformat(job['created_at'])

        # Calculate duration if completed
        duration = None
        if job.get('completed_at'):
            completed = datetime.fromisoformat(job['completed_at'])
            duration = (completed - created).total_seconds() / 60  # minutes

        job_data.append({
            'Job ID': job['job_id'][:8] + '...',
            'Project': job['project_name'],
            'Type': job['job_type'].replace('_', ' ').title(),
            'Status': job['status'],
            'Created': created.strftime('%Y-%m-%d %H:%M'),
            'Duration (min)': f"{duration:.1f}" if duration else '-',
            'Targets': ', '.join(job['target_proteins']),
            '_full_job': job  # Store full job for actions
        })

    df = pd.DataFrame(job_data)

    # Display table
    st.dataframe(
        df.drop('_full_job', axis=1),
        use_container_width=True,
        hide_index=True
    )

    # Job actions
    st.markdown("---")
    st.subheader("Job Actions")

    # Select job
    job_options = {f"{j['Job ID']} - {j['Project']}": j['_full_job'] for j in job_data}

    # Create unique key for selectbox using tab_id to ensure uniqueness across tabs
    selectbox_key = f"job_selector_{tab_id}"

    selected_job_label = st.selectbox(
        "Select job for actions:",
        options=list(job_options.keys()),
        key=selectbox_key
    )

    if selected_job_label:
        selected_job = job_options[selected_job_label]

        col1, col2, col3 = st.columns(3)

        # Create unique button keys using tab_id
        with col1:
            if st.button("📊 View Details", use_container_width=True, key=f"view_{tab_id}"):
                display_job_details(selected_job)

        with col2:
            if selected_job['status'] == JobStatus.QUEUED.value:
                if st.button("❌ Cancel Job", use_container_width=True, key=f"cancel_{tab_id}"):
                    try:
                        queue.cancel_job(selected_job['job_id'])
                        st.success("Job cancelled")
                        st.rerun()
                    except Exception as e:
                        st.error(f"Error: {e}")

        with col3:
            if selected_job['status'] in [JobStatus.COMPLETE.value, JobStatus.FAILED.value, JobStatus.CANCELLED.value]:
                if st.button("🗑️ Delete Job", use_container_width=True, key=f"delete_{tab_id}"):
                    if st.session_state.get(f'confirm_delete_{tab_id}') == selected_job['job_id']:
                        try:
                            queue.delete_job(selected_job['job_id'])
                            st.success("Job deleted")
                            st.session_state.pop(f'confirm_delete_{tab_id}')
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error: {e}")
                    else:
                        st.session_state[f'confirm_delete_{tab_id}'] = selected_job['job_id']
                        st.warning("Click again to confirm deletion")


def display_job_details(job: dict):
    """Display detailed information about a job"""

    st.markdown("---")
    st.subheader(f"Job Details: {job['project_name']}")

    col1, col2 = st.columns(2)

    with col1:
        st.write("**Job Information:**")
        st.write(f"- Job ID: `{job['job_id']}`")
        st.write(f"- Type: {job['job_type'].replace('_', ' ').title()}")
        st.write(f"- Status: {job['status']}")
        st.write(f"- Client: {job['client_name']}")

    with col2:
        st.write("**Timestamps:**")
        st.write(f"- Created: {job['created_at']}")
        if job.get('started_at'):
            st.write(f"- Started: {job['started_at']}")
        if job.get('completed_at'):
            st.write(f"- Completed: {job['completed_at']}")

    st.write("**Parameters:**")
    st.write(f"- Targets: {', '.join(job['target_proteins'])}")
    st.write(f"- Library: `{job['ligand_library_path']}`")
    st.write(f"- Threshold: {job['affinity_threshold']} kcal/mol")
    st.write(f"- Workers: {job['max_workers']}")

    # Results
    if job.get('result_summary'):
        st.write("**Results:**")
        st.json(job['result_summary'])

    # Error message
    if job.get('error_message'):
        st.error(f"**Error:** {job['error_message']}")

    # Download results
    if job['status'] == JobStatus.COMPLETE.value:
        output_dir = Path(job['output_dir'])

        if output_dir.exists():
            st.markdown("**Download Results:**")

            # Find PDF reports
            pdf_files = list(output_dir.rglob("*_report.pdf"))

            for pdf_file in pdf_files:
                with open(pdf_file, 'rb') as f:
                    st.download_button(
                        f"📄 {pdf_file.name}",
                        data=f,
                        file_name=pdf_file.name,
                        mime="application/pdf"
                    )
