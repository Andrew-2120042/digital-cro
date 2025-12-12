"""
Molecule Libraries Management Page

Save, browse, and manage molecule libraries.
"""

import streamlit as st
from pathlib import Path
import sys
import pandas as pd
from typing import Dict

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.utils.library_manager import LibraryManager


def show():
    """Display libraries management page"""

    st.markdown('<p class="main-header">📚 My Molecule Libraries</p>', unsafe_allow_html=True)

    st.write("""
    Save and organize your molecule libraries for quick reuse in screenings.
    """)

    manager = LibraryManager()

    # Statistics
    stats = manager.get_statistics()

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Total Libraries", stats['total_libraries'])

    with col2:
        st.metric("Total Molecules", f"{stats['total_molecules']:,}")

    with col3:
        if stats['most_used']:
            st.metric("Most Used", stats['most_used'])
        else:
            st.metric("Most Used", "N/A")

    with col4:
        if stats['largest']:
            st.metric("Largest Library", stats['largest'])
        else:
            st.metric("Largest Library", "N/A")

    st.markdown("---")

    # Tabs
    tab1, tab2 = st.tabs(["📚 Browse Libraries", "➕ Add New Library"])

    with tab1:
        browse_libraries(manager)

    with tab2:
        add_library(manager)


def browse_libraries(manager: LibraryManager):
    """Browse and manage existing libraries"""

    libraries = manager.list_libraries()

    if not libraries:
        st.info("No libraries saved yet. Add your first library in the '➕ Add New Library' tab!")
        return

    # Search and filter
    col1, col2 = st.columns([2, 1])

    with col1:
        search_query = st.text_input(
            "🔍 Search libraries",
            placeholder="Search by name or description...",
            key="library_search"
        )

    with col2:
        sort_by = st.selectbox(
            "Sort by",
            ["Name", "Size", "Recently Used", "Use Count"],
            key="library_sort"
        )

    # Apply search
    if search_query:
        libraries = manager.search_libraries(query=search_query)

    # Sort
    if sort_by == "Name":
        libraries.sort(key=lambda x: x['name'])
    elif sort_by == "Size":
        libraries.sort(key=lambda x: x['num_molecules'], reverse=True)
    elif sort_by == "Recently Used":
        libraries.sort(key=lambda x: x.get('last_used', ''), reverse=True)
    elif sort_by == "Use Count":
        libraries.sort(key=lambda x: x.get('use_count', 0), reverse=True)

    st.write(f"**{len(libraries)} libraries found**")

    # Display libraries as cards
    for lib in libraries:
        with st.expander(f"📚 {lib['name']} ({lib['num_molecules']:,} molecules)"):
            display_library_card(lib, manager)


def display_library_card(lib: Dict, manager: LibraryManager):
    """Display library details as a card"""

    col1, col2 = st.columns([2, 1])

    with col1:
        st.write(f"**Description:** {lib.get('description', 'No description')}")

        if lib.get('tags'):
            tags_html = ' '.join([f'<span style="background-color: #1f4788; color: white; padding: 2px 8px; border-radius: 3px; margin-right: 5px;">{tag}</span>' for tag in lib['tags']])
            st.markdown(tags_html, unsafe_allow_html=True)

        st.write(f"**File:** `{lib['file_name']}`")
        st.write(f"**Molecules:** {lib['num_molecules']:,}")

        if lib.get('stats'):
            stats = lib['stats']
            if 'mean_mw' in stats:
                st.write(f"**Avg. MW:** {stats['mean_mw']:.1f} Da")
            if 'mean_logp' in stats:
                st.write(f"**Avg. LogP:** {stats['mean_logp']:.2f}")

        st.write(f"**Created:** {lib['created_at'][:10]}")
        if lib.get('last_used'):
            st.write(f"**Last Used:** {lib['last_used'][:10]}")
        st.write(f"**Use Count:** {lib.get('use_count', 0)}")

    with col2:
        # Actions
        if st.button("🔬 Use in Screening", key=f"use_{lib['library_id']}", use_container_width=True):
            st.session_state['selected_library'] = lib['library_id']
            st.success("Library selected! Go to 'New Screening' page to use it.")

        if st.button("🗑️ Delete", key=f"delete_lib_{lib['library_id']}", use_container_width=True):
            if st.session_state.get(f"confirm_delete_{lib['library_id']}"):
                manager.delete_library(lib['library_id'])
                st.success("Library deleted!")
                st.session_state.pop(f"confirm_delete_{lib['library_id']}")
                st.rerun()
            else:
                st.session_state[f"confirm_delete_{lib['library_id']}"] = True
                st.warning("Click again to confirm deletion")


def add_library(manager: LibraryManager):
    """Add new library form"""

    st.subheader("➕ Add New Molecule Library")

    with st.form("add_library_form"):
        library_name = st.text_input(
            "Library Name *",
            placeholder="e.g., FDA Approved Drugs"
        )

        description = st.text_area(
            "Description",
            placeholder="Optional: Describe what this library contains"
        )

        tags_input = st.text_input(
            "Tags (comma-separated)",
            placeholder="e.g., kinase-inhibitors, FDA-approved, natural-products"
        )

        library_file = st.file_uploader(
            "Upload Library File *",
            type=['csv', 'smi', 'txt'],
            help="CSV with 'smiles' column or SMI format"
        )

        compute_stats = st.checkbox(
            "Compute molecular properties",
            value=True,
            help="Calculate average MW, LogP, etc. (takes a few seconds)"
        )

        submitted = st.form_submit_button("💾 Save Library", use_container_width=True, type="primary")

    if submitted:
        # Validation
        if not library_name:
            st.error("Please enter a library name")
            return

        if not library_file:
            st.error("Please upload a library file")
            return

        # Parse tags
        tags = [t.strip() for t in tags_input.split(',') if t.strip()] if tags_input else []

        # Save file temporarily
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir) / library_file.name

            with open(tmp_path, 'wb') as f:
                f.write(library_file.read())

            try:
                # Save library
                library_id = manager.save_library(
                    name=library_name,
                    file_path=str(tmp_path),
                    description=description,
                    tags=tags,
                    compute_stats=compute_stats
                )

                st.success(f"""
                ✅ **Library Saved!**

                Name: {library_name}
                ID: `{library_id}`

                You can now use this library in your screenings!
                """)

            except Exception as e:
                st.error(f"Error saving library: {e}")
