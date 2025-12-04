"""
Digital CRO Platform - Streamlit Web Interface

Main application with:
- Home page with project overview
- Screening workflow interface
- Results browser
- About/documentation
"""

import streamlit as st
from pathlib import Path
import sys

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

# Page configuration
st.set_page_config(
    page_title="Digital CRO - Drug Discovery Platform",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS - Dark Theme
st.markdown("""
<style>
    /* Main app background */
    .main {
        background-color: #0e1117;
    }

    /* Headers */
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #4da6ff;
        text-align: center;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #66b3ff;
        text-align: center;
        margin-bottom: 2rem;
    }

    /* Dark metric cards */
    .metric-card {
        background-color: #1e2130;
        padding: 1.5rem;
        border-radius: 0.5rem;
        border-left: 4px solid #4da6ff;
        color: #e0e0e0;
    }
    .metric-card h3 {
        color: #4da6ff;
        margin-bottom: 0.5rem;
    }
    .metric-card p {
        color: #b0b0b0;
    }

    /* Success and warning boxes */
    .success-box {
        background-color: #1a3d2e;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #28a745;
        color: #a8e6a8;
    }
    .warning-box {
        background-color: #3d3520;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #ffc107;
        color: #ffe066;
    }

    /* Buttons */
    .stButton>button {
        background-color: #4da6ff;
        color: #0e1117;
        font-weight: bold;
        border-radius: 0.3rem;
        padding: 0.5rem 2rem;
        border: none;
    }
    .stButton>button:hover {
        background-color: #66b3ff;
        color: #0e1117;
    }

    /* Forms - Dark theme */
    .stForm {
        border: 1px solid #2d3748;
        border-radius: 0.5rem;
        padding: 1.5rem;
        background-color: #1a1f2e;
    }

    /* Form inputs */
    .stTextInput > div > div > input,
    .stNumberInput > div > div > input,
    .stSelectbox > div > div > select {
        background-color: #2d3748;
        color: #e0e0e0;
        border-color: #4a5568;
    }

    /* File uploader */
    .stFileUploader {
        background-color: #1a1f2e;
        border: 2px dashed #4a5568;
        border-radius: 0.5rem;
        padding: 1rem;
    }

    /* Data tables */
    .dataframe {
        font-size: 0.9rem;
        background-color: #1a1f2e;
        color: #e0e0e0;
    }

    /* Expander */
    .streamlit-expanderHeader {
        background-color: #1a1f2e;
        color: #e0e0e0;
    }

    /* Info boxes */
    .stInfo {
        background-color: #1a2b3d;
        color: #a8d5ff;
    }

    /* Metric widgets */
    [data-testid="stMetricValue"] {
        color: #4da6ff;
    }

    /* Progress bar */
    .stProgress > div > div > div {
        background-color: #4da6ff;
    }
</style>
""", unsafe_allow_html=True)

# Initialize session state for navigation
if 'navigate_to' not in st.session_state:
    st.session_state['navigate_to'] = None

# Sidebar navigation
st.sidebar.title("🧬 Digital CRO Platform")
st.sidebar.markdown("---")

# Navigation options
nav_options = ["🏠 Home", "🔬 New Screening", "🎯 Multi-Target", "📋 My Jobs", "📚 My Libraries", "📊 Results Browser", "📖 Documentation"]

# Handle navigation from buttons
if st.session_state['navigate_to']:
    default_index = nav_options.index(st.session_state['navigate_to'])
    st.session_state['navigate_to'] = None
else:
    default_index = 0

page = st.sidebar.radio(
    "Navigation",
    nav_options,
    index=default_index
)

st.sidebar.markdown("---")

# Sidebar info
st.sidebar.markdown("### ℹ️ About")
st.sidebar.info(
    "**Digital CRO Platform**\n\n"
    "AI-powered drug discovery and molecular screening service.\n\n"
    "**Version:** 1.0.0\n\n"
    "**Capabilities:**\n"
    "- Molecular docking\n"
    "- ADMET predictions\n"
    "- PDF reports\n"
    "- Batch processing"
)

st.sidebar.markdown("---")

# Quick stats in sidebar
st.sidebar.markdown("### 📊 Quick Stats")
results_dir = Path("data/outputs")
if results_dir.exists():
    completed_runs = len(list(results_dir.glob("*/final_results.csv")))
    st.sidebar.metric("Completed Screenings", completed_runs)
else:
    st.sidebar.metric("Completed Screenings", 0)

# Footer
st.sidebar.markdown("---")
st.sidebar.markdown(
    '<p style="text-align: center; color: #888; font-size: 0.8rem;">'
    'Digital CRO Platform v1.0.0<br>'
    '© 2024 | <a href="mailto:support@digitalcro.com">Support</a>'
    '</p>',
    unsafe_allow_html=True
)

# Main content based on selected page
if page == "🏠 Home":
    from pages import home
    home.show()

elif page == "🔬 New Screening":
    from pages import screening
    screening.show()

elif page == "🎯 Multi-Target":
    from pages import multi_target
    multi_target.show()

elif page == "📋 My Jobs":
    from pages import jobs
    jobs.show()

elif page == "📚 My Libraries":
    from pages import libraries
    libraries.show()

elif page == "📊 Results Browser":
    from pages import results
    results.show()

elif page == "📖 Documentation":
    from pages import documentation
    documentation.show()
