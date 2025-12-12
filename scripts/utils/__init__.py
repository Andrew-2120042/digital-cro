"""
Digital CRO - Utility modules for computational drug discovery.

This package contains utility functions for:
- Molecule operations (loading, conversion, visualization)
- Property calculations (Lipinski rules, molecular descriptors)
- Protein operations (PDB handling, cleaning, analysis)
- Pocket detection (binding site identification, scoring)
- Library filtering (large-scale filtering, PAINS, diversity)
- Ligand/receptor preparation (PDBQT conversion)
- Molecular docking (AutoDock Vina interface)
- ADMET predictions (drug-likeness, BBB penetration, bioavailability)
- Molecular visualization (2D/3D structures, radar plots)
- PDF report generation (professional client-ready reports)
"""

# Phase 1: Molecule operations
from .molecule_ops import (
    load_smiles_file,
    smiles_to_mol,
    draw_molecule,
    draw_molecule_grid
)

from .properties import (
    calculate_lipinski_properties,
    passes_lipinski,
    calculate_all_properties
)

# Phase 2: Protein operations
from .protein_ops import (
    download_pdb,
    load_pdb,
    clean_pdb,
    save_pdb,
    get_protein_info,
    calculate_center_of_mass
)

# Phase 2: Pocket detection
from .pocket_detection import (
    detect_pockets_from_ligand,
    detect_pockets_geometric,
    detect_pockets_auto,
    get_pocket_residues,
    create_docking_box,
    visualize_pocket,
    score_pocket_druggability
)

# Phase 2: Library filtering
from .library_filter import (
    filter_large_library,
    apply_custom_filters,
    detect_pains_substructures,
    detect_reactive_groups,
    calculate_diversity_subset,
    calculate_library_diversity
)

# Phase 3: Ligand preparation
from .ligand_prep import (
    smiles_to_3d_mol,
    add_gasteiger_charges,
    identify_rotatable_bonds,
    write_pdbqt_file,
    prepare_ligand_for_docking,
    batch_prepare_ligands
)

# Phase 3: Receptor preparation
from .receptor_prep import (
    add_polar_hydrogens_simple,
    write_receptor_pdbqt,
    prepare_receptor_for_docking
)

# Phase 3: AutoDock Vina wrapper
from .vina_wrapper import (
    check_vina_installed,
    get_vina_box_from_pocket,
    parse_vina_output,
    dock_single_ligand,
    dock_ligand_with_pocket
)

# Phase 3: Batch docking
from .docking_batch import (
    run_batch_docking,
    dock_ligand_worker,
    create_docking_summary,
    prepare_and_dock_from_smiles
)

# Phase 3: Result analysis
from .result_analysis import (
    rank_docking_results,
    filter_hits_by_score,
    select_diverse_hits,
    cluster_hits_by_similarity,
    generate_hit_report
)

# Phase 4: ADMET predictions
from .admet_predictions import (
    calculate_lipinski_properties as calculate_lipinski_admet,
    predict_admet_properties,
    calculate_sa_score,
    calculate_qed,
    predict_bbb_penetration,
    predict_oral_bioavailability,
    analyze_molecule_complete,
    batch_admet_analysis
)

# Phase 4: Molecular visualization
from .molecular_viz import (
    draw_molecules_grid,
    generate_3d_conformer,
    visualize_docking_results,
    create_admet_radar_plot,
    plot_property_distribution,
    plot_admet_summary
)

# Phase 5: PDF Report Generation
from .report_generator import (
    DrugDiscoveryReport,
    generate_complete_report
)

# Feature 7A: Multi-Target Screening
from .multi_target import (
    run_multi_target_screening,
    create_comparative_table,
    analyze_selectivity
)

__all__ = [
    # Molecule operations
    'load_smiles_file',
    'smiles_to_mol',
    'draw_molecule',
    'draw_molecule_grid',
    # Properties
    'calculate_lipinski_properties',
    'passes_lipinski',
    'calculate_all_properties',
    # Protein operations
    'download_pdb',
    'load_pdb',
    'clean_pdb',
    'save_pdb',
    'get_protein_info',
    'calculate_center_of_mass',
    # Pocket detection
    'detect_pockets_from_ligand',
    'detect_pockets_geometric',
    'detect_pockets_auto',
    'get_pocket_residues',
    'create_docking_box',
    'visualize_pocket',
    'score_pocket_druggability',
    # Library filtering
    'filter_large_library',
    'apply_custom_filters',
    'detect_pains_substructures',
    'detect_reactive_groups',
    'calculate_diversity_subset',
    'calculate_library_diversity',
    # Ligand preparation
    'smiles_to_3d_mol',
    'add_gasteiger_charges',
    'identify_rotatable_bonds',
    'write_pdbqt_file',
    'prepare_ligand_for_docking',
    'batch_prepare_ligands',
    # Receptor preparation
    'add_polar_hydrogens_simple',
    'write_receptor_pdbqt',
    'prepare_receptor_for_docking',
    # Vina wrapper
    'check_vina_installed',
    'get_vina_box_from_pocket',
    'parse_vina_output',
    'dock_single_ligand',
    'dock_ligand_with_pocket',
    # Batch docking
    'run_batch_docking',
    'dock_ligand_worker',
    'create_docking_summary',
    'prepare_and_dock_from_smiles',
    # Result analysis
    'rank_docking_results',
    'filter_hits_by_score',
    'select_diverse_hits',
    'cluster_hits_by_similarity',
    'generate_hit_report',
    # ADMET predictions
    'calculate_lipinski_admet',
    'predict_admet_properties',
    'calculate_sa_score',
    'calculate_qed',
    'predict_bbb_penetration',
    'predict_oral_bioavailability',
    'analyze_molecule_complete',
    'batch_admet_analysis',
    # Molecular visualization
    'draw_molecules_grid',
    'generate_3d_conformer',
    'visualize_docking_results',
    'create_admet_radar_plot',
    'plot_property_distribution',
    'plot_admet_summary',
    # PDF Report Generation
    'DrugDiscoveryReport',
    'generate_complete_report',
    # Multi-Target Screening
    'run_multi_target_screening',
    'create_comparative_table',
    'analyze_selectivity'
]
