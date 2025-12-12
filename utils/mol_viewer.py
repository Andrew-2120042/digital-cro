"""
3D Molecular Visualization Module

Utilities for displaying:
* Protein structures
* Docked ligand poses
* Protein-ligand complexes
* Binding interactions
"""

import py3Dmol
from pathlib import Path
from typing import Optional, List, Dict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def create_protein_viewer(
    pdb_path: str,
    width: int = 800,
    height: int = 600,
    style: str = 'cartoon'
) -> py3Dmol.view:
    """
    Create 3D viewer for protein structure.

    Args:
        pdb_path: Path to PDB file
        width: Viewer width in pixels
        height: Viewer height in pixels
        style: Protein representation ('cartoon', 'stick', 'sphere', 'line')

    Returns:
        py3Dmol viewer object
    """
    try:
        # Read PDB file
        with open(pdb_path, 'r') as f:
            pdb_data = f.read()

        # Create viewer
        view = py3Dmol.view(width=width, height=height)
        view.addModel(pdb_data, 'pdb')

        # Set style
        if style == 'cartoon':
            view.setStyle({'cartoon': {'color': 'spectrum'}})
        elif style == 'stick':
            view.setStyle({'stick': {}})
        elif style == 'sphere':
            view.setStyle({'sphere': {}})
        else:
            view.setStyle({'line': {}})

        # Center and zoom
        view.zoomTo()

        logger.info(f"Created protein viewer for {pdb_path}")

        return view

    except Exception as e:
        logger.error(f"Error creating protein viewer: {e}")
        raise


def create_ligand_viewer(
    pdbqt_path: str,
    width: int = 400,
    height: int = 400
) -> py3Dmol.view:
    """
    Create 3D viewer for ligand structure.

    Args:
        pdbqt_path: Path to PDBQT file
        width: Viewer width in pixels
        height: Viewer height in pixels

    Returns:
        py3Dmol viewer object
    """
    try:
        # Read PDBQT file
        with open(pdbqt_path, 'r') as f:
            pdbqt_data = f.read()

        # Create viewer
        view = py3Dmol.view(width=width, height=height)
        view.addModel(pdbqt_data, 'pdbqt')

        # Set style - stick with colors
        view.setStyle({'stick': {'colorscheme': 'default'}})

        # Center and zoom
        view.zoomTo()

        return view

    except Exception as e:
        logger.error(f"Error creating ligand viewer: {e}")
        raise


def create_complex_viewer(
    protein_pdb_path: str,
    ligand_pdbqt_path: str,
    pocket_center: tuple = None,
    pocket_size: tuple = None,
    width: int = 800,
    height: int = 600,
    show_pocket: bool = True
) -> py3Dmol.view:
    """
    Create 3D viewer for protein-ligand complex.

    This is the KILLER FEATURE - shows actual docking results.

    Args:
        protein_pdb_path: Path to protein PDB file
        ligand_pdbqt_path: Path to docked ligand PDBQT file
        pocket_center: (x, y, z) coordinates of binding pocket
        pocket_size: (x, y, z) dimensions of binding box
        width: Viewer width in pixels
        height: Viewer height in pixels
        show_pocket: Whether to highlight binding pocket

    Returns:
        py3Dmol viewer object
    """
    try:
        # Read files
        with open(protein_pdb_path, 'r') as f:
            protein_data = f.read()

        with open(ligand_pdbqt_path, 'r') as f:
            ligand_data = f.read()

        # Create viewer
        view = py3Dmol.view(width=width, height=height)

        # Add protein
        view.addModel(protein_data, 'pdb')

        # Style protein - cartoon with semi-transparent surface
        view.setStyle(
            {'model': 0},
            {
                'cartoon': {'color': 'lightgray', 'opacity': 0.7},
            }
        )

        # Add ligand
        view.addModel(ligand_data, 'pdbqt')

        # Style ligand - colorful sticks
        view.setStyle(
            {'model': 1},
            {
                'stick': {
                    'colorscheme': 'greenCarbon',
                    'radius': 0.15
                }
            }
        )

        # Highlight binding pocket if coordinates provided
        if show_pocket and pocket_center:
            # Get pocket residues (within 5Å of pocket center)
            x, y, z = pocket_center

            # Highlight pocket residues
            view.addStyle(
                {
                    'model': 0,
                    'within': {
                        'distance': 8,
                        'sel': {'model': 1}
                    }
                },
                {
                    'stick': {
                        'colorscheme': 'orangeCarbon',
                        'radius': 0.2
                    }
                }
            )

            # Add pocket box (optional)
            if pocket_size:
                sx, sy, sz = pocket_size
                view.addBox({
                    'center': {'x': x, 'y': y, 'z': z},
                    'dimensions': {'w': sx, 'h': sy, 'd': sz},
                    'color': 'blue',
                    'opacity': 0.2
                })

        # Center on ligand
        view.zoomTo({'model': 1})

        logger.info(f"Created complex viewer: {protein_pdb_path} + {ligand_pdbqt_path}")

        return view

    except Exception as e:
        logger.error(f"Error creating complex viewer: {e}")
        raise


def create_multi_pose_viewer(
    protein_pdb_path: str,
    ligand_pdbqt_paths: List[str],
    width: int = 1000,
    height: int = 600,
    max_poses: int = 9
) -> py3Dmol.view:
    """
    Create viewer showing multiple binding poses for one molecule.

    Args:
        protein_pdb_path: Path to protein PDB file
        ligand_pdbqt_paths: List of paths to different poses (PDBQT files)
        width: Viewer width in pixels
        height: Viewer height in pixels
        max_poses: Maximum number of poses to show

    Returns:
        py3Dmol viewer object with all poses overlaid
    """
    try:
        # Read protein
        with open(protein_pdb_path, 'r') as f:
            protein_data = f.read()

        # Create viewer
        view = py3Dmol.view(width=width, height=height)

        # Add protein
        view.addModel(protein_data, 'pdb')
        view.setStyle(
            {'model': 0},
            {'cartoon': {'color': 'lightgray', 'opacity': 0.5}}
        )

        # Color palette for different poses
        colors = [
            'greenCarbon', 'cyanCarbon', 'magentaCarbon',
            'yellowCarbon', 'orangeCarbon', 'purpleCarbon',
            'blueCarbon', 'redCarbon', 'whiteCarbon'
        ]

        # Add ligand poses
        for idx, pose_path in enumerate(ligand_pdbqt_paths[:max_poses]):
            with open(pose_path, 'r') as f:
                ligand_data = f.read()

            view.addModel(ligand_data, 'pdbqt')

            # Style each pose with different color
            view.setStyle(
                {'model': idx + 1},
                {
                    'stick': {
                        'colorscheme': colors[idx % len(colors)],
                        'radius': 0.15,
                        'opacity': 0.7
                    }
                }
            )

        # Center view
        view.zoomTo()

        logger.info(f"Created multi-pose viewer with {len(ligand_pdbqt_paths)} poses")

        return view

    except Exception as e:
        logger.error(f"Error creating multi-pose viewer: {e}")
        raise


def add_distance_labels(
    view: py3Dmol.view,
    atom_pairs: List[Dict]
) -> py3Dmol.view:
    """
    Add distance measurements between atoms.

    Args:
        view: py3Dmol viewer object
        atom_pairs: List of dicts with keys:
            - 'atom1': {'serial': int, 'model': int}
            - 'atom2': {'serial': int, 'model': int}
            - 'label': str (optional)

    Returns:
        Updated viewer object
    """
    for pair in atom_pairs:
        view.addLabel(
            pair.get('label', 'Distance'),
            {
                'position': pair['atom1'],
                'backgroundColor': 'white',
                'fontColor': 'black'
            }
        )

        # Draw line between atoms
        view.addCylinder({
            'start': pair['atom1'],
            'end': pair['atom2'],
            'radius': 0.1,
            'color': 'yellow',
            'dashed': True
        })

    return view


def get_viewer_html(view: py3Dmol.view) -> str:
    """
    Get HTML representation of viewer for embedding.

    Args:
        view: py3Dmol viewer object

    Returns:
        HTML string
    """
    return view._make_html()
