"""
Binding pocket detection module for Digital CRO.

This module implements multiple methods for detecting and analyzing
binding pockets in protein structures:
- Ligand-based detection (most accurate when ligands present)
- Geometric grid-based detection (fallback method)
- Auto-detection (intelligent method selection)
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from scipy.spatial.distance import cdist
from scipy.cluster.hierarchy import fclusterdata
from Bio.PDB import Structure
import warnings


# ============================================================================
# LIGAND-BASED POCKET DETECTION (Primary Method)
# ============================================================================

def is_relevant_ligand(residue) -> bool:
    """
    Determine if heteroatom is a relevant ligand (not water/ion/buffer).

    Excludes:
        - HOH, WAT, H2O (water)
        - SO4, PO4, PO2, POO (phosphate/sulfate buffers)
        - GOL, EDO, PGE (cryoprotectants)
        - ACT, FMT, DMS (common crystallization additives)
        - Single-atom residues (metal ions)

    Includes:
        - Organic molecules with >3 heavy atoms
        - Known drug-like ligands

    Args:
        residue: Biopython Residue object

    Returns:
        bool: True if likely a relevant ligand

    Examples:
        >>> for residue in chain:
        ...     if is_relevant_ligand(residue):
        ...         print(f"Found ligand: {residue.resname}")
    """
    # Exclude common non-ligand heteroatoms
    excluded = {
        'HOH', 'WAT', 'H2O',  # Water
        'SO4', 'PO4', 'PO2', 'POO', 'SUL',  # Buffers
        'GOL', 'EDO', 'PGE', 'PG4', 'PEG',  # Cryoprotectants
        'ACT', 'FMT', 'DMS', 'BME', 'ACE',  # Additives
        'CL', 'NA', 'MG', 'CA', 'ZN', 'FE', 'K', 'MN',  # Ions
        'CD', 'CO', 'CU', 'NI'  # Metal ions
    }

    resname = residue.resname.strip()

    if resname in excluded:
        return False

    # Count heavy atoms (non-hydrogen)
    heavy_atoms = [atom for atom in residue.get_atoms() if atom.element != 'H']

    # Ligands should have at least 4 heavy atoms
    if len(heavy_atoms) < 4:
        return False

    return True


def calculate_ligand_center(residue) -> Tuple[float, float, float]:
    """
    Calculate geometric center of ligand atoms.

    Args:
        residue: Biopython Residue object

    Returns:
        Tuple[float, float, float]: (x, y, z) coordinates in Angstroms

    Examples:
        >>> center = calculate_ligand_center(ligand_residue)
        >>> print(f"Ligand center: {center}")
    """
    coords = np.array([atom.coord for atom in residue.get_atoms()])
    center = coords.mean(axis=0)
    return tuple(center)


def estimate_pocket_volume_from_ligand(residue) -> float:
    """
    Estimate pocket volume based on ligand size.

    Uses a simple approximation:
    - Calculate bounding box of ligand atoms
    - Estimate pocket as 1.5x ligand volume (pocket is larger than ligand)

    Args:
        residue: Biopython Residue object

    Returns:
        float: Estimated volume in Ų

    Notes:
        Typical drug-like pockets: 200-2000 Ų
        Small molecule ligands: 100-500 Ų
        Pocket volume ≈ 1.5-2x ligand volume

    Examples:
        >>> volume = estimate_pocket_volume_from_ligand(ligand)
        >>> print(f"Estimated pocket volume: {volume:.1f} Ų")
    """
    coords = np.array([atom.coord for atom in residue.get_atoms()])

    if len(coords) == 0:
        return 0.0

    # Calculate bounding box
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    dimensions = max_coords - min_coords

    # Bounding box volume
    bbox_volume = np.prod(dimensions)

    # Estimate pocket as ~1.5x ligand bounding box
    # (pocket needs to accommodate ligand plus some space)
    estimated_pocket = bbox_volume * 1.5

    return float(estimated_pocket)


def detect_pockets_from_ligand(
    structure: Structure.Structure,
    ligand_names: Optional[List[str]] = None
) -> List[Dict]:
    """
    Detect binding pocket based on co-crystallized ligand location.

    This is the MOST RELIABLE method - if a ligand is present, its
    location defines the binding site by definition.

    Args:
        structure (Structure): Biopython Structure object
        ligand_names (List[str], optional): Specific heteroatom IDs to search for.
            If None, finds all relevant heteroatoms (excluding water/ions).
            Common examples: ['MK1', 'JE2', 'N3', 'LIG', 'UNK', 'INH']

    Returns:
        List[Dict]: List of pockets, each containing:
            - center (Tuple[float, float, float]): Ligand center coordinates
            - volume (float): Estimated pocket volume in Ų
            - score (float): 1.0 (highest - this is a known site)
            - residues (List[str]): Nearby residues within 10Å
            - ligand_name (str): Heteroatom residue name
            - num_atoms (int): Number of ligand atoms
            - method (str): 'ligand-based'

    Algorithm:
        1. Scan all chains for HETATM records
        2. Filter to relevant ligands (exclude water, ions, buffers)
        3. For each ligand:
           - Calculate geometric center
           - Estimate pocket volume from ligand size
           - Find residues within 10Å
        4. Return pockets sorted by ligand size (larger = more significant)

    Examples:
        >>> pockets = detect_pockets_from_ligand(structure)
        >>> if pockets:
        ...     print(f"Found binding site at {pockets[0]['center']}")
        ...     print(f"Ligand: {pockets[0]['ligand_name']}")

    Notes:
        - Returns empty list if no ligands found
        - Multiple ligands will generate multiple pockets
        - Accuracy: Typically <1Å from true binding site
    """
    pockets = []
    model = structure[0]

    for chain in model:
        for residue in chain:
            hetero_flag = residue.id[0]

            # Only process heteroatoms (not standard amino acids)
            if hetero_flag == ' ':
                continue

            resname = residue.resname.strip()

            # If specific ligand names provided, filter to those
            if ligand_names is not None:
                if resname not in ligand_names:
                    continue
            else:
                # Otherwise, use intelligent filtering
                if not is_relevant_ligand(residue):
                    continue

            # Calculate ligand properties
            center = calculate_ligand_center(residue)
            volume = estimate_pocket_volume_from_ligand(residue)
            num_atoms = len(list(residue.get_atoms()))

            # Get nearby residues
            residues = get_pocket_residues(structure, center, radius=10.0)

            # Create pocket entry
            pocket = {
                'center': center,
                'volume': volume,
                'score': 1.0,  # Maximum score - known binding site
                'residues': residues,
                'ligand_name': resname,
                'num_atoms': num_atoms,
                'method': 'ligand-based'
            }

            pockets.append(pocket)

    # Sort by ligand size (larger ligands = more significant pockets)
    pockets.sort(key=lambda p: p['num_atoms'], reverse=True)

    return pockets


def detect_pockets_auto(
    structure: Structure.Structure,
    prefer_ligand: bool = True,
    **kwargs
) -> List[Dict]:
    """
    MAIN FUNCTION: Automatically detect pockets using best available method.

    Decision tree:
    1. Check for co-crystallized ligands
    2. If ligands found AND prefer_ligand=True → use ligand-based (highest accuracy)
    3. If no ligands → fall back to geometric method
    4. Apply volume filtering to geometric results (200-2000 Ų)
    5. Return best pockets (max 5)

    Args:
        structure (Structure): Biopython Structure object
        prefer_ligand (bool): If True, prefer ligand method when available. Default: True
        **kwargs: Additional arguments passed to geometric detection if used

    Returns:
        List[Dict]: Detected pockets with 'method' field indicating detection method:
            - 'ligand-based': From co-crystallized ligand (most reliable)
            - 'geometric': From grid-based algorithm (fallback)

    This is the recommended function for client code to call.

    Examples:
        >>> # Automatic detection
        >>> pockets = detect_pockets_auto(structure)
        >>> print(f"Method: {pockets[0]['method']}")
        >>> print(f"Accuracy: {'High' if pockets[0]['method'] == 'ligand-based' else 'Moderate'}")

        >>> # Force geometric method
        >>> pockets = detect_pockets_auto(structure, prefer_ligand=False)

    Performance:
        - Ligand-based: <0.1 seconds, <1Å accuracy
        - Geometric: 3-5 seconds, ~5-10Å accuracy
    """
    # Try ligand-based detection first
    if prefer_ligand:
        ligand_pockets = detect_pockets_from_ligand(structure)

        if ligand_pockets:
            # Found ligands - use those (most reliable)
            return ligand_pockets[:5]  # Return up to 5

    # No ligands or prefer_ligand=False - use geometric method
    # Apply strict volume filtering for druggable pockets
    geometric_pockets = detect_pockets_geometric(
        structure,
        filter_by_volume=True,
        min_volume=200.0,
        max_volume=2000.0,
        **kwargs
    )

    # Add method tag
    for pocket in geometric_pockets:
        pocket['method'] = 'geometric'

    return geometric_pockets[:5]


def detect_pockets_geometric(
    structure: Structure.Structure,
    grid_spacing: float = 1.0,
    probe_radius: float = 1.4,
    min_volume: float = 100.0,
    max_pockets: int = 10,
    filter_by_volume: bool = False,
    max_volume: float = 2000.0
) -> List[Dict]:
    """
    Detect pockets using geometric grid-based method.

    Algorithm:
    1. Create 3D grid around protein
    2. For each grid point, check if it's in a cavity:
       - Must be at least probe_radius away from any atom (accessible)
       - Must be within detection range of protein surface
    3. Cluster nearby cavity points
    4. Calculate volume and score for each pocket
    5. Rank pockets by druggability score
    6. Optionally filter by volume range (druggable pockets: 200-2000 Ų)

    Args:
        structure (Structure): Biopython Structure object
        grid_spacing (float): Distance between grid points (Angstroms). Default: 1.0
        probe_radius (float): Probe sphere radius (Angstroms). Default: 1.4 (water)
        min_volume (float): Minimum pocket volume (Ų). Default: 100.0
        max_pockets (int): Maximum number of pockets to return. Default: 10
        filter_by_volume (bool): If True, only return pockets in druggable range. Default: False
        max_volume (float): Maximum pocket volume (Ų). Default: 2000.0

    Returns:
        List[Dict]: List of pockets, each containing:
            - center (Tuple[float, float, float]): Pocket center coordinates
            - volume (float): Estimated volume in ų
            - score (float): Druggability score (0-1)
            - residues (List[str]): Nearby residues
            - num_points (int): Number of grid points in pocket
            - method (str): 'geometric' (added by detect_pockets_auto)

    Examples:
        >>> # Basic detection
        >>> pockets = detect_pockets_geometric(structure)
        >>> best = pockets[0]
        >>> print(f"Center: {best['center']}, Score: {best['score']:.2f}")

        >>> # Filter to druggable pockets only
        >>> pockets = detect_pockets_geometric(structure, filter_by_volume=True)
        >>> print(f"All pockets are between 200-2000 Ų")
    """
    # Get all atom coordinates from first model
    model = structure[0]
    atoms = list(model.get_atoms())

    if len(atoms) == 0:
        raise ValueError("Structure contains no atoms")

    atom_coords = np.array([atom.coord for atom in atoms])

    # Calculate bounding box with padding
    padding = 5.0  # Angstroms
    min_coords = atom_coords.min(axis=0) - padding
    max_coords = atom_coords.max(axis=0) + padding

    # Create 3D grid
    x_range = np.arange(min_coords[0], max_coords[0], grid_spacing)
    y_range = np.arange(min_coords[1], max_coords[1], grid_spacing)
    z_range = np.arange(min_coords[2], max_coords[2], grid_spacing)

    # Generate grid points
    grid_points = []
    for x in x_range:
        for y in y_range:
            for z in z_range:
                grid_points.append([x, y, z])

    grid_points = np.array(grid_points)

    if len(grid_points) == 0:
        return []

    # Calculate distances from grid points to all atoms
    distances = cdist(grid_points, atom_coords)
    min_distances = distances.min(axis=1)

    # Identify cavity points
    # Point is in cavity if:
    # 1. At least probe_radius away from any atom (accessible)
    # 2. Within detection range (not too far from protein)
    max_distance = 8.0  # Maximum distance from protein surface
    cavity_mask = (min_distances >= probe_radius) & (min_distances <= max_distance)
    cavity_points = grid_points[cavity_mask]

    if len(cavity_points) < 10:
        warnings.warn("Very few cavity points detected. Protein may be too small or compact.")
        return []

    # Cluster cavity points to identify distinct pockets
    # Use hierarchical clustering with distance threshold
    cluster_distance = 3.0  # Points within 3Å are in same pocket
    try:
        clusters = fclusterdata(cavity_points, cluster_distance, criterion='distance')
    except Exception as e:
        warnings.warn(f"Clustering failed: {e}")
        return []

    # Analyze each cluster
    pockets = []
    unique_clusters = np.unique(clusters)

    for cluster_id in unique_clusters:
        cluster_mask = clusters == cluster_id
        cluster_pts = cavity_points[cluster_mask]

        # Calculate pocket volume (approximate)
        volume = len(cluster_pts) * (grid_spacing ** 3)

        # Filter by minimum volume
        if volume < min_volume:
            continue

        # Calculate pocket center
        center = cluster_pts.mean(axis=0)

        # Get nearby residues
        residues = get_pocket_residues(structure, tuple(center), radius=10.0)

        # Calculate druggability score
        pocket_data = {
            'center': tuple(center),
            'volume': volume,
            'num_points': len(cluster_pts),
            'residues': residues,
            'grid_points': cluster_pts  # For advanced analysis
        }

        score = score_pocket_druggability(pocket_data, structure)
        pocket_data['score'] = score

        # Remove grid_points from output (too large)
        del pocket_data['grid_points']

        pockets.append(pocket_data)

    # Sort by score (best first)
    pockets.sort(key=lambda p: p['score'], reverse=True)

    # Apply volume filtering if requested
    if filter_by_volume:
        # Filter to druggable size range
        filtered_pockets = [
            p for p in pockets
            if min_volume <= p['volume'] <= max_volume
        ]

        if not filtered_pockets:
            # If no pockets in range, return top 3 anyway but warn
            warnings.warn(
                f"No pockets in volume range {min_volume}-{max_volume} Ų. "
                f"Returning top {min(3, len(pockets))} pockets anyway."
            )
            if pockets:
                print(f"  Volume range detected: {pockets[0]['volume']:.1f} - {pockets[-1]['volume']:.1f} Ų")
            return pockets[:3]

        return filtered_pockets[:max_pockets]

    # Return top pockets without filtering
    return pockets[:max_pockets]


def get_pocket_residues(
    structure: Structure.Structure,
    center: Tuple[float, float, float],
    radius: float = 10.0
) -> List[str]:
    """
    Get all residues within radius of pocket center.

    Args:
        structure (Structure): Biopython Structure object
        center (Tuple[float, float, float]): Pocket center coordinates
        radius (float): Distance in Angstroms. Default: 10.0

    Returns:
        List[str]: Residue identifiers (e.g., ['ALA_A_23', 'TRP_A_45'])

    Examples:
        >>> residues = get_pocket_residues(structure, (2.0, 16.0, 25.0))
        >>> print(residues[:3])
        ['PRO_A_1', 'GLN_A_2', 'ILE_A_3']
    """
    model = structure[0]
    center = np.array(center)
    nearby_residues = []

    for chain in model:
        for residue in chain:
            # Skip non-protein residues
            if residue.id[0] != ' ':
                continue

            # Calculate distance from residue to pocket center
            # Use CA atom (alpha carbon) as residue position
            if 'CA' in residue:
                ca_coord = residue['CA'].coord
                distance = np.linalg.norm(ca_coord - center)

                if distance <= radius:
                    res_name = residue.resname
                    chain_id = chain.id
                    res_num = residue.id[1]
                    residue_id = f"{res_name}_{chain_id}_{res_num}"
                    nearby_residues.append(residue_id)

    return nearby_residues


def create_docking_box(
    pocket_center: Tuple[float, float, float],
    box_size: float = 20.0,
    size_x: Optional[float] = None,
    size_y: Optional[float] = None,
    size_z: Optional[float] = None
) -> Dict[str, float]:
    """
    Create AutoDock Vina search box parameters around pocket.

    Args:
        pocket_center (Tuple[float, float, float]): Pocket center coordinates
        box_size (float): Default edge length of cubic box (Angstroms). Default: 20.0
        size_x (float, optional): Custom X dimension
        size_y (float, optional): Custom Y dimension
        size_z (float, optional): Custom Z dimension

    Returns:
        dict: Vina docking box parameters:
            - center_x, center_y, center_z: Box center coordinates
            - size_x, size_y, size_z: Box dimensions

    Examples:
        >>> box = create_docking_box((2.0, 16.0, 25.0), box_size=22.5)
        >>> print(box)
        {'center_x': 2.0, 'center_y': 16.0, 'center_z': 25.0,
         'size_x': 22.5, 'size_y': 22.5, 'size_z': 22.5}
    """
    box_params = {
        'center_x': float(pocket_center[0]),
        'center_y': float(pocket_center[1]),
        'center_z': float(pocket_center[2]),
        'size_x': float(size_x if size_x is not None else box_size),
        'size_y': float(size_y if size_y is not None else box_size),
        'size_z': float(size_z if size_z is not None else box_size)
    }

    return box_params


def visualize_pocket(
    structure: Structure.Structure,
    pocket: Dict,
    output_path: Optional[str] = None,
    show_residues: bool = False
) -> None:
    """
    Create 3D visualization of pocket location on protein.

    Args:
        structure (Structure): Biopython Structure object
        pocket (Dict): Pocket data from detect_pockets_geometric
        output_path (str, optional): Path to save PNG. If None, displays plot.
        show_residues (bool): Label nearby residues. Default: False

    Generates:
        - Matplotlib 3D scatter plot
        - Protein CA atoms (gray points)
        - Pocket center (red star)
        - Pocket sphere (transparent red)

    Examples:
        >>> visualize_pocket(structure, pocket, 'pocket_viz.png')
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # Extract protein CA coordinates
    model = structure[0]
    ca_coords = []

    for chain in model:
        for residue in chain:
            if residue.id[0] == ' ' and 'CA' in residue:
                ca_coords.append(residue['CA'].coord)

    ca_coords = np.array(ca_coords)

    # Create figure
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot protein backbone (CA atoms)
    ax.scatter(
        ca_coords[:, 0],
        ca_coords[:, 1],
        ca_coords[:, 2],
        c='gray',
        alpha=0.3,
        s=20,
        label='Protein (CA atoms)'
    )

    # Plot pocket center
    center = pocket['center']
    ax.scatter(
        center[0],
        center[1],
        center[2],
        c='red',
        marker='*',
        s=500,
        label=f"Pocket (score={pocket['score']:.2f})",
        edgecolors='black',
        linewidths=2
    )

    # Draw pocket sphere
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    radius = (pocket['volume'] / (4/3 * np.pi)) ** (1/3)  # Approximate radius

    x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z, color='red', alpha=0.2)

    # Labels and formatting
    ax.set_xlabel('X (Å)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y (Å)', fontsize=12, fontweight='bold')
    ax.set_zlabel('Z (Å)', fontsize=12, fontweight='bold')
    ax.set_title(
        f'Binding Pocket Visualization\n'
        f'Volume: {pocket["volume"]:.1f} Ų, '
        f'Druggability: {pocket["score"]:.2f}',
        fontsize=14,
        fontweight='bold'
    )
    ax.legend(fontsize=10)

    # Equal aspect ratio
    max_range = np.array([
        ca_coords[:, 0].max() - ca_coords[:, 0].min(),
        ca_coords[:, 1].max() - ca_coords[:, 1].min(),
        ca_coords[:, 2].max() - ca_coords[:, 2].min()
    ]).max() / 2.0

    mid_x = (ca_coords[:, 0].max() + ca_coords[:, 0].min()) * 0.5
    mid_y = (ca_coords[:, 1].max() + ca_coords[:, 1].min()) * 0.5
    mid_z = (ca_coords[:, 2].max() + ca_coords[:, 2].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def score_pocket_druggability(pocket_data: Dict, structure: Structure.Structure) -> float:
    """
    Score pocket druggability (0-1 scale).

    Factors considered:
    - Volume: 200-2000 Ų optimal (weight: 0.4)
    - Depth: Measured by distance from surface (weight: 0.3)
    - Shape: Sphericity vs flatness (weight: 0.3)

    Args:
        pocket_data (Dict): Pocket data with 'volume', 'center', 'grid_points'
        structure (Structure): Biopython Structure object

    Returns:
        float: Druggability score (0.0-1.0, higher is better)
              > 0.7: Highly druggable
              0.5-0.7: Moderately druggable
              < 0.5: Poorly druggable

    Examples:
        >>> score = score_pocket_druggability(pocket_data, structure)
        >>> if score > 0.7:
        ...     print("Highly druggable pocket")
    """
    score = 0.0

    # 1. Volume score (40% weight)
    # Optimal volume range: 200-2000 Ų
    volume = pocket_data['volume']
    if 200 <= volume <= 2000:
        volume_score = 1.0
    elif volume < 200:
        volume_score = max(0, volume / 200)
    else:  # volume > 2000
        volume_score = max(0, 1.0 - (volume - 2000) / 2000)

    score += 0.4 * volume_score

    # 2. Depth score (30% weight)
    # Deeper pockets are more druggable
    # Estimate depth by looking at grid points
    if 'grid_points' in pocket_data:
        grid_points = pocket_data['grid_points']

        # Get protein surface (approximate)
        model = structure[0]
        atom_coords = np.array([atom.coord for atom in model.get_atoms()])

        # Calculate distances from pocket points to nearest atom
        if len(grid_points) > 0 and len(atom_coords) > 0:
            distances = cdist(grid_points, atom_coords)
            avg_depth = distances.min(axis=1).mean()

            # Optimal depth: 3-6 Å from surface
            if 3.0 <= avg_depth <= 6.0:
                depth_score = 1.0
            elif avg_depth < 3.0:
                depth_score = avg_depth / 3.0
            else:
                depth_score = max(0, 1.0 - (avg_depth - 6.0) / 4.0)

            score += 0.3 * depth_score
        else:
            score += 0.15  # Default partial score
    else:
        score += 0.15  # Default partial score if no grid points

    # 3. Shape score (30% weight)
    # More spherical pockets score higher
    if 'grid_points' in pocket_data:
        grid_points = pocket_data['grid_points']

        if len(grid_points) >= 10:
            # Calculate variance in each dimension
            variances = grid_points.var(axis=0)
            max_var = variances.max()
            min_var = variances.min()

            # Sphericity: ratio of min to max variance (1.0 = perfect sphere)
            if max_var > 0:
                sphericity = min_var / max_var
                shape_score = sphericity
            else:
                shape_score = 0.5

            score += 0.3 * shape_score
        else:
            score += 0.15  # Default partial score
    else:
        score += 0.15  # Default partial score

    return min(1.0, score)


def find_known_binding_site(
    structure: Structure.Structure,
    ligand_name: str = 'MK1'
) -> Optional[Tuple[float, float, float]]:
    """
    Find known binding site by locating co-crystallized ligand.

    Useful for validation - if PDB has a ligand, find its center.

    Args:
        structure (Structure): Biopython Structure object
        ligand_name (str): Ligand residue name (e.g., 'MK1', 'ATP')

    Returns:
        Optional[Tuple[float, float, float]]: Ligand center coordinates, or None if not found

    Examples:
        >>> known_site = find_known_binding_site(structure, 'MK1')
        >>> if known_site:
        ...     print(f"Known ligand at {known_site}")
    """
    model = structure[0]
    ligand_atoms = []

    for chain in model:
        for residue in chain:
            if residue.resname == ligand_name:
                ligand_atoms.extend([atom.coord for atom in residue.get_atoms()])

    if ligand_atoms:
        ligand_coords = np.array(ligand_atoms)
        center = ligand_coords.mean(axis=0)
        return tuple(center)

    return None
