"""
AutoDock Vina wrapper module for Digital CRO.

This module provides a Python interface to AutoDock Vina through subprocess calls.

Requirements:
    - AutoDock Vina must be installed and accessible in PATH
    - Install: conda install -c conda-forge autodock-vina
    - Or download from: https://vina.scripps.edu/

Key functions:
    - check_vina_installed(): Verify Vina is available
    - dock_single_ligand(): Run Vina docking
    - parse_vina_output(): Extract binding affinities
    - get_vina_box_from_pocket(): Configure search box
"""

import os
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def check_vina_installed() -> Dict:
    """
    Check if AutoDock Vina is installed and accessible.

    Returns:
        dict: Status dictionary containing:
            - installed (bool): Whether Vina is available
            - version (str): Vina version string
            - path (str): Path to vina executable
            - error (str): Error message if not installed

    Examples:
        >>> status = check_vina_installed()
        >>> if status['installed']:
        ...     print(f"Vina {status['version']} found at {status['path']}")
        ... else:
        ...     print(f"Error: {status['error']}")
    """
    result = {
        'installed': False,
        'version': None,
        'path': None,
        'error': None
    }

    try:
        # Try to run vina --version
        process = subprocess.run(
            ['vina', '--version'],
            capture_output=True,
            text=True,
            timeout=5
        )

        if process.returncode == 0:
            # Parse version from output
            version_match = re.search(r'AutoDock Vina (\d+\.\d+\.\d+)', process.stdout)
            if version_match:
                result['version'] = version_match.group(1)
            else:
                result['version'] = 'unknown'

            # Get path to vina
            which_process = subprocess.run(
                ['which', 'vina'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if which_process.returncode == 0:
                result['path'] = which_process.stdout.strip()

            result['installed'] = True

        else:
            result['error'] = 'Vina command failed'

    except FileNotFoundError:
        result['error'] = 'Vina not found in PATH. Install with: conda install -c conda-forge autodock-vina'
    except subprocess.TimeoutExpired:
        result['error'] = 'Vina command timed out'
    except Exception as e:
        result['error'] = str(e)

    return result


def get_vina_box_from_pocket(pocket: Dict, padding: float = 5.0) -> Dict:
    """
    Convert pocket information to Vina search box parameters.

    Args:
        pocket (dict): Pocket dictionary with 'center' and 'volume' keys
        padding (float): Extra space around pocket (Ångstroms). Default: 5.0

    Returns:
        dict: Vina box parameters:
            - center_x, center_y, center_z (float): Box center coordinates
            - size_x, size_y, size_z (float): Box dimensions

    Examples:
        >>> pocket = {'center': [10.0, 20.0, 30.0], 'volume': 1000.0}
        >>> box = get_vina_box_from_pocket(pocket, padding=5.0)
        >>> print(f"Box center: {box['center_x']}, {box['center_y']}, {box['center_z']}")
    """
    center = pocket['center']
    volume = pocket.get('volume', 1000.0)

    # Estimate box size from volume
    # Assume cubic box: volume = side^3
    side = volume ** (1/3)

    # Add padding
    size_x = size_y = size_z = side + 2 * padding

    # Increase box size by 25% for large flexible molecules
    size_x *= 1.25
    size_y *= 1.25
    size_z *= 1.25

    # Minimum box size (Vina recommendation: at least 22.5 Å to fit ligand)
    min_size = 22.5
    size_x = max(size_x, min_size)
    size_y = max(size_y, min_size)
    size_z = max(size_z, min_size)

    # Maximum box size (increased to ~34 Å for large HIV protease inhibitors)
    max_size = 34.0
    size_x = min(size_x, max_size)
    size_y = min(size_y, max_size)
    size_z = min(size_z, max_size)

    return {
        'center_x': float(center[0]),
        'center_y': float(center[1]),
        'center_z': float(center[2]),
        'size_x': float(size_x),
        'size_y': float(size_y),
        'size_z': float(size_z)
    }


def parse_vina_output(log_text: str) -> List[Dict]:
    """
    Parse Vina output log to extract binding affinities.

    Vina output format:
        mode |   affinity | dist from best mode
             | (kcal/mol) | rmsd l.b.| rmsd u.b.
        -----+------------+----------+----------
           1       -8.2      0.000      0.000
           2       -7.9      1.234      2.567

    Args:
        log_text (str): Vina log output text

    Returns:
        List[dict]: List of docking modes with:
            - mode (int): Mode number
            - affinity (float): Binding affinity (kcal/mol)
            - rmsd_lb (float): RMSD lower bound
            - rmsd_ub (float): RMSD upper bound

    Examples:
        >>> log = "mode |   affinity | dist from best mode\\n   1       -8.2      0.000      0.000"
        >>> modes = parse_vina_output(log)
        >>> print(f"Best affinity: {modes[0]['affinity']} kcal/mol")
    """
    modes = []

    # Find the results table in log
    lines = log_text.split('\n')
    in_results = False

    for line in lines:
        # Start of results table
        if '-----+------------+----------+----------' in line:
            in_results = True
            continue

        # Parse result lines
        if in_results:
            # Match pattern: "   1       -8.2      0.000      0.000"
            match = re.match(r'\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', line)
            if match:
                mode = int(match.group(1))
                affinity = float(match.group(2))
                rmsd_lb = float(match.group(3))
                rmsd_ub = float(match.group(4))

                modes.append({
                    'mode': mode,
                    'affinity': affinity,
                    'rmsd_lb': rmsd_lb,
                    'rmsd_ub': rmsd_ub
                })
            else:
                # End of results table
                if modes:  # Only break if we've found at least one result
                    break

    return modes


def dock_single_ligand(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    output_pdbqt: str,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    exhaustiveness: int = 32,
    num_modes: int = 20,
    energy_range: float = 3.0
) -> Dict:
    """
    Run AutoDock Vina docking for a single ligand.

    Args:
        receptor_pdbqt (str): Path to receptor PDBQT file
        ligand_pdbqt (str): Path to ligand PDBQT file
        output_pdbqt (str): Path to save docked poses
        center (Tuple[float, float, float]): Search box center (x, y, z)
        size (Tuple[float, float, float]): Search box size (x, y, z)
        exhaustiveness (int): Exhaustiveness of search. Default: 32 (higher = more thorough, 4x improvement)
        num_modes (int): Number of binding modes to generate. Default: 20
        energy_range (float): Maximum energy difference (kcal/mol). Default: 3.0

    Returns:
        dict: Docking result containing:
            - success (bool): Whether docking succeeded
            - output_pdbqt (str): Path to docked poses
            - ligand_pdbqt (str): Path to input ligand
            - modes (List[dict]): Binding modes with affinities
            - best_affinity (float): Best binding affinity (kcal/mol)
            - log (str): Full Vina output log
            - error (str): Error message if failed

    Examples:
        >>> result = dock_single_ligand(
        ...     receptor_pdbqt='receptor.pdbqt',
        ...     ligand_pdbqt='ligand.pdbqt',
        ...     output_pdbqt='docked.pdbqt',
        ...     center=(10.0, 20.0, 30.0),
        ...     size=(25.0, 25.0, 25.0)
        ... )
        >>> if result['success']:
        ...     print(f"Best affinity: {result['best_affinity']} kcal/mol")
    """
    result = {
        'success': False,
        'output_pdbqt': output_pdbqt,
        'ligand_pdbqt': ligand_pdbqt,
        'modes': [],
        'best_affinity': None,
        'log': None,
        'error': None
    }

    try:
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_pdbqt) if os.path.dirname(output_pdbqt) else '.', exist_ok=True)

        # Build Vina command
        cmd = [
            'vina',
            '--receptor', receptor_pdbqt,
            '--ligand', ligand_pdbqt,
            '--out', output_pdbqt,
            '--center_x', str(center[0]),
            '--center_y', str(center[1]),
            '--center_z', str(center[2]),
            '--size_x', str(size[0]),
            '--size_y', str(size[1]),
            '--size_z', str(size[2]),
            '--exhaustiveness', str(exhaustiveness),
            '--num_modes', str(num_modes),
            '--energy_range', str(energy_range)
        ]

        # Run Vina
        process = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )

        # Store log
        result['log'] = process.stdout + '\n' + process.stderr

        # Check for errors
        if process.returncode != 0:
            result['error'] = f"Vina failed with return code {process.returncode}"
            return result

        # Parse output
        modes = parse_vina_output(result['log'])

        if not modes:
            result['error'] = 'No binding modes found in Vina output'
            return result

        result['modes'] = modes
        result['best_affinity'] = modes[0]['affinity']
        result['success'] = True

    except subprocess.TimeoutExpired:
        result['error'] = 'Vina docking timed out (>5 minutes)'
    except FileNotFoundError:
        result['error'] = 'Vina not found. Install with: conda install -c conda-forge autodock-vina'
    except Exception as e:
        result['error'] = str(e)

    return result


def dock_ligand_with_pocket(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    output_pdbqt: str,
    pocket: Dict,
    padding: float = 5.0,
    **kwargs
) -> Dict:
    """
    Convenience function: dock ligand using pocket information.

    Args:
        receptor_pdbqt (str): Path to receptor PDBQT file
        ligand_pdbqt (str): Path to ligand PDBQT file
        output_pdbqt (str): Path to save docked poses
        pocket (dict): Pocket dictionary with 'center' and 'volume'
        padding (float): Extra space around pocket. Default: 5.0
        **kwargs: Additional arguments passed to dock_single_ligand()

    Returns:
        dict: Docking result (same as dock_single_ligand)

    Examples:
        >>> pocket = {'center': [10.0, 20.0, 30.0], 'volume': 1000.0}
        >>> result = dock_ligand_with_pocket(
        ...     receptor_pdbqt='receptor.pdbqt',
        ...     ligand_pdbqt='ligand.pdbqt',
        ...     output_pdbqt='docked.pdbqt',
        ...     pocket=pocket
        ... )
    """
    # Get box parameters from pocket
    box = get_vina_box_from_pocket(pocket, padding=padding)

    # Run docking
    return dock_single_ligand(
        receptor_pdbqt=receptor_pdbqt,
        ligand_pdbqt=ligand_pdbqt,
        output_pdbqt=output_pdbqt,
        center=(box['center_x'], box['center_y'], box['center_z']),
        size=(box['size_x'], box['size_y'], box['size_z']),
        **kwargs
    )
