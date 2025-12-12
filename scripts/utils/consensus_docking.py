"""
Consensus Docking Module

Run multiple docking programs and combine results for higher confidence.

Supported methods:
- AutoDock Vina (default)
- Smina (Vina fork with better scoring)
- LeDock (optional, if installed)
"""

from pathlib import Path
from typing import List, Dict, Tuple, Optional
import logging
import subprocess
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ConsensusDocking:
    """
    Multi-method consensus docking for higher confidence results.
    """

    def __init__(self):
        self.available_methods = self._check_available_methods()
        logger.info(f"Available docking methods: {self.available_methods}")

    def _check_available_methods(self) -> List[str]:
        """Check which docking programs are installed"""
        methods = ['vina']  # Vina always available (we installed it)

        # Check for Smina
        if shutil.which('smina'):
            methods.append('smina')

        # Check for LeDock
        if shutil.which('ledock'):
            methods.append('ledock')

        return methods

    def dock_with_vina(
        self,
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float],
        output_pdbqt: str,
        exhaustiveness: int = 8
    ) -> float:
        """Run AutoDock Vina docking"""

        cmd = [
            'vina',
            '--receptor', receptor_pdbqt,
            '--ligand', ligand_pdbqt,
            '--center_x', str(center[0]),
            '--center_y', str(center[1]),
            '--center_z', str(center[2]),
            '--size_x', str(size[0]),
            '--size_y', str(size[1]),
            '--size_z', str(size[2]),
            '--exhaustiveness', str(exhaustiveness),
            '--out', output_pdbqt
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )

            # Parse affinity from output
            for line in result.stdout.split('\n'):
                if 'REMARK VINA RESULT:' in line or line.strip().startswith('1'):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            return float(parts[1])
                        except (ValueError, IndexError):
                            continue

            logger.warning("Could not parse Vina output")
            return 0.0

        except Exception as e:
            logger.error(f"Vina docking failed: {e}")
            return 0.0

    def dock_with_smina(
        self,
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float],
        output_pdbqt: str,
        exhaustiveness: int = 8
    ) -> float:
        """Run Smina docking (Vina fork with improved scoring)"""

        if 'smina' not in self.available_methods:
            logger.warning("Smina not installed, skipping")
            return None

        cmd = [
            'smina',
            '--receptor', receptor_pdbqt,
            '--ligand', ligand_pdbqt,
            '--center_x', str(center[0]),
            '--center_y', str(center[1]),
            '--center_z', str(center[2]),
            '--size_x', str(size[0]),
            '--size_y', str(size[1]),
            '--size_z', str(size[2]),
            '--exhaustiveness', str(exhaustiveness),
            '--out', output_pdbqt
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )

            # Parse affinity
            for line in result.stdout.split('\n'):
                if line.strip().startswith('1'):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            return float(parts[1])
                        except (ValueError, IndexError):
                            continue

            return 0.0

        except Exception as e:
            logger.error(f"Smina docking failed: {e}")
            return None

    def consensus_dock(
        self,
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        ligand_id: str,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float],
        output_dir: str,
        methods: Optional[List[str]] = None
    ) -> Dict:
        """
        Run consensus docking with multiple methods.

        Args:
            receptor_pdbqt: Receptor PDBQT file
            ligand_pdbqt: Ligand PDBQT file
            ligand_id: Ligand identifier
            center: Docking box center (x, y, z)
            size: Docking box size (x, y, z)
            output_dir: Output directory
            methods: List of methods to use (default: all available)

        Returns:
            Dict with scores from each method and consensus metrics
        """
        if methods is None:
            methods = self.available_methods

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        results = {
            'ligand_id': ligand_id,
            'methods': {}
        }

        # Run each method
        for method in methods:
            output_pdbqt = output_dir / f"{ligand_id}_{method}.pdbqt"

            if method == 'vina':
                score = self.dock_with_vina(
                    receptor_pdbqt, ligand_pdbqt,
                    center, size, str(output_pdbqt)
                )
            elif method == 'smina':
                score = self.dock_with_smina(
                    receptor_pdbqt, ligand_pdbqt,
                    center, size, str(output_pdbqt)
                )
            else:
                score = None

            if score is not None:
                results['methods'][method] = {
                    'score': score,
                    'output': str(output_pdbqt)
                }

        # Calculate consensus metrics
        scores = [v['score'] for v in results['methods'].values()]

        if len(scores) >= 2:
            results['consensus_score'] = np.mean(scores)
            results['score_std'] = np.std(scores)
            results['score_range'] = max(scores) - min(scores)
            results['agreement'] = results['score_std'] < 2.0  # Good agreement if stdev < 2
            results['confidence'] = 'high' if results['agreement'] else 'low'
        else:
            results['consensus_score'] = scores[0] if scores else 0.0
            results['score_std'] = 0.0
            results['score_range'] = 0.0
            results['agreement'] = True
            results['confidence'] = 'medium'

        return results

    def batch_consensus_dock(
        self,
        receptor_pdbqt: str,
        ligand_pdbqts: List[str],
        ligand_ids: List[str],
        center: Tuple[float, float, float],
        size: Tuple[float, float, float],
        output_dir: str,
        max_workers: int = 4,
        methods: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Run consensus docking for multiple ligands in parallel.

        Returns:
            DataFrame with consensus results
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        results = []

        # Parallel execution
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}

            for ligand_pdbqt, ligand_id in zip(ligand_pdbqts, ligand_ids):
                ligand_output = output_dir / ligand_id

                future = executor.submit(
                    self.consensus_dock,
                    receptor_pdbqt,
                    ligand_pdbqt,
                    ligand_id,
                    center,
                    size,
                    str(ligand_output),
                    methods
                )

                futures[future] = ligand_id

            # Collect results
            for future in as_completed(futures):
                ligand_id = futures[future]

                try:
                    result = future.result()
                    results.append(result)
                    logger.info(f"✓ Consensus docking complete: {ligand_id}")
                except Exception as e:
                    logger.error(f"✗ Consensus docking failed for {ligand_id}: {e}")

        # Convert to DataFrame
        df_data = []

        for result in results:
            row = {
                'ligand_id': result['ligand_id'],
                'consensus_score': result['consensus_score'],
                'score_std': result['score_std'],
                'score_range': result['score_range'],
                'confidence': result['confidence'],
                'agreement': result['agreement']
            }

            # Add individual method scores
            for method, data in result['methods'].items():
                row[f'{method}_score'] = data['score']

            df_data.append(row)

        df = pd.DataFrame(df_data)

        # Sort by consensus score
        df = df.sort_values('consensus_score')

        return df

    def analyze_consensus(self, df: pd.DataFrame) -> Dict:
        """
        Analyze consensus docking results.

        Returns:
            Statistics and insights
        """
        n_total = len(df)
        n_high_confidence = (df['confidence'] == 'high').sum()
        n_low_confidence = (df['confidence'] == 'low').sum()

        high_conf_rate = n_high_confidence / n_total if n_total > 0 else 0

        # Find discrepant molecules (methods disagree)
        discrepant = df[df['score_range'] > 3.0]

        analysis = {
            'total_molecules': n_total,
            'high_confidence': n_high_confidence,
            'low_confidence': n_low_confidence,
            'high_confidence_rate': high_conf_rate,
            'mean_consensus_score': df['consensus_score'].mean(),
            'mean_score_std': df['score_std'].mean(),
            'discrepant_molecules': len(discrepant),
            'discrepant_ids': discrepant['ligand_id'].tolist() if not discrepant.empty else []
        }

        return analysis


def install_smina():
    """
    Install Smina (Vina fork) for consensus docking.

    Note: This is a helper function. User should run manually.
    """
    import platform

    system = platform.system()

    if system == 'Darwin':  # macOS
        print("Installing Smina on macOS:")
        print("brew install smina")
        print("\nOr download from: https://sourceforge.net/projects/smina/")

    elif system == 'Linux':
        print("Installing Smina on Linux:")
        print("conda install -c bioconda smina")
        print("\nOr download from: https://sourceforge.net/projects/smina/")

    else:
        print(f"Unsupported system: {system}")
        print("Visit: https://sourceforge.net/projects/smina/")
