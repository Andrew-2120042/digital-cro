"""
Protein Pocket Comparison Module

Compare binding pockets across proteins to:
- Predict drug repurposing opportunities
- Identify off-target binding risks
- Explain selectivity patterns
"""

from pathlib import Path
from typing import List, Dict, Tuple, Optional
import logging
import json
import numpy as np
from collections import Counter

from Bio.PDB import PDBParser

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PocketAnalyzer:
    """
    Analyze and compare protein binding pockets.
    """

    def __init__(self):
        self.parser = PDBParser(QUIET=True)

    def extract_pocket_residues(
        self,
        pdb_file: str,
        pocket_center: Tuple[float, float, float],
        pocket_radius: float = 8.0
    ) -> Dict:
        """
        Extract residues within pocket radius.

        Args:
            pdb_file: Path to protein PDB
            pocket_center: (x, y, z) coordinates of pocket center
            pocket_radius: Radius to define pocket (Angstroms)

        Returns:
            Dict with pocket residues and properties
        """
        structure = self.parser.get_structure('protein', pdb_file)

        pocket_residues = []
        pocket_atoms = []

        cx, cy, cz = pocket_center

        for model in structure:
            for chain in model:
                for residue in chain:
                    # Skip water and hetero atoms
                    if residue.id[0] != ' ':
                        continue

                    # Check if any atom is within pocket radius
                    for atom in residue:
                        coord = atom.get_coord()
                        distance = np.sqrt(
                            (coord[0] - cx)**2 +
                            (coord[1] - cy)**2 +
                            (coord[2] - cz)**2
                        )

                        if distance <= pocket_radius:
                            pocket_residues.append(residue)
                            pocket_atoms.extend(list(residue.get_atoms()))
                            break

        # Get residue composition
        residue_names = [res.get_resname() for res in pocket_residues]
        residue_counts = Counter(residue_names)

        # Calculate pocket properties
        properties = self._calculate_pocket_properties(pocket_residues, pocket_atoms)

        return {
            'residues': pocket_residues,
            'residue_names': residue_names,
            'residue_counts': dict(residue_counts),
            'num_residues': len(pocket_residues),
            'properties': properties
        }

    def _calculate_pocket_properties(
        self,
        residues: List,
        atoms: List
    ) -> Dict:
        """Calculate chemical properties of pocket"""

        # Classify residues by property
        hydrophobic = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO']
        polar = ['SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN']
        positive = ['LYS', 'ARG', 'HIS']
        negative = ['ASP', 'GLU']

        residue_names = [res.get_resname() for res in residues]

        counts = {
            'hydrophobic': sum(1 for r in residue_names if r in hydrophobic),
            'polar': sum(1 for r in residue_names if r in polar),
            'positive': sum(1 for r in residue_names if r in positive),
            'negative': sum(1 for r in residue_names if r in negative)
        }

        total = len(residue_names)

        fractions = {
            f'{k}_fraction': v / total if total > 0 else 0
            for k, v in counts.items()
        }

        # Calculate pocket volume (approximate using atoms)
        if len(atoms) > 0:
            coords = np.array([atom.get_coord() for atom in atoms])

            # Bounding box volume
            min_coords = coords.min(axis=0)
            max_coords = coords.max(axis=0)
            dimensions = max_coords - min_coords
            volume = np.prod(dimensions)
        else:
            volume = 0

        return {
            **counts,
            **fractions,
            'volume': volume,
            'total_residues': total
        }

    def compare_pockets(
        self,
        pocket1: Dict,
        pocket2: Dict
    ) -> Dict:
        """
        Compare two binding pockets.

        Returns:
            Similarity scores and analysis
        """
        # 1. Residue composition similarity (Tanimoto-like)
        residues1 = set(pocket1['residue_names'])
        residues2 = set(pocket2['residue_names'])

        intersection = len(residues1 & residues2)
        union = len(residues1 | residues2)

        residue_similarity = intersection / union if union > 0 else 0

        # 2. Property similarity (chemical character)
        props1 = pocket1['properties']
        props2 = pocket2['properties']

        property_similarity = self._calculate_property_similarity(props1, props2)

        # 3. Size similarity
        size1 = props1['total_residues']
        size2 = props2['total_residues']

        size_similarity = min(size1, size2) / max(size1, size2) if max(size1, size2) > 0 else 0

        # 4. Overall similarity (weighted average)
        overall_similarity = (
            0.4 * residue_similarity +
            0.4 * property_similarity +
            0.2 * size_similarity
        )

        # Interpretation
        interpretation = self._interpret_similarity(overall_similarity)

        return {
            'overall_similarity': overall_similarity,
            'residue_similarity': residue_similarity,
            'property_similarity': property_similarity,
            'size_similarity': size_similarity,
            'interpretation': interpretation,
            'pocket1_size': size1,
            'pocket2_size': size2,
            'common_residues': list(residues1 & residues2),
            'unique_to_pocket1': list(residues1 - residues2),
            'unique_to_pocket2': list(residues2 - residues1)
        }

    def _calculate_property_similarity(
        self,
        props1: Dict,
        props2: Dict
    ) -> float:
        """Calculate similarity of chemical properties"""

        # Compare property fractions
        properties = ['hydrophobic_fraction', 'polar_fraction', 'positive_fraction', 'negative_fraction']

        differences = []

        for prop in properties:
            val1 = props1.get(prop, 0)
            val2 = props2.get(prop, 0)

            # Normalize difference (0 = identical, 1 = completely different)
            diff = abs(val1 - val2)
            differences.append(diff)

        # Average difference
        avg_diff = np.mean(differences)

        # Convert to similarity (1 = identical, 0 = completely different)
        similarity = 1 - avg_diff

        return similarity

    def _interpret_similarity(self, similarity: float) -> str:
        """Interpret similarity score"""

        if similarity >= 0.8:
            return "Very High - Strong repurposing potential, expect cross-reactivity"
        elif similarity >= 0.6:
            return "High - Moderate repurposing potential, some cross-reactivity expected"
        elif similarity >= 0.4:
            return "Moderate - Limited repurposing potential, selectivity likely"
        elif similarity >= 0.2:
            return "Low - Different pockets, high selectivity expected"
        else:
            return "Very Low - Completely different pockets, no cross-reactivity expected"

    def analyze_selectivity(
        self,
        pocket_comparisons: List[Dict],
        target_names: List[str]
    ) -> Dict:
        """
        Analyze selectivity pattern across multiple targets.

        Args:
            pocket_comparisons: List of pairwise comparisons
            target_names: Names of targets being compared

        Returns:
            Selectivity analysis and recommendations
        """
        # Find most similar pocket pairs
        sorted_comparisons = sorted(
            pocket_comparisons,
            key=lambda x: x['overall_similarity'],
            reverse=True
        )

        # Most similar pair (highest cross-reactivity risk)
        most_similar = sorted_comparisons[0] if sorted_comparisons else None

        # Most different pair (best selectivity)
        least_similar = sorted_comparisons[-1] if sorted_comparisons else None

        # Average similarity
        avg_similarity = np.mean([c['overall_similarity'] for c in pocket_comparisons]) if pocket_comparisons else 0

        # Selectivity assessment
        if avg_similarity >= 0.7:
            selectivity_level = "Low - High cross-reactivity risk"
        elif avg_similarity >= 0.5:
            selectivity_level = "Moderate - Some cross-reactivity expected"
        else:
            selectivity_level = "High - Good target selectivity"

        return {
            'avg_similarity': avg_similarity,
            'most_similar_pair': most_similar,
            'least_similar_pair': least_similar,
            'selectivity_level': selectivity_level,
            'num_comparisons': len(pocket_comparisons)
        }

    def generate_repurposing_recommendations(
        self,
        pocket_comparisons: List[Dict],
        target_names: List[str],
        binding_data: Optional[Dict] = None
    ) -> List[str]:
        """
        Generate drug repurposing recommendations.

        Args:
            pocket_comparisons: Pairwise pocket comparisons
            target_names: Target protein names
            binding_data: Optional binding affinity data per target

        Returns:
            List of recommendations
        """
        recommendations = []

        # Find high similarity pairs (> 0.7)
        similar_pairs = [
            c for c in pocket_comparisons
            if c['overall_similarity'] > 0.7
        ]

        if similar_pairs:
            recommendations.append(
                f"**Drug Repurposing Opportunity**: Found {len(similar_pairs)} pocket pair(s) "
                f"with >70% similarity. Molecules binding one target likely bind the other."
            )

            for pair in similar_pairs[:3]:  # Top 3
                recommendations.append(
                    f"  • Similarity: {pair['overall_similarity']:.2f} - {pair['interpretation']}"
                )

        # Find low similarity pairs (< 0.3)
        different_pairs = [
            c for c in pocket_comparisons
            if c['overall_similarity'] < 0.3
        ]

        if different_pairs:
            recommendations.append(
                f"**Selectivity Opportunity**: Found {len(different_pairs)} highly different pocket pair(s). "
                f"Design selective molecules targeting unique residues."
            )

        # Binding data insights
        if binding_data:
            recommendations.append(
                "**Cross-Activity Analysis**: Combine pocket similarity with binding data "
                "to identify true repurposing candidates."
            )

        return recommendations


def compare_multi_target_pockets(
    multi_target_dir: str,
    output_dir: Optional[str] = None
) -> Dict:
    """
    Compare pockets from multi-target screening results.

    Args:
        multi_target_dir: Directory with multi-target results
        output_dir: Where to save analysis

    Returns:
        Complete pocket comparison analysis
    """
    multi_target_dir = Path(multi_target_dir)

    # Load manifest to get pocket info
    manifest_path = multi_target_dir / "multi_target_summary.json"

    if not manifest_path.exists():
        raise FileNotFoundError("Multi-target summary not found")

    with open(manifest_path, 'r') as f:
        summary = json.load(f)

    analyzer = PocketAnalyzer()

    # Extract pockets for each target
    pockets = {}
    target_names = []

    for target_id in summary['targets'].keys():
        target_dir = multi_target_dir / f"target_{target_id}"

        # Find protein PDB
        protein_files = list(target_dir.glob("proteins/*_cleaned.pdb"))

        if not protein_files:
            logger.warning(f"No protein found for {target_id}")
            continue

        protein_pdb = protein_files[0]

        # Get pocket info from manifest
        target_manifest = target_dir / "manifest.json"

        if target_manifest.exists():
            with open(target_manifest, 'r') as f:
                target_info = json.load(f)
                pocket_center = tuple(target_info['pocket_center'])
                pocket_size = target_info['pocket_size']
        else:
            # Default
            pocket_center = (0, 0, 0)
            pocket_size = [20, 20, 20]

        # Extract pocket
        pocket = analyzer.extract_pocket_residues(
            str(protein_pdb),
            pocket_center,
            pocket_radius=8.0
        )

        pockets[target_id] = pocket
        target_names.append(target_id)

    logger.info(f"✓ Extracted {len(pockets)} pockets")

    # Pairwise comparisons
    comparisons = []
    comparison_matrix = {}

    for i, target1 in enumerate(target_names):
        for j, target2 in enumerate(target_names):
            if i >= j:
                continue

            comparison = analyzer.compare_pockets(
                pockets[target1],
                pockets[target2]
            )

            comparison['target1'] = target1
            comparison['target2'] = target2

            comparisons.append(comparison)
            comparison_matrix[f"{target1}_vs_{target2}"] = comparison['overall_similarity']

    # Selectivity analysis
    selectivity = analyzer.analyze_selectivity(comparisons, target_names)

    # Recommendations
    recommendations = analyzer.generate_repurposing_recommendations(
        comparisons,
        target_names
    )

    results = {
        'pockets': {k: {
            'num_residues': v['num_residues'],
            'properties': v['properties'],
            'residue_counts': v['residue_counts']
        } for k, v in pockets.items()},
        'comparisons': comparisons,
        'comparison_matrix': comparison_matrix,
        'selectivity': selectivity,
        'recommendations': recommendations,
        'target_names': target_names
    }

    # Save results
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save JSON
        with open(output_dir / 'pocket_comparison.json', 'w') as f:
            # Remove non-serializable objects
            save_results = {
                'comparisons': [{
                    'target1': c['target1'],
                    'target2': c['target2'],
                    'overall_similarity': c['overall_similarity'],
                    'residue_similarity': c['residue_similarity'],
                    'property_similarity': c['property_similarity'],
                    'interpretation': c['interpretation']
                } for c in comparisons],
                'selectivity': selectivity,
                'recommendations': recommendations
            }
            json.dump(save_results, f, indent=2)

        # Save recommendations as text
        with open(output_dir / 'repurposing_recommendations.txt', 'w') as f:
            f.write("Protein Pocket Comparison Analysis\n")
            f.write("=" * 50 + "\n\n")

            f.write(f"Targets Analyzed: {', '.join(target_names)}\n")
            f.write(f"Average Pocket Similarity: {selectivity['avg_similarity']:.2f}\n")
            f.write(f"Selectivity Level: {selectivity['selectivity_level']}\n\n")

            f.write("Recommendations:\n")
            f.write("-" * 50 + "\n")

            for i, rec in enumerate(recommendations, 1):
                f.write(f"{i}. {rec}\n\n")

        logger.info(f"✓ Pocket comparison saved to {output_dir}")

    return results
