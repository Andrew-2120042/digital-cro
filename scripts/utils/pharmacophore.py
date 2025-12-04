"""
Pharmacophore Analysis Module

Generate pharmacophore hypotheses from top-binding molecules.
Identifies common 3D features (H-bond donors, acceptors, aromatic rings, etc).
"""

from pathlib import Path
from typing import List, Dict, Tuple, Optional
import logging
import pandas as pd
from collections import Counter

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PharmacophoreAnalyzer:
    """
    Analyze molecules to identify pharmacophore features.
    """

    def __init__(self):
        # Load feature definitions
        fdef_path = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        self.factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)

    def extract_features(self, mol: Chem.Mol) -> List[Dict]:
        """
        Extract pharmacophore features from molecule.

        Returns list of features with:
        - type: Feature type (Donor, Acceptor, Aromatic, etc)
        - position: 3D coordinates
        - atoms: Atom indices
        """
        if mol is None:
            return []

        # Generate 3D coordinates if not present
        if mol.GetNumConformers() == 0:
            try:
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol)
            except:
                # If 3D generation fails, continue without coordinates
                pass

        features = []
        raw_features = self.factory.GetFeaturesForMol(mol)

        for feat in raw_features:
            feature_dict = {
                'type': feat.GetFamily(),
                'atoms': feat.GetAtomIds(),
            }

            # Get 3D position
            if mol.GetNumConformers() > 0:
                try:
                    conf = mol.GetConformer()
                    pos = feat.GetPos()
                    feature_dict['position'] = (pos.x, pos.y, pos.z)
                except:
                    pass

            features.append(feature_dict)

        return features

    def analyze_molecules(
        self,
        smiles_list: List[str],
        mol_ids: List[str]
    ) -> pd.DataFrame:
        """
        Extract features from multiple molecules.

        Returns DataFrame with:
        - mol_id
        - smiles
        - num_donors
        - num_acceptors
        - num_aromatic
        - num_hydrophobic
        - features (detailed list)
        """
        results = []

        for smiles, mol_id in zip(smiles_list, mol_ids):
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                continue

            features = self.extract_features(mol)

            # Count feature types
            feature_counts = Counter(f['type'] for f in features)

            results.append({
                'mol_id': mol_id,
                'smiles': smiles,
                'num_donors': feature_counts.get('Donor', 0),
                'num_acceptors': feature_counts.get('Acceptor', 0),
                'num_aromatic': feature_counts.get('Aromatic', 0),
                'num_hydrophobic': feature_counts.get('Hydrophobe', 0),
                'num_positive': feature_counts.get('PosIonizable', 0),
                'num_negative': feature_counts.get('NegIonizable', 0),
                'total_features': len(features),
                'features': features
            })

        return pd.DataFrame(results)

    def find_common_features(
        self,
        feature_df: pd.DataFrame,
        min_occurrence: float = 0.7
    ) -> Dict[str, int]:
        """
        Find features that occur in most molecules.

        Args:
            feature_df: DataFrame from analyze_molecules
            min_occurrence: Minimum fraction of molecules (0.7 = 70%)

        Returns:
            Dict of feature type -> minimum count
        """
        n_molecules = len(feature_df)
        threshold = int(n_molecules * min_occurrence)

        common_features = {}

        # For each feature type
        for feat_type in ['donors', 'acceptors', 'aromatic', 'hydrophobic', 'positive', 'negative']:
            col_name = f'num_{feat_type}'

            if col_name not in feature_df.columns:
                continue

            # Count molecules with at least N of this feature
            for min_count in range(1, 6):  # Check 1-5
                n_with_feature = (feature_df[col_name] >= min_count).sum()

                if n_with_feature >= threshold:
                    if feat_type not in common_features:
                        common_features[feat_type] = min_count

        return common_features

    def generate_hypothesis(
        self,
        feature_df: pd.DataFrame,
        min_occurrence: float = 0.7
    ) -> Dict:
        """
        Generate pharmacophore hypothesis from features.

        Returns:
            - common_features: Dict of essential features
            - hypothesis_text: Human-readable description
            - n_molecules: Number of molecules analyzed
        """
        common = self.find_common_features(feature_df, min_occurrence)

        # Generate hypothesis text
        if not common:
            hypothesis = "No common pharmacophore features found across the majority of molecules."
        else:
            parts = []

            feature_names = {
                'donors': 'hydrogen bond donor(s)',
                'acceptors': 'hydrogen bond acceptor(s)',
                'aromatic': 'aromatic ring(s)',
                'hydrophobic': 'hydrophobic region(s)',
                'positive': 'positive ionizable group(s)',
                'negative': 'negative ionizable group(s)'
            }

            for feat_type, count in common.items():
                name = feature_names.get(feat_type, feat_type)
                parts.append(f"{count} {name}")

            hypothesis = f"Pharmacophore hypothesis: Active molecules contain at least " + ", ".join(parts) + "."

        return {
            'common_features': common,
            'hypothesis_text': hypothesis,
            'n_molecules': len(feature_df),
            'min_occurrence': min_occurrence
        }

    def calculate_feature_importance(
        self,
        feature_df: pd.DataFrame,
        binding_affinities: List[float]
    ) -> pd.DataFrame:
        """
        Correlate features with binding affinity.

        Returns DataFrame with:
        - feature_type
        - correlation (with binding affinity)
        - mean_count (average count in molecules)
        """
        feature_df = feature_df.copy()
        feature_df['affinity'] = binding_affinities

        correlations = []

        for feat_col in ['num_donors', 'num_acceptors', 'num_aromatic', 'num_hydrophobic']:
            if feat_col not in feature_df.columns:
                continue

            corr = feature_df[feat_col].corr(feature_df['affinity'])
            mean_count = feature_df[feat_col].mean()

            correlations.append({
                'feature_type': feat_col.replace('num_', ''),
                'correlation': corr if not pd.isna(corr) else 0.0,
                'mean_count': mean_count
            })

        return pd.DataFrame(correlations)


def analyze_pharmacophore_from_results(
    results_csv: str,
    top_n: int = 20,
    min_occurrence: float = 0.7,
    output_dir: Optional[str] = None
) -> Dict:
    """
    Analyze pharmacophore from screening results.

    Args:
        results_csv: Path to final_results.csv
        top_n: Number of top hits to analyze
        min_occurrence: Minimum occurrence fraction for common features
        output_dir: Where to save results

    Returns:
        Dict with hypothesis and feature analysis
    """
    df = pd.read_csv(results_csv)

    # Get top N hits
    df_top = df.nsmallest(top_n, 'binding_affinity')

    # Extract features
    analyzer = PharmacophoreAnalyzer()

    logger.info(f"Analyzing pharmacophore for top {len(df_top)} molecules...")

    feature_df = analyzer.analyze_molecules(
        smiles_list=df_top['smiles'].tolist(),
        mol_ids=df_top['ligand_id'].tolist()
    )

    # Generate hypothesis
    hypothesis = analyzer.generate_hypothesis(feature_df, min_occurrence)

    # Calculate feature importance
    importance = analyzer.calculate_feature_importance(
        feature_df,
        df_top['binding_affinity'].tolist()
    )

    results = {
        'hypothesis': hypothesis,
        'feature_df': feature_df,
        'importance': importance,
        'top_molecules': df_top
    }

    # Save results
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save feature data
        feature_df.drop('features', axis=1).to_csv(
            output_dir / 'pharmacophore_features.csv',
            index=False
        )

        # Save importance
        importance.to_csv(
            output_dir / 'feature_importance.csv',
            index=False
        )

        # Save hypothesis as text
        with open(output_dir / 'pharmacophore_hypothesis.txt', 'w') as f:
            f.write(hypothesis['hypothesis_text'])
            f.write(f"\n\nAnalyzed {hypothesis['n_molecules']} molecules")
            f.write(f"\nMinimum occurrence: {hypothesis['min_occurrence']*100}%")
            f.write(f"\n\nCommon features:\n")
            for feat, count in hypothesis['common_features'].items():
                f.write(f"  - {feat}: {count}\n")

        logger.info(f"✓ Pharmacophore analysis saved to {output_dir}")

    return results
