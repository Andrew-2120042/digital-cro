"""
SAR Analysis Module

Analyze structure-activity relationships:
- Identify molecular scaffolds
- Detect R-group variations
- Analyze substituent effects
- Generate matched molecular pairs
"""

from pathlib import Path
from typing import List, Dict, Tuple, Optional
import logging
import pandas as pd
import numpy as np
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SARAnalyzer:
    """
    Analyze structure-activity relationships in molecule sets.
    """

    def __init__(self):
        self.scaffolds = {}
        self.molecule_data = {}

    def identify_scaffolds(
        self,
        molecules: List[Chem.Mol],
        mol_ids: List[str],
        min_molecules_per_scaffold: int = 2
    ) -> Dict[str, List[str]]:
        """
        Identify Murcko scaffolds in molecule set.

        Returns:
            Dict mapping scaffold SMILES -> list of molecule IDs
        """
        scaffold_to_mols = defaultdict(list)

        for mol, mol_id in zip(molecules, mol_ids):
            if mol is None:
                continue

            try:
                # Get Murcko scaffold
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smiles = Chem.MolToSmiles(scaffold)

                scaffold_to_mols[scaffold_smiles].append(mol_id)

            except Exception as e:
                logger.warning(f"Could not get scaffold for {mol_id}: {e}")
                continue

        # Filter scaffolds with minimum molecules
        filtered = {
            scaffold: mols
            for scaffold, mols in scaffold_to_mols.items()
            if len(mols) >= min_molecules_per_scaffold
        }

        logger.info(f"✓ Found {len(filtered)} scaffolds with ≥{min_molecules_per_scaffold} molecules")

        return filtered

    def analyze_scaffold_series(
        self,
        df: pd.DataFrame,
        scaffold_smiles: str,
        mol_ids: List[str]
    ) -> pd.DataFrame:
        """
        Analyze a series of molecules with the same scaffold.

        Returns DataFrame with:
        - mol_id
        - smiles
        - binding_affinity
        - mw (molecular weight)
        - logp
        - num_rotatable_bonds
        - delta_affinity (vs best in series)
        """
        # Get data for this scaffold
        scaffold_df = df[df['ligand_id'].isin(mol_ids)].copy()

        # Calculate molecular properties
        properties = []

        for _, row in scaffold_df.iterrows():
            mol = Chem.MolFromSmiles(row['smiles'])

            if mol:
                properties.append({
                    'mol_id': row['ligand_id'],
                    'smiles': row['smiles'],
                    'binding_affinity': row['binding_affinity'],
                    'mw': Descriptors.MolWt(mol),
                    'logp': Descriptors.MolLogP(mol),
                    'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                    'num_hbd': Descriptors.NumHDonors(mol),
                    'num_hba': Descriptors.NumHAcceptors(mol),
                    'tpsa': Descriptors.TPSA(mol)
                })

        prop_df = pd.DataFrame(properties)

        # Calculate delta from best binder
        if len(prop_df) > 0:
            best_affinity = prop_df['binding_affinity'].min()
            prop_df['delta_affinity'] = prop_df['binding_affinity'] - best_affinity

            # Sort by affinity
            prop_df = prop_df.sort_values('binding_affinity')

        return prop_df

    def find_matched_molecular_pairs(
        self,
        molecules: List[Chem.Mol],
        mol_ids: List[str],
        affinities: List[float],
        max_structural_difference: int = 3
    ) -> List[Dict]:
        """
        Find matched molecular pairs (MMPs).

        MMPs are molecules that differ by a single structural change.

        Returns list of MMPs with:
        - mol1_id, mol2_id
        - mol1_affinity, mol2_affinity
        - delta_affinity
        - structural_change description
        """
        mmps = []

        n = len(molecules)

        for i in range(n):
            for j in range(i + 1, n):
                mol1, mol2 = molecules[i], molecules[j]

                if mol1 is None or mol2 is None:
                    continue

                # Calculate Tanimoto similarity
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)

                similarity = Chem.DataStructs.TanimotoSimilarity(fp1, fp2)

                # Consider as MMP if highly similar (>0.7)
                if similarity > 0.7:
                    # Calculate difference in heavy atoms
                    diff = abs(mol1.GetNumHeavyAtoms() - mol2.GetNumHeavyAtoms())

                    if diff <= max_structural_difference:
                        delta_affinity = affinities[j] - affinities[i]

                        mmps.append({
                            'mol1_id': mol_ids[i],
                            'mol2_id': mol_ids[j],
                            'mol1_affinity': affinities[i],
                            'mol2_affinity': affinities[j],
                            'delta_affinity': delta_affinity,
                            'similarity': similarity,
                            'heavy_atom_diff': diff,
                            'improvement': 'Better' if delta_affinity < 0 else 'Worse'
                        })

        # Sort by absolute delta affinity (biggest changes first)
        mmps.sort(key=lambda x: abs(x['delta_affinity']), reverse=True)

        logger.info(f"✓ Found {len(mmps)} matched molecular pairs")

        return mmps

    def analyze_property_correlations(
        self,
        df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Analyze correlations between molecular properties and binding.

        Returns DataFrame with:
        - property
        - correlation
        - p_value (if significant)
        """
        from scipy import stats

        properties = ['mw', 'logp', 'num_rotatable_bonds', 'num_hbd', 'num_hba', 'tpsa']

        correlations = []

        for prop in properties:
            if prop not in df.columns:
                continue

            # Calculate Pearson correlation
            corr, p_value = stats.pearsonr(df[prop], df['binding_affinity'])

            correlations.append({
                'property': prop.replace('_', ' ').title(),
                'correlation': corr,
                'p_value': p_value,
                'significant': 'Yes' if p_value < 0.05 else 'No'
            })

        corr_df = pd.DataFrame(correlations)
        corr_df = corr_df.sort_values('correlation', key=abs, ascending=False)

        return corr_df

    def generate_sar_recommendations(
        self,
        scaffold_analyses: List[pd.DataFrame],
        mmps: List[Dict],
        correlations: pd.DataFrame
    ) -> List[str]:
        """
        Generate SAR-based recommendations for optimization.
        """
        recommendations = []

        # 1. Best scaffolds
        if scaffold_analyses:
            best_scaffold = min(scaffold_analyses, key=lambda df: df['binding_affinity'].min() if len(df) > 0 else float('inf'))
            if len(best_scaffold) > 0:
                best_affinity = best_scaffold['binding_affinity'].min()

                recommendations.append(
                    f"**Best Scaffold**: Focus on scaffold with best binder ({best_affinity:.2f} kcal/mol). "
                    f"This core structure shows the most promise."
                )

        # 2. Property-based recommendations
        if not correlations.empty:
            # Find strongest correlations
            strong_corr = correlations[correlations['significant'] == 'Yes'].head(3)

            for _, row in strong_corr.iterrows():
                prop = row['property']
                corr = row['correlation']

                if corr < -0.3:
                    recommendations.append(
                        f"**Increase {prop}**: Strong negative correlation ({corr:.2f}) suggests "
                        f"increasing {prop} improves binding."
                    )
                elif corr > 0.3:
                    recommendations.append(
                        f"**Decrease {prop}**: Positive correlation ({corr:.2f}) suggests "
                        f"decreasing {prop} improves binding."
                    )

        # 3. MMP-based recommendations
        if mmps:
            # Find best improvements
            improvements = [mmp for mmp in mmps if mmp['delta_affinity'] < -1.0]

            if improvements:
                best_mmp = improvements[0]
                recommendations.append(
                    f"**Structural Optimization**: Transformation from {best_mmp['mol1_id']} to "
                    f"{best_mmp['mol2_id']} improved binding by {abs(best_mmp['delta_affinity']):.2f} kcal/mol. "
                    f"Apply similar modifications to other leads."
                )

        # 4. General recommendations
        if recommendations:
            recommendations.append(
                "**Next Steps**: Synthesize 5-10 analogs based on above recommendations and re-screen. "
                "Focus on maintaining the core scaffold while optimizing identified properties."
            )
        else:
            recommendations.append(
                "**Continue Exploration**: No strong SAR trends identified yet. "
                "Consider expanding chemical diversity or testing more analogs of promising scaffolds."
            )

        return recommendations


def analyze_sar_from_results(
    results_csv: str,
    min_molecules_per_scaffold: int = 2,
    output_dir: Optional[str] = None
) -> Dict:
    """
    Perform complete SAR analysis from screening results.

    Args:
        results_csv: Path to final_results.csv
        min_molecules_per_scaffold: Minimum molecules to analyze scaffold
        output_dir: Where to save results

    Returns:
        Dict with SAR analysis results
    """
    df = pd.read_csv(results_csv)

    logger.info(f"Analyzing SAR for {len(df)} molecules...")

    # Convert SMILES to molecules
    molecules = [Chem.MolFromSmiles(smi) for smi in df['smiles']]
    mol_ids = df['ligand_id'].tolist()
    affinities = df['binding_affinity'].tolist()

    analyzer = SARAnalyzer()

    # 1. Identify scaffolds
    scaffolds = analyzer.identify_scaffolds(
        molecules,
        mol_ids,
        min_molecules_per_scaffold
    )

    # 2. Analyze each scaffold series
    scaffold_analyses = []

    for scaffold_smiles, scaffold_mols in list(scaffolds.items())[:10]:  # Top 10 scaffolds
        analysis = analyzer.analyze_scaffold_series(df, scaffold_smiles, scaffold_mols)

        if len(analysis) > 0:
            scaffold_analyses.append({
                'scaffold_smiles': scaffold_smiles,
                'n_molecules': len(scaffold_mols),
                'best_affinity': analysis['binding_affinity'].min(),
                'mean_affinity': analysis['binding_affinity'].mean(),
                'affinity_range': analysis['binding_affinity'].max() - analysis['binding_affinity'].min(),
                'analysis': analysis
            })

    # Sort by best affinity
    scaffold_analyses.sort(key=lambda x: x['best_affinity'])

    # 3. Find matched molecular pairs
    mmps = analyzer.find_matched_molecular_pairs(
        molecules,
        mol_ids,
        affinities,
        max_structural_difference=3
    )

    # 4. Property correlations (using all molecules)
    # Add properties to dataframe
    prop_list = []

    for mol, mol_id in zip(molecules, mol_ids):
        if mol:
            prop_list.append({
                'mol_id': mol_id,
                'mw': Descriptors.MolWt(mol),
                'logp': Descriptors.MolLogP(mol),
                'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                'num_hbd': Descriptors.NumHDonors(mol),
                'num_hba': Descriptors.NumHAcceptors(mol),
                'tpsa': Descriptors.TPSA(mol)
            })

    prop_df = pd.DataFrame(prop_list)
    df_with_props = df.merge(prop_df, left_on='ligand_id', right_on='mol_id', how='left')

    correlations = analyzer.analyze_property_correlations(df_with_props)

    # 5. Generate recommendations
    recommendations = analyzer.generate_sar_recommendations(
        [s['analysis'] for s in scaffold_analyses],
        mmps,
        correlations
    )

    results = {
        'scaffolds': scaffold_analyses,
        'mmps': mmps[:20],  # Top 20 pairs
        'correlations': correlations,
        'recommendations': recommendations,
        'n_molecules': len(df),
        'n_scaffolds': len(scaffolds)
    }

    # Save results
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save scaffold analysis
        if scaffold_analyses:
            scaffold_summary = pd.DataFrame([
                {
                    'scaffold_smiles': s['scaffold_smiles'],
                    'n_molecules': s['n_molecules'],
                    'best_affinity': s['best_affinity'],
                    'mean_affinity': s['mean_affinity'],
                    'affinity_range': s['affinity_range']
                }
                for s in scaffold_analyses
            ])

            scaffold_summary.to_csv(output_dir / 'scaffold_summary.csv', index=False)

        # Save MMPs
        if mmps:
            pd.DataFrame(mmps).to_csv(output_dir / 'matched_pairs.csv', index=False)

        # Save correlations
        correlations.to_csv(output_dir / 'property_correlations.csv', index=False)

        # Save recommendations
        with open(output_dir / 'sar_recommendations.txt', 'w') as f:
            f.write("SAR Analysis Recommendations\n")
            f.write("=" * 50 + "\n\n")

            for i, rec in enumerate(recommendations, 1):
                f.write(f"{i}. {rec}\n\n")

        logger.info(f"✓ SAR analysis saved to {output_dir}")

    return results
