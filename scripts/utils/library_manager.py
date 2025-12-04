"""
Molecule Library Manager

Save, organize, and reuse molecule libraries.
"""

import json
import shutil
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class LibraryManager:
    """
    Manage saved molecule libraries.

    Libraries are stored as:
    - Original file (CSV/SMI)
    - Metadata JSON (name, size, stats, created date)
    """

    def __init__(self, libraries_dir: str = "data/molecule_libraries"):
        self.libraries_dir = Path(libraries_dir)
        self.libraries_dir.mkdir(parents=True, exist_ok=True)

        self.metadata_file = self.libraries_dir / "libraries.json"

        # Initialize metadata file
        if not self.metadata_file.exists():
            self._save_metadata({})

    def _load_metadata(self) -> Dict:
        """Load all library metadata"""
        try:
            with open(self.metadata_file, 'r') as f:
                return json.load(f)
        except:
            return {}

    def _save_metadata(self, metadata: Dict):
        """Save library metadata"""
        with open(self.metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)

    def save_library(
        self,
        name: str,
        file_path: str,
        description: str = "",
        tags: List[str] = None,
        compute_stats: bool = True
    ) -> str:
        """
        Save a molecule library.

        Args:
            name: Library name (must be unique)
            file_path: Path to CSV/SMI file
            description: Optional description
            tags: Optional tags for categorization
            compute_stats: If True, compute molecular properties

        Returns:
            library_id: Unique identifier for this library
        """
        metadata = self._load_metadata()

        # Check if name exists
        for lib_id, lib in metadata.items():
            if lib['name'] == name:
                raise ValueError(f"Library '{name}' already exists")

        # Generate library ID
        library_id = name.lower().replace(' ', '_').replace('-', '_')
        counter = 1
        original_id = library_id

        while library_id in metadata:
            library_id = f"{original_id}_{counter}"
            counter += 1

        # Create library directory
        lib_dir = self.libraries_dir / library_id
        lib_dir.mkdir(exist_ok=True)

        # Copy file
        file_path = Path(file_path)
        saved_file = lib_dir / file_path.name
        shutil.copy(file_path, saved_file)

        # Read molecules
        df = self._read_library_file(saved_file)

        # Compute statistics
        stats = {}
        if compute_stats and df is not None:
            stats = self._compute_library_stats(df)

        # Create metadata
        lib_metadata = {
            'library_id': library_id,
            'name': name,
            'description': description,
            'tags': tags or [],
            'file_path': str(saved_file),
            'file_name': file_path.name,
            'file_type': file_path.suffix.lower(),
            'num_molecules': len(df) if df is not None else 0,
            'created_at': datetime.now().isoformat(),
            'last_used': None,
            'use_count': 0,
            'stats': stats
        }

        # Save metadata
        metadata[library_id] = lib_metadata
        self._save_metadata(metadata)

        logger.info(f"✓ Library saved: {name} ({len(df)} molecules)")

        return library_id

    def _read_library_file(self, file_path: Path) -> Optional[pd.DataFrame]:
        """Read molecule library file with proper column normalization"""
        try:
            if file_path.suffix.lower() == '.csv':
                # Read CSV with proper parsing
                df = pd.read_csv(file_path, dtype=str)

                # Normalize column names to lowercase
                df.columns = df.columns.str.lower().str.strip()

                # Handle different column name variations
                column_mapping = {}

                # Find SMILES column
                smiles_candidates = ['smiles', 'smile', 'smi', 'structure']
                for col in df.columns:
                    if col in smiles_candidates:
                        column_mapping[col] = 'smiles'
                        break

                # Find ID column
                id_candidates = ['id', 'name', 'mol_id', 'molecule_id', 'compound_id']
                for col in df.columns:
                    if col in id_candidates and col not in column_mapping:
                        column_mapping[col] = 'id'
                        break

                # Rename columns
                if column_mapping:
                    df = df.rename(columns=column_mapping)

                # Ensure required columns exist
                if 'smiles' not in df.columns:
                    # If only 2 columns and no header match, assume first is ID, second is SMILES
                    if len(df.columns) == 2:
                        df.columns = ['id', 'smiles']
                    else:
                        logger.error(f"Could not find SMILES column in {file_path}")
                        return None

                if 'id' not in df.columns:
                    # Generate IDs if missing
                    df['id'] = [f'mol_{i:06d}' for i in range(len(df))]

                # Return only id and smiles columns
                return df[['id', 'smiles']].copy()

            elif file_path.suffix.lower() in ['.smi', '.txt']:
                # Read SMI file (tab-separated: SMILES<tab>ID)
                lines = file_path.read_text().strip().split('\n')

                smiles_list = []
                id_list = []

                for idx, line in enumerate(lines):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    if '\t' in line:
                        parts = line.split('\t')
                        smiles_list.append(parts[0].strip())
                        id_list.append(parts[1].strip() if len(parts) > 1 else f'mol_{idx:06d}')
                    else:
                        # No tab, just SMILES
                        smiles_list.append(line)
                        id_list.append(f'mol_{idx:06d}')

                return pd.DataFrame({'id': id_list, 'smiles': smiles_list})
        except Exception as e:
            logger.error(f"Error reading library file {file_path}: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None

    def _compute_library_stats(self, df: pd.DataFrame) -> Dict:
        """Compute statistics for library"""
        stats = {}

        if 'smiles' not in df.columns:
            return stats

        # Basic stats
        stats['num_molecules'] = len(df)

        # Molecular properties
        mw_list = []
        logp_list = []

        for smiles in df['smiles'].head(100):  # Sample first 100 for speed
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mw_list.append(Descriptors.MolWt(mol))
                    logp_list.append(Descriptors.MolLogP(mol))
            except:
                continue

        if mw_list:
            stats['mean_mw'] = sum(mw_list) / len(mw_list)
            stats['min_mw'] = min(mw_list)
            stats['max_mw'] = max(mw_list)

        if logp_list:
            stats['mean_logp'] = sum(logp_list) / len(logp_list)

        return stats

    def list_libraries(self) -> List[Dict]:
        """
        List all saved libraries.

        Returns sorted by name.
        """
        metadata = self._load_metadata()
        libraries = list(metadata.values())
        libraries.sort(key=lambda x: x['name'])
        return libraries

    def get_library(self, library_id: str) -> Optional[Dict]:
        """Get library metadata by ID"""
        metadata = self._load_metadata()
        return metadata.get(library_id)

    def get_library_path(self, library_id: str) -> Optional[str]:
        """
        Get file path for library.

        IMPORTANT: This returns the path to a .smi file (tab-separated format)
        that the docking workflow expects: SMILES<tab>ID

        If the original library was CSV, this method converts it to .smi format.
        """
        lib = self.get_library(library_id)
        if not lib:
            return None

        original_path = Path(lib['file_path'])

        # If it's already a .smi file, return it
        if original_path.suffix.lower() in ['.smi', '.txt']:
            return str(original_path)

        # If it's CSV, convert to .smi format
        if original_path.suffix.lower() == '.csv':
            # Create/update .smi version
            smi_path = original_path.parent / f"{original_path.stem}.smi"

            # Read the CSV properly
            df = self._read_library_file(original_path)

            if df is not None and 'id' in df.columns and 'smiles' in df.columns:
                # Write as tab-separated .smi file
                with open(smi_path, 'w') as f:
                    for _, row in df.iterrows():
                        f.write(f"{row['smiles']}\t{row['id']}\n")

                return str(smi_path)

        # Fallback to original path
        return str(original_path)

    def update_usage(self, library_id: str):
        """Update usage statistics when library is used"""
        metadata = self._load_metadata()

        if library_id in metadata:
            metadata[library_id]['last_used'] = datetime.now().isoformat()
            metadata[library_id]['use_count'] = metadata[library_id].get('use_count', 0) + 1
            self._save_metadata(metadata)

    def delete_library(self, library_id: str):
        """Delete a library"""
        metadata = self._load_metadata()

        if library_id not in metadata:
            raise ValueError(f"Library {library_id} not found")

        lib = metadata[library_id]

        # Delete files
        lib_dir = Path(lib['file_path']).parent
        if lib_dir.exists():
            shutil.rmtree(lib_dir)

        # Remove from metadata
        del metadata[library_id]
        self._save_metadata(metadata)

        logger.info(f"✓ Library deleted: {lib['name']}")

    def search_libraries(
        self,
        query: str = "",
        tags: List[str] = None,
        min_size: int = None,
        max_size: int = None
    ) -> List[Dict]:
        """
        Search libraries by name, tags, or size.

        Args:
            query: Search in name/description
            tags: Filter by tags
            min_size: Minimum number of molecules
            max_size: Maximum number of molecules

        Returns:
            Filtered list of libraries
        """
        libraries = self.list_libraries()

        # Filter by query
        if query:
            query = query.lower()
            libraries = [
                lib for lib in libraries
                if query in lib['name'].lower() or query in lib.get('description', '').lower()
            ]

        # Filter by tags
        if tags:
            libraries = [
                lib for lib in libraries
                if any(tag in lib.get('tags', []) for tag in tags)
            ]

        # Filter by size
        if min_size is not None:
            libraries = [
                lib for lib in libraries
                if lib['num_molecules'] >= min_size
            ]

        if max_size is not None:
            libraries = [
                lib for lib in libraries
                if lib['num_molecules'] <= max_size
            ]

        return libraries

    def get_statistics(self) -> Dict:
        """Get overall library statistics"""
        libraries = self.list_libraries()

        total_molecules = sum(lib['num_molecules'] for lib in libraries)

        most_used = None
        largest = None

        if libraries:
            most_used = max(libraries, key=lambda x: x.get('use_count', 0))['name']
            largest = max(libraries, key=lambda x: x['num_molecules'])['name']

        return {
            'total_libraries': len(libraries),
            'total_molecules': total_molecules,
            'most_used': most_used,
            'largest': largest
        }
