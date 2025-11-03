"""
ML-ready feature table utilities for the SRA to Features Pipeline.
"""


from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, auc, confusion_matrix, f1_score, make_scorer,
    precision_score, recall_score, roc_auc_score, roc_curve,
    ConfusionMatrixDisplay, RocCurveDisplay
)
from sklearn.model_selection import (
    train_test_split, GridSearchCV, StratifiedKFold
)
from sklearn.preprocessing import (
    RobustScaler, LabelBinarizer, LabelEncoder
)
from typing import List, Dict, Any
import itertools
import json
import logging
import matplotlib
matplotlib.use('Agg')
# Supress font messages
logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import seaborn as sns
import structlog

from ..models.features import FeatureSet


class FeatureTable:
    """Manage pipeline features and convert to tabular format."""

    def __init__(self, logger: structlog.BoundLogger):
        """
        Initialize feature table converter.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger
        self.feature_dicts: List[Dict] = []
    
    def add_feature_dict(self, feature_dict: Dict):
        """Add a feature set to the collection."""
        self.feature_dicts.append(feature_dict)
        self.logger.info("Added feature set", 
                         sra_id=feature_dict['sra_id'])
        return self

    def add_feature_dicts_from_directory(self, directory: Path):
        """
        Load feature dicts from a directory containing 
        features.json files.
        """
        self.logger.info("Loading feature dicts from directory", 
                         directory=str(directory))
        
        # Find all features.json files
        feature_files = list(directory.rglob("features.json"))
        
        for feature_file in feature_files:
            try:
                with open(feature_file, 'r') as f:
                    data = json.load(f)
                
                feature_dict = self._process_features(data)
                self.add_feature_dict(feature_dict)
                
            except Exception as e:
                self.logger.warning("Failed to load feature set", 
                                  file=str(feature_file), error=str(e))
        return self

    def create_sample_features_table(self) -> pd.DataFrame:
        """
        Create a samples × features table for ML training.
        
        Returns:
            DataFrame with samples as rows and features as columns
        """
        if not self.feature_dicts:
            raise ValueError("No feature dicts available")
        
        self.logger.info("Creating ML feature table", 
                         sample_count=len(self.feature_dicts))
        
        # Convert to DataFrame
        df = pd.DataFrame(self.feature_dicts)

        # Transpose table and set indices
        df = df.T
        df.columns = df.iloc[0]
        df = df[1:]

        # Remove rows that sum 0
        df = df[df.sum(axis=1)!=0]
        # Remove columns that sum 0
        df = df.loc[:, (df.sum(axis=0)!=0)]

        # Add 0 values to cnv features
        df.loc[df.index.astype(str).str.startswith("cnv_length")] = \
            df.loc[df.index.astype(str).str.startswith("cnv_length")].fillna(0)

        # Drop NA values so all samples have the same features
        df = df.dropna()

        # Reorder by columns
        df = df.sort_index(axis=1)

        # Reset index to make feature name a column
        df = df.reset_index()

        self.logger.info("ML feature table created", 
                        shape=df.shape, 
                        feature_count=len(df.columns))
        
        return df

    def get_feature_summary(self) -> Dict[str, Any]:
        """Get summary statistics of the feature table."""
        df = self.create_sample_features_table()
        
        return {
            'sample_count': len(df),
            'feature_count': len(df.columns),
            'numeric_features': len(df.select_dtypes(include=[np.number]).columns),
            'categorical_features': len(df.select_dtypes(include=['object']).columns),
            'missing_values': df.isnull().sum().sum(),
            'feature_names': list(df.columns),
        }

    def save_feature_table(self, output_path: Path, format: str='csv'):
        """
        Save the feature table to file.
        
        Args:
            output_path: Output file path
            format: Output format ('csv', 'tsv', 'parquet', 'json')
        """
        # Get sample_features dataframe
        df: pd.DataFrame = self.create_sample_features_table()

        # Define if saving index is needed
        save_index = False

        # Make sure output_path exists
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Save in required format
        if format.lower() == 'csv':
            df.to_csv(output_path, index=save_index)
        elif format.lower() == 'tsv':
            df.to_csv(output_path, index=save_index, sep='\t')
        elif format.lower() == 'parquet':
            df.to_parquet(output_path, index=save_index)
        elif format.lower() == 'json':
            df.to_json(output_path, orient='records', indent=2)
        else:
            raise ValueError(f"Unsupported format: {format}")
        
        self.logger.info("Feature table saved", 
                        path=str(output_path), 
                        format=format,
                        shape=df.shape)

        return self
    
    def _process_features(self, data) -> Dict:
        """
        Process raw feature data into a structured dictionary.
        
        Args:
            data: Raw feature data from JSON file
        
        Returns:
            Processed feature dictionary
        """
        # Initialize feature_dict
        feature_dict = {
            'sra_id': data.get('sra_id', '')
        }
        # Get different sub-elements from data
        genomic_bins: List[Dict] = data.get('genomic_bins', [])
        gene_stats: List[Dict] = data.get('gene_stats', [])
        cnv_regions: List[Dict] = data.get('cnv_regions', [])
        fragment_stats: Dict = data.get('fragment_stats', {})

        # Go through genomic_bins
        for i in range(len(genomic_bins)):
            curr_data = genomic_bins[i]
            chr_n = curr_data['chromosome']
            start = curr_data['start']
            end = curr_data['end']
            # Define key and value
            key = f'bin_gvs_{chr_n}:{start}-{end}'
            val = curr_data['variant_count']
            feature_dict[key] = int(val)
        
        # Go through gene_stats
        for i in range(len(gene_stats)):
            curr_data = gene_stats[i]
            gene_name = curr_data['gene_name'].split(':')[-1]
            start = curr_data['start']
            end = curr_data['end']
            total_gv = int(curr_data['total_variants'])
            dn_ds = curr_data['dn_ds_ratio']
            # Define key and value for GVs
            key_gv = f'gene_gvs_{gene_name}'
            feature_dict[key_gv] = int(total_gv)
            # Check if there are GVs to calculate dn/ds
            if total_gv!=0 and dn_ds is not None:
                key_dn_ds = f'gene_dn_ds_{gene_name}'
                feature_dict[key_dn_ds] = float(dn_ds)
        
        # Go through cnv_regions
        l_cnv_keys = []
        for i in range(len(cnv_regions)):
            curr_data = cnv_regions[i]
            chr_n = curr_data['chromosome']
            cnv_len = abs(int(curr_data['end']) - int(curr_data['start']))
            confidence = curr_data['confidence']
            copy_number = curr_data['copy_number']
            cnv_type = curr_data['type']
            # Filter by confidence
            if confidence > 0.999999:
                # Define key
                key = f'cnv_length_{cnv_type}_chr{chr_n}'
                # Make sure that key exists
                if not (key in feature_dict.keys()):
                    l_cnv_keys.append(str(key))
                    feature_dict[key] = 0
                # Define value to add
                if cnv_type=='gain':
                    copy_add = round(copy_number - 1)
                    val = cnv_len * max(copy_add, 1)
                elif cnv_type=='loss':
                    copy_loss = round(2.0 - (copy_number*2))
                    val = cnv_len * max(copy_loss, 1)
                else:
                    raise ValueError(f'Unrecognised cnv_type: {cnv_type}')
                # Add value to feature_dict[key]
                feature_dict[key] += val
        
        # Go through fragment_stats
        for key, val in fragment_stats.items():
            if key in ['max', 'min', 'mean', 'median', 'std']:
                feature_dict[f'fragment_length_{key}'] = float(val)
            else:
                feature_dict[f'fragment_{key}'] = float(val)
        
        return feature_dict


class MLFeatureTable:
    """Convert pipeline features to ML-ready tabular format."""
    
    def __init__(self, logger: structlog.BoundLogger):
        """
        Initialize ML feature table converter.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger
        self.feature_sets: List[FeatureSet] = []
    
    def add_feature_set(self, feature_set: FeatureSet):
        """Add a feature set to the collection."""
        self.feature_sets.append(feature_set)
        self.logger.info("Added feature set", sra_id=feature_set.sra_id)
    
    def add_feature_sets_from_directory(self, directory: Path):
        """Load feature sets from a directory containing features.json files."""
        self.logger.info("Loading feature sets from directory", directory=str(directory))
        
        # Find all features.json files
        feature_files = list(directory.rglob("features.json"))
        
        for feature_file in feature_files:
            try:
                with open(feature_file, 'r') as f:
                    data = json.load(f)
                
                feature_set = FeatureSet(**data)
                self.add_feature_set(feature_set)
                
            except Exception as e:
                self.logger.warning("Failed to load feature set", 
                                  file=str(feature_file), error=str(e))
    
    def create_sample_features_table(self) -> pd.DataFrame:
        """
        Create a samples × features table for ML training.
        
        Returns:
            DataFrame with samples as rows and features as columns
        """
        if not self.feature_sets:
            raise ValueError("No feature sets available")
        
        self.logger.info("Creating ML feature table", sample_count=len(self.feature_sets))
        
        # Initialize feature dictionary
        features_dict = {}
        
        for feature_set in self.feature_sets:
            sample_id = feature_set.sra_id
            
            # Basic sample features
            sample_features = {
                'sample_id': sample_id,
                'sample_name': feature_set.sample_name or sample_id,
                'pipeline_version': feature_set.pipeline_version,
                'processing_time': feature_set.processing_time,
            }
            
            # Quality metrics features
            qm = feature_set.quality_metrics
            quality_features = {
                'total_reads': qm.total_reads,
                'mapped_reads': qm.mapped_reads,
                'mapping_rate': qm.mapping_rate,
                'mean_coverage': qm.mean_coverage,
                'coverage_std': qm.coverage_std,
                'gc_content': qm.gc_content,
                'duplication_rate': qm.duplication_rate,
            }
            
            # Add quality scores
            for score_name, score_value in qm.quality_scores.items():
                quality_features[f'quality_score_{score_name}'] = score_value
            
            # Fragment length features
            fragment_features = {}
            if feature_set.fragment_stats:
                fs = feature_set.fragment_stats
                fragment_features = {
                    'fragment_mean': fs.mean,
                    'fragment_median': fs.median,
                    'fragment_std': fs.std,
                    'fragment_min': fs.min,
                    'fragment_max': fs.max,
                    'fragment_count': fs.count,
                }
            
            # Genomic variant features
            genomic_features = self._extract_genomic_features(feature_set.genomic_bins)
            
            # Gene-level features
            gene_features = self._extract_gene_features(feature_set.gene_stats)
            
            # CNV features
            cnv_features = self._extract_cnv_features(feature_set.cnv_regions)
            
            # Combine all features
            all_features = {
                **sample_features,
                **quality_features,
                **fragment_features,
                **genomic_features,
                **gene_features,
                **cnv_features,
            }
            
            features_dict[sample_id] = all_features
        
        # Convert to DataFrame
        df = pd.DataFrame.from_dict(features_dict, orient='index')
        
        # Reset index to make sample_id a column
        df = df.reset_index(drop=True)
        
        self.logger.info("ML feature table created", 
                        shape=df.shape, 
                        feature_count=len(df.columns))
        
        return df
    
    def _extract_genomic_features(
        self,
        genomic_bins: List
    ) -> Dict[str, Any]:
        """Extract genomic variant features."""
        if not genomic_bins:
            return {}
        
        # Calculate summary statistics
        variant_counts = [bin.variant_count for bin in genomic_bins]
        
        return {
            'total_genomic_variants': sum(variant_counts),
            'mean_variants_per_bin': np.mean(variant_counts),
            'std_variants_per_bin': np.std(variant_counts),
            'max_variants_per_bin': max(variant_counts),
            'min_variants_per_bin': min(variant_counts),
            'bins_with_variants': sum(
                1 for count in variant_counts if count > 0
            ),
            'total_genomic_bins': len(genomic_bins),
            'variant_density': sum(variant_counts) / len(genomic_bins),
        }
    
    def _extract_gene_features(
        self,
        gene_stats: List
    ) -> Dict[str, Any]:
        """Extract gene-level variant features."""
        if not gene_stats:
            return {}
        
        # Calculate summary statistics
        total_variants = [gene.total_variants for gene in gene_stats]
        synonymous_variants = [gene.synonymous_variants \
                               for gene in gene_stats]
        nonsynonymous_variants = [gene.nonsynonymous_variants \
                                  for gene in gene_stats]
        dn_ds_ratios = [gene.dn_ds_ratio for gene in gene_stats \
                        if gene.dn_ds_ratio is not None]
        
        return {
            'total_genes_with_variants': len(gene_stats),
            'total_gene_variants': sum(total_variants),
            'mean_variants_per_gene': np.mean(total_variants),
            'std_variants_per_gene': np.std(total_variants),
            'total_synonymous_variants': sum(synonymous_variants),
            'total_nonsynonymous_variants': sum(nonsynonymous_variants),
            'mean_dn_ds_ratio': np.mean(dn_ds_ratios) \
                if dn_ds_ratios else 0,
            'std_dn_ds_ratio': np.std(dn_ds_ratios) \
                if dn_ds_ratios else 0,
            'genes_with_high_variants': sum(
                1 for count in total_variants if count > 10
            ),
        }
    
    def _extract_cnv_features(
        self,
        cnv_regions: List
    ) -> Dict[str, Any]:
        """Extract CNV features."""
        if not cnv_regions:
            return {}
        
        # Calculate summary statistics
        copy_numbers = [region.copy_number for region in cnv_regions]
        confidences = [region.confidence for region in cnv_regions]
        gains = sum(1 for region in cnv_regions \
                    if region.type == 'gain')
        losses = sum(1 for region in cnv_regions \
                     if region.type == 'loss')
        
        return {
            'total_cnv_regions': len(cnv_regions),
            'cnv_gains': gains,
            'cnv_losses': losses,
            'mean_copy_number': np.mean(copy_numbers),
            'std_copy_number': np.std(copy_numbers),
            'mean_cnv_confidence': np.mean(confidences),
            'high_confidence_cnvs': sum(
                1 for conf in confidences if conf > 0.9
            ),
            'cnv_burden': len(cnv_regions),
        }
    
    def save_feature_table(
        self,
        output_path: Path,
        format: str = 'csv'
    ):
        """
        Save the feature table to file.
        
        Args:
            output_path: Output file path
            format: Output format ('csv', 'tsv', 'parquet', 'json')
        """
        df = self.create_sample_features_table()
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if format.lower() == 'csv':
            df.to_csv(output_path, index=False)
        elif format.lower() == 'tsv':
            df.to_csv(output_path, index=False, sep='\t')
        elif format.lower() == 'parquet':
            df.to_parquet(output_path, index=False)
        elif format.lower() == 'json':
            df.to_json(output_path, orient='records', indent=2)
        else:
            raise ValueError(f"Unsupported format: {format}")
        
        self.logger.info("Feature table saved", 
                        path=str(output_path), 
                        format=format,
                        shape=df.shape)
    
    def get_feature_summary(self) -> Dict[str, Any]:
        """Get summary statistics of the feature table."""
        df = self.create_sample_features_table()
        
        return {
            'sample_count': len(df),
            'feature_count': len(df.columns),
            'numeric_features': len(df.select_dtypes(include=[np.number]).columns),
            'categorical_features': len(df.select_dtypes(include=['object']).columns),
            'missing_values': df.isnull().sum().sum(),
            'feature_names': list(df.columns),
        }


def analyze_feature_groups_cvrf(
    df: pd.DataFrame,
    output_folder: Path,
    features_to_analyze: List[str],
    feature_n: int,
    target_column: str,
    logger: structlog.BoundLogger,
    split_n: int=5,
    rand_seed: int=None
) -> pd.DataFrame:
    """
    Performs a classification analysis for all combinations of N
    features from a given list and ranks them by AUC score.

    Performs cross-validated Random Forest classification.

    Args:
        df (pd.DataFrame): The DataFrame containing features and the 
                           target.
        output_folder (Path): Folder where output files are created.
        features_to_analyze (list): A list of feature names to be 
                                    grouped.
        feature_n (int): Number of features to group together.
        target_column (str): The name of the target variable column.
        logger: Logger instance.
        split_n (int): Number of splits for cross-validation.
        rand_seed (int): Random seed for repeatability.
    
    Returns:
        pd.DataFrame: A DataFrame of feature pairs grouped by AUC.
    """
    if len(features_to_analyze) < feature_n:
        err = f"Please provide at least {feature_n} "+\
            "features for grouping."
        logger.error(err, feature_n=feature_n)
        return pd.DataFrame()

    # Make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)
    # Define an output folder for models
    models_out = output_folder / 'saved_models'
    # Make sure models folder exists
    models_out.mkdir(parents=True, exist_ok=True)

    # Define X and y
    X = df[features_to_analyze]
    y = df[target_column]

    # Generate all unique combinations of the selected features
    feature_groups = list(itertools.combinations(
        features_to_analyze, feature_n
    ))

    logger.info(
        str(f"Generated {len(feature_groups)} unique feature groups "+\
            f"from the list of {len(features_to_analyze)} features."),
        feature_n=feature_n
    )

    # Initialize results dictionary
    results = {}

    # Initialize CV and Model
    cv = StratifiedKFold(n_splits=split_n, shuffle=True, 
                         random_state=rand_seed)
    rf_base_model = RandomForestClassifier(
        n_estimators=100,
        max_depth=5,
        class_weight='balanced',
        random_state=rand_seed,
        n_jobs=-1
    )
    # Initialize the label encoder and encode y
    label_encoder = LabelEncoder()
    y_encoded = label_encoder.fit_transform(y)
    n_classes = len(label_encoder.classes_)

    ### Display
    n_feat_pairs = len(feature_groups)
    display_range = int((10**max(2, feature_n)) / 4)
    ###

    logger.debug(f'display_range set to {display_range}',
                 feature_n=feature_n)

    # Initialize max auc score
    max_auc_score = 0

    for i, curr_group in enumerate(feature_groups):
        # Define a feature group key
        feature_group_key = f'ID_{i}'
        # Create a new DataFrame with only two features for training
        X_group = X[list(curr_group)].values

        # Initialize the list of model names
        l_model_files: List[Path] = []

        # Initialize AUC scores list and used model
        fold_auc_scores = []
        rf_model = rf_base_model

        # Perform Cross-Validated Random Forest
        for fold, (train_index, test_index) in \
            enumerate(cv.split(X_group, y_encoded)):
            X_train, X_test = X_group[train_index], \
                X_group[test_index]
            y_train, y_test = y_encoded[train_index], \
                y_encoded[test_index]

            try:
                rf_model.fit(X_train, y_train)
                y_proba = rf_model.predict_proba(X_test)
                
                auc_score = roc_auc_score(
                    y_test, 
                    y_proba, 
                    multi_class='ovr', 
                    average='macro',
                    labels=np.arange(n_classes)
                )
                fold_auc_scores.append(auc_score)

                # Pickle the model
                model_path = models_out / \
                    f'cvrf_model_{feature_group_key}_fold{fold+1}.pickle'
                pickle_model(rf_model, Path(model_path))
                # Append model name to l_model_files
                l_model_files.append(model_path)
                
            except ValueError as e:
                logger.warn(
                    f"Skipping fold {fold} for group "+\
                    f"{feature_group_key} due to error: {e}"
                )

        if fold_auc_scores:
            mean_auc = np.mean(fold_auc_scores)
            results[feature_group_key] = [fi for fi in curr_group] +\
                                         [mean_auc]
        else:
            results[feature_group_key] = [fi for fi in curr_group] +\
                                         [np.nan]

        # Check if mean_auc is bigger than max_auc_score
        if mean_auc > max_auc_score:
            # Update max_auc_score
            max_auc_score = mean_auc
        else:
            # Remove all files in l_model_files
            for model_file in l_model_files:
                # Remove using Path functions
                model_file.unlink(missing_ok=True)

        ### Display
        if i==0 or (i+1)%display_range==0:
            logger.debug(f'Progress: {i+1} / {n_feat_pairs}',
                         feature_n=len(features_to_analyze))
        ###

    # Create a DataFrame and sort the results
    feature_cols = [f'Feature {i+1}' for i in range(feature_n)]
    feature_cols.append('AUC Score')
    group_auc_df = pd.DataFrame.from_dict(
        results, orient='index', columns=feature_cols
    )
    group_auc_df.index.name = 'Feature Group ID'
    group_auc_df = group_auc_df.sort_values(by='AUC Score',
                                            ascending=False)

    return group_auc_df


def analyze_feature_pairs(
    df: pd.DataFrame,
    output_folder: Path,
    features_to_analyze: List[str],
    target_column: str,
    logger: structlog.BoundLogger,
    rand_seed: int=None
) -> pd.DataFrame:
    """
    Performs a classification analysis for all combinations of two 
    features from a given list and ranks them by AUC score.

    Args:
        df (pd.DataFrame): The DataFrame containing features and the 
                           target.
        output_folder (Path): Folder where output files are created.
        features_to_analyze (list): A list of feature names to be 
                                    paired.
        target_column (str): The name of the target variable column.
        logger: Logger instance.
        rand_seed (int): Random seed for repeatability.

    Returns:
        pd.DataFrame: A DataFrame of feature pairs ranked by AUC.
    """
    if len(features_to_analyze) < 2:
        logger.error(
            "Please provide at least two features for pairing.",
            feature_n=len(features_to_analyze)
        )
        return pd.DataFrame()

    # Make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)

    X = df[features_to_analyze]
    y = df[target_column]

    # Generate all unique pairs of the selected features
    feature_pairs = list(itertools.combinations(
        features_to_analyze, 2
    ))

    logger.info(
        str(f"Generated {len(feature_pairs)} unique feature pairs "+\
            f"from the list of {len(features_to_analyze)} features."),
        feature_n=len(features_to_analyze)
    )

    pair_results = {}
    model = RandomForestClassifier(random_state=rand_seed)
    lb = LabelBinarizer()

    # Counter
    cont = 0
    ### Display
    n_feat_pairs = len(feature_pairs)
    ###

    for f1, f2 in feature_pairs:
        # Create a new DataFrame with only two features for training
        X_pair = X[[f1, f2]]

        X_train, X_test, y_train, y_test = train_test_split(
            X_pair,
            y,
            test_size=0.3,
            random_state=rand_seed,
            stratify=y
        )

        # Save train and test sets
        X_train.to_csv(f'{output_folder}/X_train_{cont}.csv')
        X_test.to_csv(f'{output_folder}/X_test_{cont}.csv')
        y_train.to_csv(f'{output_folder}/y_train_{cont}.csv')
        y_test.to_csv(f'{output_folder}/y_test_{cont}.csv')

        # Train a model
        trained_model = model.fit(X_train, y_train)

        # Pickle trained model
        pickle_model(trained_model,
                     f'{output_folder}/model_pair_ID_{cont}.pickle')

        # Handle multi-class classification for AUC calculation
        y_test_bin = lb.fit_transform(y_test)
        y_pred_proba = trained_model.predict_proba(X_test)

        if y_test_bin.shape[1] == 1:
            auc_score = roc_auc_score(y_test_bin, y_pred_proba[:, 1])
        else:
            auc_score = roc_auc_score(y_test_bin, y_pred_proba,
                                      multi_class='ovr')

        # Load the data
        pair_results[f'ID_{cont}'] = (f1, f2, auc_score)

        ### Display
        if cont==0 or (cont+1)%10==0:
            logger.debug(f'Progress: {cont+1} / {n_feat_pairs}',
                         feature_n=len(features_to_analyze))
        ###

        cont += 1
    
    # Create a DataFrame and sort the results
    pair_auc_df = pd.DataFrame.from_dict(
        pair_results, orient='index', columns=[
            'Feature 1', 'Feature 2', 'AUC Score'
        ]
    )
    pair_auc_df.index.name = 'Feature Pair ID'
    pair_auc_df = pair_auc_df.sort_values(by='AUC Score',
                                          ascending=False)

    return pair_auc_df


def analyze_feature_trios(
    df: pd.DataFrame,
    output_folder: Path,
    features_to_analyze: List[str],
    target_column: str,
    logger: structlog.BoundLogger,
    rand_seed: int=None
) -> pd.DataFrame:
    """
    Performs a classification analysis for all combinations of three 
    features from a given list and ranks them by AUC score.

    Args:
        df (pd.DataFrame): The DataFrame containing features and the 
                           target.
        output_folder (Path): Folder where output files are created.
        features_to_analyze (list): A list of feature names to be 
                                    grouped.
        target_column (str): The name of the target variable column.
        logger: Logger instance.
        rand_seed (int): Random seed for repeatability.

    Returns:
        pd.DataFrame: A DataFrame of feature pairs ranked by AUC.
    """
    if len(features_to_analyze) < 3:
        logger.error(
            "Please provide at least three features for pairing.",
            feature_n=len(features_to_analyze)
        )
        return pd.DataFrame()

    # Make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)

    X = df[features_to_analyze]
    y = df[target_column]

    # Generate all unique trios of the selected features
    feature_trios = list(itertools.combinations(
        features_to_analyze, 3
    ))

    logger.info(
        str(f"Generated {len(feature_trios)} unique feature trios "+\
            f"from the list of {len(features_to_analyze)} features."),
        feature_n=len(features_to_analyze)
    )

    trio_results = {}
    model = RandomForestClassifier(random_state=rand_seed)
    lb = LabelBinarizer()

    # Counter
    cont = 0
    ### Display
    n_feat_trios = len(feature_trios)
    ###

    for f1, f2, f3 in feature_trios:
        # Create a new DataFrame with only the features for training
        X_pair = X[[f1, f2, f3]]

        X_train, X_test, y_train, y_test = train_test_split(
            X_pair,
            y,
            test_size=0.3,
            random_state=rand_seed,
            stratify=y
        )

        # Save train and test sets
        X_train.to_csv(f'{output_folder}/X_train_{cont}.csv')
        X_test.to_csv(f'{output_folder}/X_test_{cont}.csv')
        y_train.to_csv(f'{output_folder}/y_train_{cont}.csv')
        y_test.to_csv(f'{output_folder}/y_test_{cont}.csv')

        # Train a model
        trained_model = model.fit(X_train, y_train)

        # Pickle trained model
        pickle_model(trained_model,
                     f'{output_folder}/model_trio_ID_{cont}.pickle')

        # Handle multi-class classification for AUC calculation
        y_test_bin = lb.fit_transform(y_test)
        y_pred_proba = trained_model.predict_proba(X_test)

        if y_test_bin.shape[1] == 1:
            auc_score = roc_auc_score(y_test_bin, y_pred_proba[:, 1])
        else:
            auc_score = roc_auc_score(y_test_bin, y_pred_proba, 
                                      multi_class='ovr')

        # Load the data
        trio_results[f'ID_{cont}'] = (f1, f2, f3, auc_score)

        ### Display
        if cont==0 or (cont+1)%100==0:
            logger.debug(f'Progress: {cont+1} / {n_feat_trios}',
                         feature_n=len(features_to_analyze))
        ###

        cont += 1
    
    # Create a DataFrame and sort the results
    trio_auc_df = pd.DataFrame.from_dict(
        trio_results, orient='index', columns=[
            'Feature 1', 'Feature 2', 'Feature 3', 'AUC Score'
        ]
    )
    trio_auc_df.index.name = 'Feature Trio ID'
    trio_auc_df = trio_auc_df.sort_values(by='AUC Score', 
                                          ascending=False)

    return trio_auc_df


def analyze_features_and_rank_by_auc(
    df: pd.DataFrame,
    output_folder: Path,
    target_var: str,
    model,
    logger: structlog.BoundLogger,
    rand_seed:int=None
) -> pd.DataFrame:
    """
    Performs feature-by-feature analysis and ranks them by AUC score.

    Args:
        df (pd.DataFrame): The DataFrame with features and a 
                           target variable column.
        output_folder (Path): Folder where output files are created.
        target_var (str): Target variable to be predicted with the 
                          feature classification.
        model: A sklearn model instance.
        logger: Logger instance.
        rand_seed (int): Random seed for repeatability.

    Returns:
        pd.DataFrame: A DataFrame of features ranked by AUC.
    """
    # Obtain X and y
    X = df.drop(target_var, axis=1)
    y = df[target_var]
    # Initialize results dictionary
    results = {}
    # Define a folder for feature training and testing sets
    out_sets = output_folder / 'feature_sets'
    # Make sure the folder exists
    out_sets.mkdir(parents=True, exist_ok=True)
    ### Display
    # Counter
    cont = 0
    n_col = len(X.columns)
    ###
    # Initialize max_auc_score with 0
    max_auc_score = 0.0
    for feature in X.columns:
        auc_score = get_single_feature_auc(X, y, out_sets, feature,
                                           model, rand_seed,
                                           max_auc_score)
        # Update max_auc_score if needed
        if auc_score > max_auc_score:
            max_auc_score = auc_score
        # Store the result
        results[feature] = auc_score
        ### Display
        if cont==0 or (cont+1)%1000==0:
            logger.debug(f'Progress: {cont+1} / {n_col}')
        cont += 1
        ###
    
    # Create a DataFrame to store and sort the results
    auc_df = pd.DataFrame.from_dict(
        results, orient='index', columns=['AUC Score']
    ).sort_values(by='AUC Score', ascending=False)
    
    return auc_df


def cleanup_empty_rows_cols(df: pd.DataFrame):
    """
    Performs cleanup of a DataFrame by removing empty rows
    and columns.
    """
    # Replace empty strings with NaN
    df_out = df.replace('', np.nan)
    # Remove columns with NAN
    df_out = df_out.dropna(axis=1, how='any')
    # Remove rows that sum 0
    df_out = df_out[df.sum(axis=1)!=0]
    # Remove columns that sum 0
    df_out = df_out.loc[:, (df_out.sum(axis=0)!=0)]

    return df_out


def collapse_probabilities(y_pred: np.ndarray) -> np.ndarray:
    """
    Collapse probability predictions into class labels.
    """
    # Select maximum value per row and set it to 1, the rest to 0
    for i in range(len(y_pred)):
        curr_row = y_pred[i]
        max_val = 0
        for j in range(len(curr_row)):
            curr_cell = curr_row[j]
            if curr_cell>max_val:
                max_val = curr_cell
                max_col = int(j)
        for j in range(len(curr_row)):
            if j==max_col:
                y_pred[i][j] = 1
            else:
                y_pred[i][j] = 0
    # Set type to int
    y_pred = y_pred.astype(int)
    return y_pred


def create_auc_roc_curves(
    y_true_bin: np.ndarray,
    y_pred: np.ndarray,
    out_folder: Path,
    bname: str,
    out_name: str,
    model,
    logger: structlog.BoundLogger
) -> None:
    """
    Creates a ROC curve plot for multi-class classification.
    """

    logger.info('Calculating micro-averaged AUC score for OvR strategy.')

    # Prepare for plotting the one-vs-rest ROC curves
    n_classes = len(model.classes_)
    
    # Set up the plot
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    
    # Plot chance level (diagonal line)
    plt.plot([0, 1], [0, 1], linestyle="--", lw=2, color="navy", 
             label="Chance level", alpha=0.8)

    # Loop through each class to plot its individual ROC curve
    for i in range(n_classes):
        # Create a display for each class using RocCurveDisplay
        RocCurveDisplay.from_predictions(
            y_true_bin[:, i],
            y_pred[:, i],
            name=f"ROC curve (class {i})",
            ax=ax,
        )

    out_name_spaces = out_name.replace('_', ' ').title()
    title = 'Receiver Operating Characteristic (ROC) Curve '+\
        f'for {out_name_spaces}'
    plt.title(title)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.savefig(out_folder / f'roc_curve_{bname}{out_name}.png')


def create_conf_matrix(y_true, y_pred, classes, out_path):
    # Define integer labels
    int_labels = np.arange(len(classes))
    # Plot confusion matrix
    cmb = confusion_matrix(y_true, y_pred, 
                           labels=int_labels)

    fig, a0 = plt.subplots(figsize=(8,4))
    dispb = ConfusionMatrixDisplay(confusion_matrix=cmb, 
                                   display_labels=classes)
    out_name = out_path.stem.replace('_', ' ').title()
    a0.set_title(out_name)
    dispb.plot(ax=a0,colorbar=False)
    plt.tight_layout()
    fig.savefig(out_path)


def create_feature_table_from_directory(
    input_directory: Path,
    output_path: Path,
    logger: structlog.BoundLogger,
    format: str = 'csv'
) -> pd.DataFrame:
    """
    Create a feature table from a pipeline output directory.
    
    Args:
        input_directory: Directory containing pipeline results
        output_path: Output folder path for feature table
        format: Output format ('csv', 'tsv', 'parquet', 'json')
        logger: Logger instance
        
    Returns:
        DataFrame with samples × features
    """
    ml_table = FeatureTable(logger)
    ml_table.add_feature_dicts_from_directory(input_directory)
    
    if not ml_table.feature_dicts:
        raise ValueError(f"No feature sets found in {input_directory}")
    
    # Make sure output path is a directory
    output_path.mkdir(parents=True, exist_ok=True)

    # Define output file path
    output_file = output_path / f'features_table.{format}'

    # Save feature table
    ml_table.save_feature_table(output_file, format)

    # Log summary
    summary = ml_table.get_feature_summary()
    logger.info("ML feature table created successfully", **summary)
    
    return ml_table.create_sample_features_table()


def create_heatmap(
    df:pd.DataFrame,
    title:str,
    save_fig:bool,
    out_path:str,
    logger: structlog.BoundLogger
) -> None:
    """Creates a heatmap using Seaborn and Matplotlib"""
    # Check that df is not empty
    if df.empty:
        logger.warning(
            f'Heatmap with title "{title}" could not be created.'
        )
        return None
    logger.info(f'Creating heatmap with title "{title}"...',
                title=title)
    # Define the figure size
    fig_width = 9
    fig_height = 9
    # Set the font scale to make labels readable
    sns.set_theme(font_scale=0.8)
    # Create the figure and axes objects
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    # Load the dataframe
    sns.heatmap(df, cmap="viridis", ax=ax, vmin=0) # vmin=-20, vmax=20
    # Add titles and formatting
    plt.title(f"Features Heatmap for {title}")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    # Show or save the plot
    if save_fig:
        path_out = os.path.join(out_path, f'heatmap_{title}.png')
        plt.savefig(path_out)
    return None


def create_ml_feature_table_from_directory(
    input_directory: Path,
    output_path: Path,
    logger: structlog.BoundLogger,
    format: str = 'csv'
) -> pd.DataFrame:
    """
    Create ML-ready feature table from pipeline output directory.
    
    Args:
        input_directory: Directory containing pipeline results
        output_path: Output file path for feature table
        format: Output format ('csv', 'tsv', 'parquet', 'json')
        logger: Logger instance
        
    Returns:
        DataFrame with samples × features
    """
    ml_table = MLFeatureTable(logger)
    ml_table.add_feature_sets_from_directory(input_directory)
    
    if not ml_table.feature_sets:
        raise ValueError(f"No feature sets found in {input_directory}")
    
    # Save feature table
    ml_table.save_feature_table(output_path, format)
    
    # Log summary
    summary = ml_table.get_feature_summary()
    logger.info("ML feature table created successfully", **summary)
    
    return ml_table.create_sample_features_table()


def create_macro_avg_auc_roc_curve(
    y_true_bin: np.ndarray,
    y_pred: np.ndarray,
    out_folder: Path,
    bname: str,
    out_name: str,
    logger: structlog.BoundLogger
) -> None:
    """
    Creates a macro-averaged ROC curve plot for multi-class 
    classification.
    
    Args:
        y_true_bin (np.ndarray): True binary labels (one-hot encoded).
        y_pred (np.ndarray): Predicted probabilities 
                             (n_samples, n_classes).
        out_folder (Path): Folder to save the plot.
        bname (str): Base name for the output file.
        out_name (str): Additional name component for the output and 
                        plot title.
        logger (structlog.BoundLogger): Logger instance.
    """

    logger.info('Creating macro-averaged ROC curve.')
    
    # Determine the number of classes
    n_classes = y_true_bin.shape[1]

    # Dictionaries to store per-class ROC components
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # Calculate ROC curve and AUC for each class
    for i in range(n_classes):
        # Calculate the ROC curve for class 'i' as the positive class 
        # against all others
        fpr[i], tpr[i], _ = roc_curve(y_true_bin[:, i], y_pred[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Aggregate all false positive rates (FPRs)
    all_fpr = np.unique(np.concatenate(
        [fpr[i] for i in range(n_classes)]
    ))

    # Interpolate all ROC curves at these common FPR points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        # Interpolate the TPR for class 'i' to match the 'all_fpr' points
        mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])

    # Average the interpolated TPRs
    mean_tpr /= n_classes

    # Store the final macro-average results
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    macro_roc_auc = auc(fpr["macro"], tpr["macro"])

    # Generate and save the AUC graph
    plt.figure(figsize=(8, 6))

    # Plot the macro-average curve
    label_str = f'Macro-average ROC curve (AUC = {macro_roc_auc:.2f})'
    plt.plot(fpr["macro"], tpr["macro"], label=label_str, color='blue',
             linewidth=2)
             
    # Plot the chance level line
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--',
             label='Chance level')

    # Define title of the graph
    out_name_spaces = out_name.replace('_', ' ').lstrip().title()
    title = 'Macro-averaged Receiver Operating Characteristic (ROC) ' + \
        f'Curve for {out_name_spaces}'
    # Set title and labels of the graph
    plt.title(title)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.savefig(out_folder / f'macro_avg_roc_curve_{bname}{out_name}.png')


def create_micro_avg_auc_roc_curve(
    y_true_bin: np.ndarray,
    y_pred: np.ndarray,
    out_folder: Path,
    bname: str,
    out_name: str,
    logger: structlog.BoundLogger
) -> None:
    """
    Creates a micro-averaged ROC curve plot for multi-class 
    classification.

    Args:
        y_true_bin (np.ndarray): True binary labels (one-hot encoded).
        y_pred (np.ndarray): Predicted probabilities 
                             (n_samples, n_classes).
        out_folder (Path): Folder to save the plot.
        bname (str): Base name for the output file.
        out_name (str): Additional name component for the output and 
                        plot title.
        logger (structlog.BoundLogger): Logger instance.
    """

    logger.info('Creating micro-averaged ROC curve.')
    
    # Calculate micro-averaged ROC curve
    fpr, tpr, _ = roc_curve(y_true_bin.ravel(), y_pred.ravel())
    micro_roc_auc = auc(fpr, tpr)

    # Generate and save the AUC graph
    plt.figure(figsize=(8, 6))

    # Plot the micro-average curve
    label_str = f'Micro-average ROC curve (AUC = {micro_roc_auc:.2f})'
    plt.plot(fpr, tpr, color='blue', linewidth=2, label=label_str)

    # Plot the chance level line
    plt.plot(np.arange(0, 1.1, 0.1), np.arange(0, 1.1, 0.1),
             color='navy', lw=2, linestyle='--', label='Chance level')

    # Define title of the graph
    out_name_spaces = out_name.replace('_', ' ').lstrip().title()
    title = 'Micro-averaged Receiver Operating Characteristic (ROC) '+\
        f'Curve for {out_name_spaces}'
    # Set title and labels of the graph
    plt.title(title)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.savefig(out_folder / f'micro_avg_roc_curve_{bname}{out_name}.png')


def cross_validate_rf(
    df: pd.DataFrame,
    target_var: str,
    output_folder: Path,
    logger,
    rand_seed: int|None,
    top_n: int=0,
    refit: str='roc_auc'
) -> List[Any]:
    """
    Perform cross-validation with a Random Forest classifier.

    Args:
        df (pd.DataFrame): DataFrame containing features and target 
                           variable.
        target_var (str): Name of the target variable column.
        output_folder (Path): Folder where output files are created.
        logger: Logger instance.
        rand_seed (int): Random seed for repeatability.
        refit (str): Metric to optimize during hyperparameter tuning.
                     Default is 'roc_auc'. 
                     Recommended alternative is 'f1'. 

    Returns:
        float: Mean ROC-AUC score across folds.
        Any: Best model from cross-validation.
        List[str]: List of top N feature names by importance.
    """

    X = df.drop(columns=[target_var])
    y = df[target_var]

    model = RandomForestClassifier(random_state=rand_seed)
    
    # Encode the target variable
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    
    # Store the mapping for reference
    label_mapping = dict(zip(le.classes_, le.transform(le.classes_)))
    logger.info(f"Target label mapping: {label_mapping}")

    logger.info("Splitting data into training and testing sets...")

    # Split the data into training and testing sets
    Xb_train, Xb_test, yb_train, yb_test = train_test_split(
        X, y_encoded, test_size=0.3, random_state=rand_seed
    )

    # Define hyperparameters for CV to search over
    cv_params = {
        'max_depth': [2,3],
        'min_samples_leaf': [1,2,3],
        'min_samples_split': [2,3,4],
        'max_features': [1,2,3,4],
        'n_estimators': [75, 100]
    }
    # Define scoring metrics
    scoring = {
        'accuracy': 'accuracy',
        'precision': make_scorer(precision_score, average='weighted', 
                                 zero_division=0),
        'recall': make_scorer(recall_score, average='weighted', 
                              zero_division=0),
        'f1': make_scorer(f1_score, average='weighted', 
                          zero_division=0),
        'roc_auc': 'roc_auc_ovr' 
                    # make_scorer(roc_auc_score, multi_class='ovr')
    }
    refit = 'roc_auc' # Metric to optimize

    logger.info("Instantiating cross-validated Random Forest model...")

    # Instantiate the cross validated random forest
    model_cv = GridSearchCV(model, cv_params, scoring=scoring,
                            cv=5, refit=refit, n_jobs=4)

    logger.info("Fitting model with cross-validation...")

    # Fit the model
    model_cv.fit(Xb_train, yb_train)

    logger.info("Model fitting completed.")

    logger.info("Getting best model from cross-validation...")

    # Get the best model and best parameters
    best_model = model_cv.best_estimator_
    best_params = model_cv.best_params_

    # Save best model
    pickle_model(best_model, 
                 f'{output_folder}/best_model_single.pickle')

    logger.info("Best model found", best_params=best_params)

    logger.info("Evaluating best model on test set...")

    # Evaluate the best model on the test set
    yb_pred_proba = best_model.predict_proba(Xb_test)

    logger.info("Model evaluation on test set completed.")

    logger.info("Calculating AUC score...")

    # Calculate mean score
    try:
        mean_score = roc_auc_score(yb_test, yb_pred_proba,
                                   multi_class='ovr')
    except Exception as e:
        logger.warning(f"Final ROC-AUC calculation failed: {e}")
        mean_score = 0.5

    logger.info("Cross-validation completed.", mean_auc=mean_score)

    # Get predictions on the test set
    yb_pred = best_model.predict(Xb_test)

    # Plot confusion matrix
    cmb = confusion_matrix(yb_test, yb_pred, 
                           labels=best_model.classes_)

    fig, a0 = plt.subplots(figsize=(8,4))
    dispb = ConfusionMatrixDisplay(confusion_matrix=cmb, 
                                   display_labels=best_model.classes_)
    a0.set_title("Best CVRF model (single feature)")
    dispb.plot(ax=a0,colorbar=False)
    plt.tight_layout()
    fig.savefig(
        f'{output_folder}/confusion_matrix.png'
    )

    importance = best_model.feature_importances_
    forest_importance = pd.Series(importance, 
                                  index=model_cv.feature_names_in_)

    # Get feature names - use the best model's feature names
    if hasattr(best_model, 'feature_names_in_'):
        feature_names = best_model.feature_names_in_
    else:
        # Fallback: use column names from the original DataFrame
        feature_names = X.columns.tolist()

    # Create a Series with feature importances and names
    forest_importance = pd.Series(importance, index=feature_names)

    # Sort by importance (descending)
    importance_sorted = forest_importance.sort_values(
        ascending=False
    )

    # Save the sorted importances to a CSV file
    importance_sorted.to_csv(
        os.path.join(output_folder, f'feature_importances.csv'),
        header=['Importance']
    )

    if top_n >= 2:
        # Select top n features by importance
        top_n_features = importance_sorted.head(top_n)

        # Plot the top n feature importances
        fig, ax = plt.subplots()
        top_n_features.plot.bar(ax=ax)
        ax.set_title(f"Top {top_n} feature importances")
        ax.set_ylabel("Mean decrease in impurity")
        fig.tight_layout()
        fig.savefig(
            f'{output_folder}/top_{top_n}_feature_importances.png'
        )

        # Get just the names of the top n features as a list
        top_n_feature_names = top_n_features.index.tolist()

        # Print or log the results
        print("Top 20 features by importance:")
        for i, feature_name in enumerate(top_n_feature_names, 1):
            print(f"{i}. {feature_name}:"+\
                  f" {top_n_features[feature_name]:.6f}")
    else:
        top_n_feature_names = []

    return mean_score, best_model, top_n_feature_names


def cross_validated_feature_analysis(
    table_name: str,
    output_folder: Path,
    target_var: str,
    logger: structlog.BoundLogger,
    top_n: int=20,
    split_n: int=5,
    rand_seed: int=None
) -> None:
    """
    Performs RandomForest CrossValidated feature analysis.
    """
    # Make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)
    # Define base name for outputs
    bname = str(table_name).rsplit('/')[-1].rsplit('.')[0]

    # Read feature table
    df = pd.read_csv(table_name, index_col=0)
    # Transpose df to have samples as rows and features as columns
    df = df.T

    out_feature_models = output_folder / 'feature_models'

    # Perform CVRF
    ordered_features = top_single_features_cvrf(
        df,
        target_var,
        out_feature_models,
        logger,
        split_n=split_n,
        random_state=rand_seed
    )
    # Save top features and their scores
    ordered_features.to_csv(output_folder / f"{bname}_top_features.csv",
                            header=['AUC Score'])
    
    # Pick the top_n features and get the feature names
    top_features_series = ordered_features.head(top_n)
    top_feature_names = top_features_series.index.tolist()

    # Train final models and store them
    final_models_and_results = train_and_store_final_models(
        df, 
        target_var, 
        output_folder / 'top_models',
        top_feature_names, 
        ordered_features,
        logger,
        rand_seed
    )

    # Generate ROC curves for the top features
    logger.info("Generating ROC curves for top features")
    for feature_name, (cv_auc, fitted_model) in \
        final_models_and_results.items():
        plot_feature_roc(
            df, 
            target_var, 
            output_folder / 'top_models',
            feature_name,
            fitted_model,
            logger,
            random_state = rand_seed
        )

    # Run in pairs/trios
    if top_n >= 2:

        logger.debug(f'Top {top_n} features obtained.',
                     top_features=top_feature_names)

        logger.info(f'Analyzing top {top_n} features in pairs.',
                    table_name=table_name,
                    top_n=top_n)

        # Define pairs output name
        output_folder_pairs = output_folder / "pairs_auc_out"
        output_pairs = output_folder_pairs / \
            f"{bname}_feature_pairs.csv"

        # Run the analysis
        top_feature_pairs_ranked = analyze_feature_groups_cvrf(
            df,
            output_folder_pairs,
            top_feature_names,
            2,
            target_var,
            logger,
            split_n=split_n,
            rand_seed=rand_seed
        )

        # Save the resulting dataframe
        top_feature_pairs_ranked.to_csv(output_pairs)
        
        if top_n >= 3:
            logger.info(f'Analyzing top {top_n} features in trios.',
                    table_name=table_name,
                    top_n=top_n)
            # Define trios output name
            output_folder_trios = output_folder / "trios_auc_out"
            output_trios = output_folder_trios / \
                f"{bname}_feature_trios.csv"

            # Run the analysis for trios
            top_feature_trios_ranked = analyze_feature_groups_cvrf(
                df,
                output_folder_trios,
                top_feature_names,
                3,
                target_var,
                logger,
                split_n=split_n,
                rand_seed=rand_seed
            )

            # Save the resulting dataframe
            top_feature_trios_ranked.to_csv(output_trios)
    return None


def cvrf_load_top_features(
    df: pd.DataFrame,
    target_var: str,
    output_folder: Path,
    top_features_file: Path,
    bname: str,
    group_n: int,
    logger: structlog.BoundLogger,
    top_n: int=20,
    split_n: int=10,
    rand_seed: int=None
) -> None:
    """
    Runs the cross_validated_feature_analysis pipeline by loading
    top features from a previous run.
    """
    # Open top features file to get the ordered feature scores
    ordered_features = pd.read_csv(top_features_file, index_col=0)

    # Pick the top_n features and get the feature names
    top_features_series = ordered_features.head(top_n)
    top_feature_names = top_features_series.index.tolist()

    # Train final models and store them
    final_models_and_results = train_and_store_final_models(
        df, 
        target_var, 
        output_folder / 'top_models',
        top_feature_names, 
        ordered_features,
        logger,
        rand_seed
    )

    # Generate ROC curves for the top features
    logger.info("Generating ROC curves for top features")
    for feature_name, (cv_auc, fitted_model) in \
        final_models_and_results.items():
        plot_feature_roc(
            df, 
            target_var, 
            output_folder / 'top_models',
            feature_name,
            fitted_model,
            logger,
            random_state = rand_seed
        )

    # Run in pairs/trios
    if (group_n > 1) and (top_n >= group_n):

        logger.debug(f'Top {top_n} features obtained.',
                     top_features=top_feature_names)

        logger.info(
            f'Analyzing top {top_n} features in groups of {group_n}.',
            top_n=top_n
        )

        # Define group output name
        output_folder_group = output_folder / \
            f"{group_n}_groups_auc_out"
        output_group = output_folder_group / \
            f"{bname}_feature_{group_n}_groups.csv"

        # Run the analysis
        top_feature_groups_ranked = analyze_feature_groups_cvrf(
            df,
            output_folder_group,
            top_feature_names,
            group_n,
            target_var,
            logger,
            split_n=split_n,
            rand_seed=rand_seed
        )
        # Save the resulting dataframe
        top_feature_groups_ranked.to_csv(output_group)

    return None


def get_single_feature_auc(
    X:pd.DataFrame,
    y:pd.Series,
    output_folder: Path,
    feature_name:str,
    model,
    rand_seed:int,
    max_auc_score:float=0.0
) -> float:
    """
    Trains a model using a single feature and returns the ROC-AUC 
    score.

    Args:
        X (pd.DataFrame): DataFrame of features.
        y (pd.Series): Series of target labels.
        feature_name (str): The name of the feature to use.
        model: The machine learning model to train.
        rand_seed (int): Random seed for repeatability.
        max_auc_score (float): Maximum AUC score observed so far.
                               Used to skip saving files if the new
                               score is not better.

    Returns:
        float: The ROC-AUC score.
    """
    # Get the single feature from X
    X_single = X[[feature_name]]
    # Split train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
        X_single,
        y,
        test_size=0.3,
        random_state=rand_seed,
        stratify=y
    )
    
    # Handle multi-class classification for AUC
    lb = LabelBinarizer()
    y_test_bin = lb.fit_transform(y_test)
    curr_model = model.fit(X_train, y_train)
    y_pred_proba = curr_model.predict_proba(X_test)

    # Calculate AUC score
    if y_test_bin.shape[1] == 1:
        auc_score = roc_auc_score(y_test_bin, y_pred_proba[:, 1])
    else:
        auc_score = roc_auc_score(y_test_bin, y_pred_proba, 
                                  multi_class='ovr')

    # If the AUC score is better than max_auc_score, 
    # delete previous files and save new ones
    if auc_score > max_auc_score:
        # Delete previous files
        for f in os.listdir(output_folder):
            os.remove(os.path.join(output_folder, f))
        # Save train and test sets
        feat_name = feature_name.replace(' ', '_')
        X_train.to_csv(f'{output_folder}/{feat_name}_X_train.csv')
        X_test.to_csv(f'{output_folder}/{feat_name}_X_test.csv')
        y_train.to_csv(f'{output_folder}/{feat_name}_y_train.csv')
        y_test.to_csv(f'{output_folder}/{feat_name}_y_test.csv')
        # Save the model
        pickle_model(curr_model,
                    f'{output_folder}/{feat_name}_model.pickle')

    return auc_score


def merge_with_metadata(
    feature_table: pd.DataFrame,
    metadata_file: Path,
    output_folder: Path,
    table_name: str,
    logger: structlog.BoundLogger
) -> pd.DataFrame:
    """
    Add metadata from a file to the feature table and save it.
    """
    # Make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)

    # Define output file name
    out_csv_transposed = os.path.join(output_folder, 
                                      f'{table_name}_transposed.csv')
    out_csv = os.path.join(output_folder, f'{table_name}.csv')

    logger.info('Opening metadata file.', metadata_file=metadata_file)

    # Open metadata file
    metadata_df = pd.read_csv(metadata_file, sep=',',
                              index_col=0)
    
    metadata_df = metadata_df.T

    # Modify gender to 1-0
    metadata_df.loc['Gender'] = \
        metadata_df.loc['Gender'].replace({'male': 0, 'female': 1})

    # Transpose tables
    feature_table = feature_table.T
    metadata_df = metadata_df.T

    logger.info('Attempting merge with feature table.',
                metadata_file=metadata_file)

    # Ensure both dataframes have the same index and are aligned
    if not feature_table.index.equals(metadata_df.index):
        logger.warning("DataFrames indices do not match. Merging...",
                       metadata_file=metadata_file)
        # This is a robust way to merge if indices don't match
        merged_df = pd.merge(metadata_df, feature_table,
                             left_index=True, right_index=True, 
                             how='inner')
    else:
        merged_df = pd.concat([metadata_df, feature_table], axis=1)

    logger.info('Merge completed. Saving files...',
                metadata_file=metadata_file)
    
    # Save df as csv
    merged_df.to_csv(out_csv, sep=',')
    merged_df.T.to_csv(out_csv_transposed, sep=',')
    return merged_df


def normalize_by_gv(
    df_in: pd.DataFrame,
    logger: structlog.BoundLogger
) -> pd.DataFrame:
    """
    Normalizes a pandas DataFrame by number of GVs per sample.
    """
    logger.info('Separating the different kinds of feature.')

    # Get df for the different kinds of feature
    bins_df = df_in.loc[
        df_in.index.to_series().str.startswith('bin_gvs')
    ]
    genes_df = df_in.loc[
        df_in.index.to_series().str.startswith('gene_gvs')
    ]
    dn_ds_df = df_in.loc[
        df_in.index.to_series().str.startswith('dn_ds')
    ]
    fl_df = df_in.loc[
        df_in.index.to_series().str.startswith('fragment')
    ]
    cnv_df = df_in.loc[
        df_in.index.to_series().str.startswith('cnv_length')
    ]

    logger.info('Dividing by number of GVs.')

    # Get the sum of gvs per column (from bins_df)
    gv_col_sum = (bins_df.sum())
    # Divide by total GVs
    bins_df = bins_df / gv_col_sum
    genes_df = genes_df / gv_col_sum

    logger.info('Re-joining separated dfs.')

    # Join dataframes into df
    df_out = pd.concat([dn_ds_df, fl_df, cnv_df, genes_df, bins_df], 
                        ignore_index=False)
    # Return df_out
    return df_out


def normalize_feature_table(
    input_file: Path,
    output_folder: Path,
    logger: structlog.BoundLogger,
    log_transform: bool = False,
    robust_norm: bool = False,
    rand_seed: int|None = None
) -> List[pd.DataFrame]:
    """
    Normalize a feature table and create heatmaps.
    
    Args:
        input_file (Path): Path to the feature table CSV file.
        output_folder (Path): Output directory for heatmaps and 
                              normalized tables.
        logger: Logger instance.
        log_transform (bool): Defines if logarithmic transform is
                              performed on the data.
        robust_norm (bool): Defines if robust normalization is used 
                            instead of subtracting the minimum value 
                            and dividing by the average.
    """
    # Make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)

    # Define table name from input_file
    table_name = input_file.name.rsplit('.')[0]

    logger.info('Reading feature table.', table_path=input_file)

    # Read feature table
    df = pd.read_csv(input_file, sep=',')
    df.set_index(df.columns[0], inplace=True)

    # Define output file names
    out_raw = os.path.join(output_folder, f'{table_name}_RAW.csv')
    out_csv = os.path.join(output_folder, 
                           f'{table_name}_NORMALIZED.csv')
    out_max_per_row = os.path.join(output_folder,
                                   f'{table_name}_max_col.csv')
    out_min_per_row = os.path.join(output_folder,
                                   f'{table_name}_min_per_row.csv')

    df = cleanup_empty_rows_cols(df)

    # Copy raw df before numerical modifications
    raw_df = df.copy()

    df = normalize_by_gv(df, logger)

    if log_transform:
        logger.info('Perform log transformation.', 
                    table_path=input_file)
        # Perform log transformation
        df = df.apply(np.log1p)

    # Perform normalization
    if robust_norm:
        logger.info('Performing robust normalization.',
                    table_path=input_file)
        df = robust_normalize(df)
    else:
        logger.info('Performing normalization.',
                    table_path=input_file)
        df_max_per_row = df.max(axis=1)
        df = df.div(df_max_per_row, axis=0)
        # Save df_max_per_row
        df_max_per_row.to_csv(out_max_per_row,
                              header=['max_per_row'])

    df = cleanup_empty_rows_cols(df)

    logger.info('Removing minimums to avoid negative values.',
                table_path=input_file)

    # Remove minimums per row to avoid negative values
    df_min_per_row = df.min(axis=1)
    # Save df_min_per_row
    df_min_per_row.to_csv(out_min_per_row,
                          header=['min_per_row'])
    df = df.sub(df_min_per_row, axis=0)

    df = cleanup_empty_rows_cols(df)

    logger.info('Creating heatmaps.', table_path=input_file)
    # Get sample of the dataframe
    df_sample = sample_df(df, rand_seed, n=10)
    # Change feature names for better visualization
    df_sample.index = df_sample.index.map(process_feature_names)
    # Create the sample heatmap
    create_heatmap(df_sample, table_name+'_sample', True, 
                   output_folder, logger=logger)

    logger.info('Modifying feature names and saving dataframes.',
                table_path=input_file)

    # Modify feature names
    df.index = df.index.map(process_feature_names)
    # Save df as csv
    df.to_csv(out_csv, sep=',')

    # Save raw_df too
    raw_df.index = raw_df.index.map(process_feature_names)
    raw_df.to_csv(out_raw, sep=',')
    
    # Return raw df and df
    return raw_df, df


def per_feature_analysis(
    table_name: str,
    output_folder: Path,
    target_var: str,
    logger: structlog.BoundLogger,
    top_n: int=20,
    rand_seed: int=None
) -> None:
    """Performs per-feature analysis on raw or normalized data."""
    # Make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)
    # Define base name for outputs
    bname = str(table_name).rsplit('/')[-1].rsplit('.')[0]

    # Read feature table
    df = pd.read_csv(table_name, index_col=0)
    # Transpose df to have samples as rows and features as columns
    df = df.T
    
    model = RandomForestClassifier(random_state=rand_seed)

    logger.info("Analyzing feature table.", table_name=table_name)
    
    # Process features
    auc_results = analyze_features_and_rank_by_auc(
        df, output_folder, target_var, model, logger, rand_seed
    )

    logger.info("Feature analysis completed. Saving results.",
                table_name=table_name,
                feature_count=len(auc_results))
    
    # Define output file name
    file_name = f"{bname}_feature_singles.csv"
    file_out = output_folder / file_name

    # Save output
    auc_results.to_csv(file_out)

    if top_n >= 2:
        # Define top features from auc_results
        top_features = auc_results.head(top_n).index.tolist()

        logger.debug(f'Top {top_n} features obtained.',
                     top_features=top_features)

        logger.info(f'Analyzing top {top_n} features in pairs.',
                    table_name=table_name,
                    top_n=top_n)

        # Define pairs output name
        output_folder_pairs = output_folder / "pairs_auc_out"
        output_pairs = output_folder_pairs / \
            f"{bname}_feature_pairs.csv"

        # Run the analysis
        top_feature_pairs_ranked = analyze_feature_pairs(
            df,
            output_folder_pairs,
            top_features,
            target_var,
            logger,
            rand_seed
        )

        # Save the resulting dataframe
        top_feature_pairs_ranked.to_csv(output_pairs)
        
        if top_n >= 3:
            logger.info(f'Analyzing top {top_n} features in trios.',
                    table_name=table_name, top_n=top_n)
            # Define trios output name
            output_folder_trios = output_folder / "trios_auc_out"
            output_trios = output_folder_trios / \
                f"{bname}_feature_trios.csv"

            # Run the analysis for trios
            top_feature_trios_ranked = analyze_feature_trios(
                df,
                output_folder_trios,
                top_features,
                target_var,
                logger,
                rand_seed
            )

            # Save the resulting dataframe
            top_feature_trios_ranked.to_csv(output_trios)

    return None


def pickle_load(model_path: Path):
    """Loads a pickled model from the specified input path."""
    with open(model_path, "rb") as f:
        model_file = pickle.load(f)
    return model_file


def pickle_model(model, output_path: Path) -> None:
    """Pickles a model to the specified output path."""
    with open(output_path, 'wb') as f:
        pickle.dump(model, f)
    return None


def plot_feature_roc(
    df: pd.DataFrame, 
    target_var: str, 
    output_folder: Path,
    feature_name: str,
    fitted_model: Any,
    logger: structlog.BoundLogger,
    random_state: int|None = None
) -> None:
    """
    Generates and displays a multi-class One-vs-Rest ROC curve plot 
    for a single feature using the provided model, evaluated on a 
    held-out test set.
    
    Args:
        df (pd.DataFrame): The full DataFrame.
        target_var (str): Name of the target variable column.
        output_folder (Path): Directory to store the final models.
        feature_name (str): The name of the single feature to plot.
        fitted_model: The trained RF model for this feature 
                      (fitted on full data).
        logger: Logger instance.
        random_state (int|None): Seed for reproducibility.
    """

    logger.info(f"Generating ROC curve for feature: {feature_name} ")
    
    # Make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)

    # Define X and y
    X = df.drop(columns=[target_var])
    y = df[target_var]

    # Use LabelBinarizer for easy One-vs-Rest handling
    label_binarizer = LabelBinarizer()
    label_binarizer.fit(y)
    classes = label_binarizer.classes_
    n_classes = len(classes)
    
    # Prepare Data and Split
    X_feature = X[[feature_name]]
    
    # Split the data to get a test set for plotting AUC
    X_train, X_test, y_train_orig, y_test_orig = train_test_split(
        X_feature, y, 
        test_size=0.33, 
        stratify=y, 
        random_state=random_state
    )
    
    # Binarize the test labels
    y_test_bin = label_binarizer.transform(y_test_orig)

    # Predict probabilities using the provided fitted model
    y_proba = fitted_model.predict_proba(X_test)

    # Calculate OVR AUC
    roc_auc_ovr_macro = roc_auc_score(
        y_test_orig,
        y_proba,
        multi_class="ovr",
        average="macro"
    )

    # Plot the figure
    plt.figure(figsize=(8, 6))
    
    fpr = dict()
    tpr = dict()
    roc_auc_class = dict()

    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test_bin[:, i], y_proba[:, i])
        roc_auc_class[i] = auc(fpr[i], tpr[i])
        str_label = f'ROC curve for {classes[i]} '+\
            f'(AUC = {roc_auc_class[i]:.2f})'
        plt.plot(
            fpr[i], 
            tpr[i], 
            label=str_label
        )

    plt.plot([0, 1], [0, 1], 'k--', label='Chance Level (AUC = 0.50)')
    
    plot_title = 'ROC Curve (Test Split) for Feature: ' +\
        f'{feature_name}\n(Overall Macro OVR AUC: ' +\
        f'{roc_auc_ovr_macro:.3f})'
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(plot_title)
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    feat_name = feature_name.replace(' ', '_')
    plt.savefig(output_folder / f'roc_curve_{feat_name}.png')

    logger.info(f"Generated ROC curve for feature: {feature_name} " +\
                f"(Macro AUC: {roc_auc_ovr_macro:.3f})")
    
    return None


def process_feature_names(row_name:str) -> str:
    """Processes a feature row name and formats it."""
    ret = row_name
    if row_name.startswith('bin_gvs'):
        region = row_name.rsplit('_')[-1]
        ret = f'GVs in region {region}'
    elif row_name.startswith('gene_gvs'):
        gene_name = row_name.rsplit('_')[-1]
        ret = f'GVs in gene {gene_name}'
    elif row_name.startswith('fragment'):
        feat = row_name.rsplit('_')[-1]
        if feat == 'mean':
            ret = 'Mean fragment length'
        elif feat == 'median':
            ret = 'Median fragment length'
        elif feat == 'std':
            ret = 'Fragment length standard deviation'
        elif feat == 'min':
            ret = 'Minimum fragment length'
        elif feat == 'max':
            ret = 'Maximum fragment length'
        elif feat == 'count':
            ret = 'Fragment count'
        else:
            er = f'WARNING: row_name {row_name} not processed properly.'
            print(er)
    elif row_name.startswith('cnv_length'):
        chr_n = row_name.rsplit('_')[-1]
        gain_loss = row_name.rsplit('_')[-2]
        if gain_loss == 'gain':
            ret = f'Nucleotides added by CNVs in {chr_n}'
        elif gain_loss == 'loss':
            ret = f'Nucleotides lost by CNVs in {chr_n}'
        else:
            er = f'WARNING: row_name {row_name} not processed properly.'
            print(er)
    elif row_name.startswith('dn_ds'):
        gene_name = row_name.rsplit('_')[-1]
        ret = f'dN/dS proportion in gene {gene_name}'
    else:
        print(f'row_name {row_name} not processed.')
    return ret


def robust_normalize(df:pd.DataFrame) -> pd.DataFrame:
    """Performs robust normalization on a dataframe."""
    # Define Q1 and Q3
    q1 = df.quantile(0.25, axis=1)
    q3 = df.quantile(0.75, axis=1)
    # Define rows to keep
    rows_to_keep = (q3 != 0) | (q1 != 0)

    # Filter the DataFrame
    df_filtered = df[rows_to_keep]

    print('# df_filtered:')
    print(df_filtered)

    # Perform robust scaling
    # Create the scaler object
    scaler = RobustScaler()

    # Transpose the dataframe
    df_transposed = df_filtered.T
    print('# df_transposed:')
    print(df_transposed)
    # Scale the transposed DataFrame
    scaled_data = scaler.fit_transform(df_transposed)

    print('# Scaled raw df:')
    print(scaled_data)

    # Transpose the result back and create a new DataFrame
    df_row_scaled = pd.DataFrame(scaled_data.T,
                                 columns=df_filtered.columns,
                                 index=df_filtered.index)
    return df_row_scaled


def run_model_validation_and_test(
    model_path: Path,
    data_table_path: Path,
    out_folder: Path,
    target_var: str,
    logger: structlog.BoundLogger,
    out_name: str=''
):
    """
    Runs a model on a data table and creates performance metrics.
    Includes the creation of a confusion matrix figure.
    """
    # Make sure output folder exists
    out_folder.mkdir(parents=True, exist_ok=True)

    logger.info('Loading model.', model_path=model_path)

    # Load model
    model = pickle_load(model_path)
    # Define basename from model_path
    bname = model_path.name.rsplit('.')[0]

    logger.info('Reading data table.', data_table=data_table_path)

    # Load data table
    df = pd.read_csv(data_table_path, sep=',', index_col=0)
    # Define X and y values
    X = df.T.drop(columns=[target_var])
    y = df.T[target_var]
    # Binarize y values
    ALL_POSSIBLE_CLASSES = ['Healthy', 'CRC', 'BRC'] 
    lb = LabelBinarizer()
    lb.fit(ALL_POSSIBLE_CLASSES)
    y_true_bin = lb.transform(y)
    # Select only the features used in the model
    model_features = model.feature_names_in_
    X = X[model_features]

    logger.info('Evaluating model on data table.',
                data_table=data_table_path,
                feature_count=len(model_features),
                sample_count=len(X))

    # Evaluate the model on the test set
    y_pred = model.predict_proba(X)

    # Define output confusion matrix path
    out_conf = out_folder / f'conf_matrix_{bname}{out_name}.png'
    # Transform one-hot encoded lists to labels
    y_true_labels = np.argmax(y_true_bin, axis=1)
    y_pred_labels = np.argmax(y_pred, axis=1)

    logger.info('Creating confusion matrix.',
                data_table=data_table_path,
                model_path=model_path)
    
    # Create confusion matrix figures
    create_conf_matrix(y_true_labels, y_pred_labels, 
                       ALL_POSSIBLE_CLASSES, out_conf)
    
    logger.info('Calculating performance metrics.',
                data_table=data_table_path,
                model_path=model_path)

    # Calculate AUC score
    auc_score = roc_auc_score(y_true_bin, y_pred,
                              multi_class='ovr')
    
    # Create ROC curve figure
    create_auc_roc_curves(
        y_true_bin,
        y_pred,
        out_folder,
        bname,
        out_name,
        model,
        logger,
    )

    create_micro_avg_auc_roc_curve(
        y_true_bin,
        y_pred,
        out_folder,
        bname,
        out_name,
        logger
    )

    create_macro_avg_auc_roc_curve(
        y_true_bin,
        y_pred,
        out_folder,
        bname,
        out_name,
        logger
    )

    # Obtain accuracy, precision, recall, f1
    accuracy = accuracy_score(y_true_labels, y_pred_labels)
    precision = precision_score(y_true_labels, y_pred_labels,
                                average='weighted', zero_division=0)
    recall = recall_score(y_true_labels, y_pred_labels,
                          average='weighted', zero_division=0)
    f1 = f1_score(y_true_labels, y_pred_labels,
                  average='weighted', zero_division=0)
    
    logger.info('Performance metrics calculated. Saving results.',
                data_table=data_table_path,
                model_path=model_path)
    
    # Save results to a text file
    out_txt = out_folder / f'performance_{bname}{out_name}.txt'
    with open(out_txt, 'w') as f:
        f.write(f'AUC: {auc_score}\n')
        f.write(f'Accuracy: {accuracy}\n')
        f.write(f'Precision: {precision}\n')
        f.write(f'Recall: {recall}\n')
        f.write(f'F1: {f1}\n')
    # Return the performance metrics into a dictionary
    results = {
        'AUC': auc_score,
        'Accuracy': accuracy,
        'Precision': precision,
        'Recall': recall,
        'F1': f1
    }
    return results


def sample_df(
    df: pd.DataFrame,
    rand_seed: int|None,
    n: int=10
) -> pd.DataFrame:
    """
    Takes a sample of n rows from each feature type.
    """
    # Get df for the different kinds of feature
    bins_df = df.loc[df.index.to_series().str.startswith('bin_gvs')]
    genes_df = df.loc[df.index.to_series().str.startswith('gene_gvs')]
    dn_ds_df = df.loc[df.index.to_series().str.startswith('dn_ds')]
    fl_df = df.loc[df.index.to_series().str.startswith('fragment')]
    cnv_df = df.loc[df.index.to_series().str.startswith('cnv_length')]
    # Create a heatmap with samples from subsets of features
    l_df = [
        fl_df,
        cnv_df.sample(n=n, random_state=rand_seed),
        genes_df.sample(n=n, random_state=rand_seed),
        bins_df.sample(n=n, random_state=rand_seed)
    ]
    if not dn_ds_df.empty:
        l_df.append(dn_ds_df.sample(n=n, random_state=rand_seed))
    df_sample = pd.concat(l_df, ignore_index=False)

    return df_sample


def top_single_features_cvrf(
    df: pd.DataFrame,
    target_var: str,
    output_path: Path,
    logger: structlog.BoundLogger,
    split_n: int = 5,
    random_state: int|None = None
) -> pd.Series:
    """
    Performs cross-validated Random Forest (RF) classification on each
    feature in X individually and returns the mean AUC for each feature.

    This function handles multi-class classification by using 'ovr'
    (One-vs-Rest) approach for AUC calculation.

    Returns a pandas Series mapping feature names to their mean 
    cross-validated AUC, sorted in descending order of performance.
    """
    # Make sure output folder exists
    output_path.mkdir(parents=True, exist_ok=True)
    # Define X and y
    X = df.drop(columns=[target_var])
    y = df[target_var]

    logger.info(
        f"Starting single-feature CVRF on {X.shape[1]} features..."
    )
    
    # Encode target labels
    label_encoder = LabelEncoder()
    y_encoded = label_encoder.fit_transform(y)
    n_classes = len(label_encoder.classes_)
    
    if n_classes < 2:
        logger.error(
            "Target variable must have at least two unique classes."
        )
        return pd.Series({})

    # Initialize CV and Model
    cv = StratifiedKFold(n_splits=split_n, shuffle=True,
                         random_state=random_state)
    
    # Use conservative RF parameters for the single-feature model
    rf_base_model = RandomForestClassifier(
        n_estimators=100,
        max_depth=5,
        class_weight='balanced',
        random_state=random_state,
        n_jobs=-1
    )

    results = {}

    ### Display
    # Counter
    cont = 0
    n_col = len(X.columns)
    ###

    # Iterate through each feature
    for feature_name in X.columns:
        # Extract single feature data and reshape
        X_feature = X[[feature_name]].values 
        fold_auc_scores = []
        
        # Clone the model for each feature run
        rf_model = rf_base_model

        # Perform Cross-Validation
        for fold, (train_index, test_index) in \
            enumerate(cv.split(X_feature, y_encoded)):
            X_train, X_test = \
                X_feature[train_index], X_feature[test_index]
            y_train, y_test = \
                y_encoded[train_index], y_encoded[test_index]

            try:
                # Train the model
                rf_model.fit(X_train, y_train)
                
                # Predict probabilities
                y_proba = rf_model.predict_proba(X_test)
                
                # Calculate AUC
                auc = roc_auc_score(
                    y_test, 
                    y_proba, 
                    multi_class='ovr', 
                    average='macro',
                    labels=np.arange(n_classes)
                )
                fold_auc_scores.append(auc)
                
                # Pickle the model
                feat_name = feature_name.replace(' ', '_')
                model_path = output_path / \
                    f'cvrf_model_{feat_name}_fold{fold+1}.pickle'
                pickle_model(rf_model, Path(model_path))

            except ValueError as e:
                logger.warning(f"Skipping fold {fold} for feature " +\
                               f"{feature_name} due to error: {e}")
        ### Display
        if cont==0 or (cont+1)%1000==0:
            logger.debug(f'Progress: {cont+1} / {n_col}')
        cont += 1
        ###
        if fold_auc_scores:
            mean_auc = np.mean(fold_auc_scores)
            results[feature_name] = mean_auc
        else:
            results[feature_name] = np.nan
    
    # Return sorted results
    results_series = \
        pd.Series(results).sort_values(ascending=False).dropna()
    
    logger.info("CVRF classification complete.")

    return results_series


def train_and_store_final_models(
    df: pd.DataFrame,
    target_var: str,
    output_folder: Path,
    top_feature_names: List[str],
    auc_results: pd.Series,
    logger: structlog.BoundLogger,
    random_state: int|None = None
) -> Dict:
    """
    Trains a final Random Forest model on the entire dataset for each
    of the top performing features and stores the model object.

    Args:
        df (pd.DataFrame): The full DataFrame.
        target_var (str): Name of the target variable column.
        output_folder (Path): Directory to store the final models.
        top_feature_names (List[str]): List of feature names to train
                                       final models for.
        auc_results (pd.Series): Series of CVRF AUC scores for 
                                 reporting.
        logger: Logger instance.
        random_state (int|None): Seed for reproducibility.

    Returns:
        A dictionary mapping feature name to a tuple: 
        (mean_cv_auc, fitted_model).
    """
    # Make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)

    # Define X and y
    X = df.drop(columns=[target_var])
    y = df[target_var]

    final_models = {}
    
    # Prepare base model and target encoding
    rf_model = RandomForestClassifier(
        n_estimators=100, 
        max_depth=5, 
        class_weight='balanced', 
        random_state=random_state,
        n_jobs=-1
    )
    y_encoded = LabelEncoder().fit_transform(y)
    
    logger.info("Training and storing models for the " +\
                f"top {len(top_feature_names)} features...")
    
    # Train a model on the full data for each top feature
    for feature_name in top_feature_names:
        X_feature = X[[feature_name]].values
        # Copy rf model per feature
        rf_final_model = rf_model
        # Train the model on ALL data for the final saved model
        rf_final_model.fit(X_feature, y_encoded)
        
        # Store the CV AUC score and the final fitted model object
        mean_cv_auc = auc_results.loc[feature_name]
        final_models[feature_name] = (mean_cv_auc, rf_final_model)
        
        # Pickle the feature model
        feat_name = feature_name.replace(' ', '_')
        model_path = output_folder / f'{feat_name}_model.pickle'
        pickle_model(rf_final_model, model_path)
        
    return final_models
