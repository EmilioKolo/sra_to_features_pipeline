"""
ML-ready feature table utilities for the SRA to Features Pipeline.
"""


from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import RobustScaler, LabelBinarizer
from typing import List, Dict, Any
import itertools
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
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

    def save_feature_table(self, output_path: Path, format: str = 'csv'):
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
    
    def _extract_genomic_features(self, genomic_bins: List) -> Dict[str, Any]:
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
            'bins_with_variants': sum(1 for count in variant_counts if count > 0),
            'total_genomic_bins': len(genomic_bins),
            'variant_density': sum(variant_counts) / len(genomic_bins),
        }
    
    def _extract_gene_features(self, gene_stats: List) -> Dict[str, Any]:
        """Extract gene-level variant features."""
        if not gene_stats:
            return {}
        
        # Calculate summary statistics
        total_variants = [gene.total_variants for gene in gene_stats]
        synonymous_variants = [gene.synonymous_variants for gene in gene_stats]
        nonsynonymous_variants = [gene.nonsynonymous_variants for gene in gene_stats]
        dn_ds_ratios = [gene.dn_ds_ratio for gene in gene_stats if gene.dn_ds_ratio is not None]
        
        return {
            'total_genes_with_variants': len(gene_stats),
            'total_gene_variants': sum(total_variants),
            'mean_variants_per_gene': np.mean(total_variants),
            'std_variants_per_gene': np.std(total_variants),
            'total_synonymous_variants': sum(synonymous_variants),
            'total_nonsynonymous_variants': sum(nonsynonymous_variants),
            'mean_dn_ds_ratio': np.mean(dn_ds_ratios) if dn_ds_ratios else 0,
            'std_dn_ds_ratio': np.std(dn_ds_ratios) if dn_ds_ratios else 0,
            'genes_with_high_variants': sum(1 for count in total_variants if count > 10),
        }
    
    def _extract_cnv_features(self, cnv_regions: List) -> Dict[str, Any]:
        """Extract CNV features."""
        if not cnv_regions:
            return {}
        
        # Calculate summary statistics
        copy_numbers = [region.copy_number for region in cnv_regions]
        confidences = [region.confidence for region in cnv_regions]
        gains = sum(1 for region in cnv_regions if region.type == 'gain')
        losses = sum(1 for region in cnv_regions if region.type == 'loss')
        
        return {
            'total_cnv_regions': len(cnv_regions),
            'cnv_gains': gains,
            'cnv_losses': losses,
            'mean_copy_number': np.mean(copy_numbers),
            'std_copy_number': np.std(copy_numbers),
            'mean_cnv_confidence': np.mean(confidences),
            'high_confidence_cnvs': sum(1 for conf in confidences if conf > 0.9),
            'cnv_burden': len(cnv_regions),  # Total CNV burden
        }
    
    def save_feature_table(self, output_path: Path, format: str = 'csv'):
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


def analyze_feature_pairs(
    df: pd.DataFrame,
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
        features_to_analyze (list): A list of feature names to be 
                                    paired.
        target_column (str): The name of the target variable column.

    Returns:
        pd.DataFrame: A DataFrame of feature pairs ranked by AUC.
    """
    if len(features_to_analyze) < 2:
        logger.error(
            "Please provide at least two features for pairing.",
            feature_n=len(features_to_analyze)
        )
        return pd.DataFrame()

    X = df.T[features_to_analyze]
    y = df.T[target_column]

    # Generate all unique pairs of the selected features
    feature_pairs = list(itertools.combinations(features_to_analyze, 2))
    logger.info(f"Generated {len(feature_pairs)} unique feature pairs from the list of {len(features_to_analyze)} features.",
                feature_n=len(features_to_analyze))

    pair_results = {}
    model = RandomForestClassifier(random_state=rand_seed)
    lb = LabelBinarizer()

    ### Display
    # Counter
    cont = 0
    n_feat_pairs = len(feature_pairs)
    ###

    for f1, f2 in feature_pairs:
        # Create a new DataFrame with only the two features for training
        X_pair = X[[f1, f2]]

        X_train, X_test, y_train, y_test = train_test_split(
            X_pair, y, test_size=0.3, random_state=rand_seed, stratify=y
        )

        # Handle multi-class classification for AUC calculation
        y_test_bin = lb.fit_transform(y_test)
        y_pred_proba = model.fit(X_train, y_train).predict_proba(X_test)

        if y_test_bin.shape[1] == 1:
            auc_score = roc_auc_score(y_test_bin, y_pred_proba[:, 1])
        else:
            auc_score = roc_auc_score(y_test_bin, y_pred_proba, multi_class='ovr')

        pair_results[(f1, f2)] = auc_score

        ### Display
        if cont==0 or (cont+1)%10==0:
            logger.debug(f'Progress: {cont+1} / {n_feat_pairs}',
                         feature_n=len(features_to_analyze))
        cont += 1
        ###
    
    # Create a DataFrame and sort the results
    pair_auc_df = pd.DataFrame.from_dict(
        pair_results, orient='index', columns=['AUC Score']
    )
    pair_auc_df.index.name = 'Feature Pair'
    pair_auc_df = pair_auc_df.sort_values(by='AUC Score', ascending=False)

    return pair_auc_df


def analyze_feature_trios(
    df: pd.DataFrame,
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
        features_to_analyze (list): A list of feature names to be 
                                    grouped.
        target_column (str): The name of the target variable column.

    Returns:
        pd.DataFrame: A DataFrame of feature pairs ranked by AUC.
    """
    if len(features_to_analyze) < 3:
        logger.error(
            "Please provide at least three features for pairing.",
            feature_n=len(features_to_analyze)
        )
        return pd.DataFrame()

    X = df.T[features_to_analyze]
    y = df.T[target_column]

    # Generate all unique trios of the selected features
    feature_trios = list(itertools.combinations(features_to_analyze, 3))
    logger.info(f"Generated {len(feature_trios)} unique feature trios from the list of {len(features_to_analyze)} features.",
                feature_n=len(features_to_analyze))

    pair_results = {}
    model = RandomForestClassifier(random_state=rand_seed)
    lb = LabelBinarizer()

    ### Display
    # Counter
    cont = 0
    n_feat_trios = len(feature_trios)
    ###

    for f1, f2, f3 in feature_trios:
        # Create a new DataFrame with only the features for training
        X_pair = X[[f1, f2, f3]]

        X_train, X_test, y_train, y_test = train_test_split(
            X_pair, y, test_size=0.3, random_state=rand_seed, stratify=y
        )

        # Handle multi-class classification for AUC calculation
        y_test_bin = lb.fit_transform(y_test)
        y_pred_proba = model.fit(X_train, y_train).predict_proba(X_test)

        if y_test_bin.shape[1] == 1:
            auc_score = roc_auc_score(y_test_bin, y_pred_proba[:, 1])
        else:
            auc_score = roc_auc_score(y_test_bin, y_pred_proba, 
                                      multi_class='ovr')

        pair_results[(f1, f2, f3)] = auc_score

        ### Display
        if cont==0 or (cont+1)%100==0:
            logger.debug(f'Progress: {cont+1} / {n_feat_trios}',
                         feature_n=len(features_to_analyze))
        cont += 1
        ###
    
    # Create a DataFrame and sort the results
    pair_auc_df = pd.DataFrame.from_dict(
        pair_results, orient='index', columns=['AUC Score']
    )
    pair_auc_df.index.name = 'Feature Pair'
    pair_auc_df = pair_auc_df.sort_values(by='AUC Score', 
                                          ascending=False)

    return pair_auc_df


def analyze_features_and_rank_by_auc(
    df: pd.DataFrame,
    target_var: str,
    logger: structlog.BoundLogger,
    rand_seed:int=None
) -> pd.DataFrame:
    """
    Performs feature-by-feature analysis and ranks them by AUC score.

    Args:
        df (pd.DataFrame): The DataFrame with features and a 
                           target variable column.
        target_var (str): Target variable to be predicted with the 
                          feature classification.
        output_folder (Path): Folder where output files are created.
        logger: Logger instance.
        rand_seed (int): Random seed for repeatability.

    Returns:
        pd.DataFrame: A DataFrame of features ranked by AUC.
    """
    df = df.T
    X = df.drop(target_var, axis=1)
    y = df[target_var]

    results = {}
    model = RandomForestClassifier(random_state=rand_seed)

    ### Display
    # Counter
    cont = 0
    n_col = len(X.columns)
    ###
    for feature in X.columns:
        auc_score = get_single_feature_auc(X, y, feature, model,
                                           rand_seed)
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
    out_path:str
):
    """Creates a heatmap using Seaborn and Matplotlib"""
    # Check that df is not empty
    if df.empty:
        print(f'# WARNING: Heatmap with title "{title}" could not be created.')
        return None
    print(f'# Creating heatmap with title "{title}"...')
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
    plt.show()
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


def get_single_feature_auc(
    X:pd.DataFrame,
    y:pd.Series,
    feature_name:str,
    model,
    rand_seed:int
) -> float:
    """
    Trains a model using a single feature and returns the ROC-AUC score.

    Args:
        X (pd.DataFrame): DataFrame of features.
        y (pd.Series): Series of target labels.
        feature_name (str): The name of the feature to use.
        model: The machine learning model to train.
        rand_seed (int): Random seed for repeatability.

    Returns:
        float: The ROC-AUC score.
    """
    X_single = X[[feature_name]]
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
    y_pred_proba = model.fit(X_train, y_train).predict_proba(X_test)
    
    # Calculate AUC score
    if y_test_bin.shape[1] == 1:
        auc_score = roc_auc_score(y_test_bin, y_pred_proba[:, 1])
    else:
        auc_score = roc_auc_score(y_test_bin, y_pred_proba, 
                                  multi_class='ovr')

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
        merged_df = pd.merge(feature_table, metadata_df, 
                             left_index=True, right_index=True, 
                             how='inner')
    else:
        merged_df = pd.concat([feature_table, metadata_df], axis=1)
    
    logger.info('Merge completed. Saving files...',
                metadata_file=metadata_file)
    
    # Save df as csv
    merged_df.to_csv(out_csv, sep=',')
    merged_df.T.to_csv(out_csv_transposed, sep=',')
    return merged_df


def normalize_feature_table(
    input_file: Path,
    output_folder: Path,
    logger: structlog.BoundLogger,
    log_transform: bool = False,
    robust_norm: bool = False
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

    # Remove rows that sum 0
    df = df[df.sum(axis=1)!=0]
    # Remove columns that sum 0
    df = df.loc[:, (df.sum(axis=0)!=0)]

    # Copy raw df before numerical modifications
    raw_df = df.copy()

    logger.info('Separating the different kinds of feature.', 
                table_path=input_file)

    # Get df for the different kinds of feature
    bins_df = df.loc[df.index.to_series().str.startswith('bin_gvs')]
    genes_df = df.loc[df.index.to_series().str.startswith('gene_gvs')]
    dn_ds_df = df.loc[df.index.to_series().str.startswith('dn_ds')]
    fl_df = df.loc[df.index.to_series().str.startswith('fragment')]
    cnv_df = df.loc[df.index.to_series().str.startswith('cnv_length')]

    logger.info('Dividing by number of GVs.', table_path=input_file)

    # Get the sum of gvs per column (from bins_df)
    gv_col_sum = (bins_df.sum())
    # Divide by total GVs
    bins_df = bins_df / gv_col_sum
    genes_df = genes_df / gv_col_sum

    if log_transform:
        logger.info('Perform log transformation.', 
                    table_path=input_file)
        # Perform log transformation
        dn_ds_df = dn_ds_df.apply(np.log1p)
        fl_df = fl_df.apply(np.log1p)
        cnv_df = cnv_df.apply(np.log1p)
        genes_df = genes_df.apply(np.log1p)
        bins_df = bins_df.apply(np.log1p)

    logger.info('Re-join separated dfs.', table_path=input_file)

    # Join dataframes into df
    df = pd.concat([dn_ds_df, fl_df, cnv_df, genes_df, bins_df], 
                    ignore_index=False)

    # Perform normalization
    if robust_norm:
        logger.info('Performing robust normalization.',
                table_path=input_file)
        df = robust_normalize(df)
    else:
        logger.info('Performing normalization.',
                table_path=input_file)
        df = df.div(df.max(axis=1), axis=0)

    logger.info('Re-separate the different kinds of feature.', 
                table_path=input_file)

    # Get df for the different kinds of feature again
    bins_df = df.loc[df.index.to_series().str.startswith('bin_gvs')]
    genes_df = df.loc[df.index.to_series().str.startswith('gene_gvs')]
    dn_ds_df = df.loc[df.index.to_series().str.startswith('dn_ds')]
    fl_df = df.loc[df.index.to_series().str.startswith('fragment')]
    cnv_df = df.loc[df.index.to_series().str.startswith('cnv_length')]

    logger.info('Removing minimums to avoid negative values.', 
                table_path=input_file)

    # Remove minimums to avoid negative values
    cnv_df = cnv_df - min(cnv_df.min().min(), 0)
    fl_df = fl_df - min(fl_df.min().min(), 0)
    bins_df = bins_df - min(bins_df.min().min(), 0)
    genes_df = genes_df - min(genes_df.min().min(), 0)
    if not dn_ds_df.empty:
        dn_ds_df = dn_ds_df - min(dn_ds_df.min().min(), 0)

    logger.info('Re-joining separated dfs.', table_path=input_file)

    # Join dataframes into df
    df = pd.concat([dn_ds_df, fl_df, cnv_df, genes_df, bins_df], 
                    ignore_index=False)

    logger.info('Creating heatmaps.', table_path=input_file)

    # Define if heatmaps will be created
    create_heatmaps = False
    if create_heatmaps:
        save_tables = True
        # Create heatmaps
        create_heatmap(dn_ds_df, table_name+'_dn_ds', save_tables,
                    output_folder)
        create_heatmap(fl_df, table_name+'_fl', save_tables,
                    output_folder)
        create_heatmap(cnv_df, table_name+'_cnvs', save_tables,
                    output_folder)
        create_heatmap(genes_df, table_name+'_genes', save_tables,
                    output_folder)
        create_heatmap(bins_df, table_name+'_bins', save_tables,
                    output_folder)
        create_heatmap(df, table_name, save_tables, output_folder)
    
    # Get df for the different kinds of feature
    bins_df = df.loc[df.index.to_series().str.startswith('bin_gvs')]
    genes_df = df.loc[df.index.to_series().str.startswith('gene_gvs')]
    dn_ds_df = df.loc[df.index.to_series().str.startswith('dn_ds')]
    fl_df = df.loc[df.index.to_series().str.startswith('fragment')]
    cnv_df = df.loc[df.index.to_series().str.startswith('cnv_length')]
    # Create a heatmap with samples from subsets of features
    l_df = [
        fl_df,
        cnv_df.sample(n=10, random_state=12),
        genes_df.sample(n=10, random_state=12),
        bins_df.sample(n=10, random_state=12)
    ]
    if not dn_ds_df.empty:
        l_df.append(dn_ds_df.sample(n=10, random_state=12))
    df_sample = pd.concat(l_df, ignore_index=False)
    df_sample.index = df_sample.index.map(process_feature_names)
    create_heatmap(df_sample, table_name+'_sample', True, 
                   output_folder)

    print('### Pre-feature modification df:')
    print(df)

    logger.info('Modifying feature names and saving dataframes.', 
                table_path=input_file)

    # Modify feature names
    df.index = df.index.map(process_feature_names)

    print('### Post-feature modification df:')
    print(df)

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
):
    """Performs per-feature analysis on raw or normalized data."""
    # Make sure output folder exists
    output_folder.mkdir(parents=True, exist_ok=True)

    # Process features
    logger.info("Analyzing feature table.", table_name=table_name)
    df = pd.read_csv(table_name, index_col=0)
    auc_results = analyze_features_and_rank_by_auc(
        df, target_var, logger, rand_seed
    )

    logger.info("Feature analysis completed. Saving results.",
                table_name=table_name,
                feature_count=len(auc_results))
    # Define output file name
    bname = str(table_name).rsplit('/')[-1].rsplit('.')[0]
    file_name = f"{bname}_feature_singles.csv"
    file_out = output_folder / file_name

    # Save output
    auc_results.to_csv(file_out)

    if top_n >= 2:
        # Define top features from auc_results
        top_features = auc_results.head(top_n).index.tolist()

        logger.debug(f'Top {top_n} features:',
                     top_features=top_features)

        logger.info(f'Analyzing top {top_n} features in pairs.',
                    table_name=table_name,
                    top_n=top_n)

        # Define pairs output name
        output_pairs = output_folder / f"{bname}_feature_pairs.csv"

        # Run the analysis
        top_feature_pairs_ranked = analyze_feature_pairs(
            df, top_features, target_var, logger, rand_seed
        )

        # Save the resulting dataframe
        top_feature_pairs_ranked.to_csv(output_pairs)
        
        if top_n >= 3:
            logger.info(f'Analyzing top {top_n} features in trios.',
                    table_name=table_name,
                    top_n=top_n)
            # Define trios output name
            output_trios = output_folder / f"{bname}_feature_trios.csv"

            # Run the analysis for trios
            top_feature_trios_ranked = analyze_feature_trios(
                df, top_features, target_var, logger, rand_seed
            )

            # Save the resulting dataframe
            top_feature_trios_ranked.to_csv(output_trios)

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
        print(f'WARNING: row_name {row_name} not processed properly.')
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
