"""
ML-ready feature table utilities for the SRA to Features Pipeline.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Any, Optional, Union
import structlog
import json

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
            len_bin = abs(int(end) - int(start))
            # Define key and value
            key = f'bin_gvs_{chr_n}:{start}-{end}'
            val = curr_data['variant_count']
            feature_dict[key] = float(val)/len_bin if len_bin>0 else 0.0
        
        # Go through gene_stats
        for i in range(len(gene_stats)):
            curr_data = gene_stats[i]
            gene_name = curr_data['gene_name'].split(':')[-1]
            start = curr_data['start']
            end = curr_data['end']
            len_gene = abs(int(end) - int(start))
            total_gv = int(curr_data['total_variants'])
            dn_ds = curr_data['dn_ds_ratio']
            # Define key and value for GVs
            key_gv = f'gene_gvs_{gene_name}'
            feature_dict[key_gv] = float(total_gv) / len_gene if len_gene>0 else 0.0
            # Check if there are GVs to calculate dn/ds
            if total_gv!=0 and dn_ds is not None:
                key_dn_ds = f'gene_dn_ds_{gene_name}'
                feature_dict[key_dn_ds] = min(float(dn_ds) / 1000, 1)
        
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
        # Normalize cnv lengths
        for cnv_key in l_cnv_keys:
            feature_dict[cnv_key] = min(float(feature_dict[cnv_key]) / (100*1000*1000), 1)
        
        # Go through fragment_stats
        for key, val in fragment_stats.items():
            if key in ['max', 'min', 'mean', 'median', 'std']:
                feature_dict[f'fragment_length_{key}'] = float(val) / 1000
            else:
                feature_dict[f'fragment_{key}'] = min(float(val) / (100*1000*1000), 1)
        
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
        output_path: Output file path for feature table
        format: Output format ('csv', 'tsv', 'parquet', 'json')
        logger: Logger instance
        
    Returns:
        DataFrame with samples × features
    """
    ml_table = FeatureTable(logger)
    ml_table.add_feature_dicts_from_directory(input_directory)
    
    if not ml_table.feature_dicts:
        raise ValueError(f"No feature sets found in {input_directory}")
    
    # Save feature table
    ml_table.save_feature_table(output_path, format)
    
    # Log summary
    summary = ml_table.get_feature_summary()
    logger.info("ML feature table created successfully", **summary)
    
    return ml_table.create_sample_features_table()


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