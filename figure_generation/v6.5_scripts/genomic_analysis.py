#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GenomicAnalysis class for performing comparative analyses between genomic samples.
@author: michelemeline
"""

import numpy as np
import pandas as pd
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt
from matplotlib import patheffects as pe
from scipy.stats import linregress
from scipy import stats
import seaborn as sns


class GenomicAnalysis:
    """
    Class to perform comparative analyses between genomic samples.
    """
    
    def __init__(self, samples):
        """
        Initialize analysis with a dictionary of samples.
        
        Parameters:
        -----------
        samples : dict
            Dictionary mapping sample names to GenomicSample objects
        """
        self.samples = samples
    
    def create_venn2(self, sample1_name, sample2_name, feature_type='enhancers', 
                     output_file=None):
        """
        Create a 2-way Venn diagram comparing two samples.
        Uses colors assigned to each sample.
        
        Parameters:
        -----------
        sample1_name, sample2_name : str
            Names of samples to compare
        feature_type : str
            Feature type to compare ('enhancers', 'super_enhancers', etc.)
        output_file : str, optional
            Path to save the figure
        """
        s1 = self.samples[sample1_name]
        s2 = self.samples[sample2_name]
        
        bed1 = getattr(s1, feature_type)
        bed2 = getattr(s2, feature_type)
        
        overlap = bed1.intersect(bed2, f=0.5, r=True, wo=True)
        
        only_s1 = bed1.count() - overlap.count()
        only_s2 = bed2.count() - overlap.count()
        both = overlap.count()
        
        plt.figure()
        plt.rcParams.update({'font.size': 8})
        venn2(subsets=(only_s1, only_s2, both), 
              set_labels=(f"{sample1_name} {feature_type}", f"{sample2_name} {feature_type}"),
              set_colors=(s1.color, s2.color))
        
        if output_file:
            plt.savefig(output_file)
        plt.show()
        
        return overlap
    
    def create_venn3(self, sample1_name, sample2_name, sample3_name, 
                     feature_type='enhancers',
                     output_file=None):
        """
        Create a 3-way Venn diagram comparing three samples.
        Uses colors assigned to each sample.
        
        Parameters:
        -----------
        sample1_name, sample2_name, sample3_name : str
            Names of samples to compare
        feature_type : str
            Feature type to compare
        output_file : str, optional
            Path to save the figure
        """
        s1 = self.samples[sample1_name]
        s2 = self.samples[sample2_name]
        s3 = self.samples[sample3_name]
        
        bed1 = getattr(s1, feature_type)
        bed2 = getattr(s2, feature_type)
        bed3 = getattr(s3, feature_type)
        
        # Calculate overlaps
        overlap12 = bed1.intersect(bed2, f=0.5, r=True, wo=True)
        overlap13 = bed1.intersect(bed3, f=0.5, r=True, wo=True)
        overlap23 = bed2.intersect(bed3, f=0.5, r=True, wo=True)
        overlap123 = bed1.intersect(overlap23, f=0.5, r=True, wo=True)
        
        ABC = overlap123.count()
        ABc = overlap12.count() - ABC
        AbC = overlap13.count() - ABC
        Abc = bed1.count() - (ABc + AbC + ABC)
        aBC = overlap23.count() - ABC
        aBc = bed2.count() - (ABc + aBC + ABC)
        abC = bed3.count() - (AbC + aBC + ABC)
        
        plt.figure()
        plt.rcParams.update({'font.size': 6})
        venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC),
              set_labels=(f"{sample1_name} {feature_type}", 
                         f"{sample2_name} {feature_type}",
                         f"{sample3_name} {feature_type}"),
              set_colors=(s1.color, s2.color, s3.color))
        
        if output_file:
            plt.savefig(output_file)
        plt.show()
    
    def create_scatter_comparison(self, sample1_name, sample2_name, 
                                  feature_type='enhancers',
                                  output_file=None):
        """
        Create a scatter plot comparing odds ratios between two samples.
        Uses the color assigned to sample1.
        
        Parameters:
        -----------
        sample1_name, sample2_name : str
            Names of samples to compare
        feature_type : str
            Feature type to compare
        output_file : str, optional
            Path to save the figure
        
        Returns:
        --------
        dict : Statistics including R-squared and Spearman correlation
        """
        s1 = self.samples[sample1_name]
        s2 = self.samples[sample2_name]
        
        bed1 = getattr(s1, feature_type)
        bed2 = getattr(s2, feature_type)
        
        # Perform intersection
        overlap = bed1.intersect(bed2, f=0.5, r=True, wo=True)
        overlap_df = pd.read_csv(overlap.fn, sep='\t', header=None, low_memory=False)
        overlap_df.insert(0, 'ID', range(len(overlap_df)))
        
        # Filter for positive odds ratios
        filtered = overlap_df[(overlap_df.iloc[:, 5] > 0) & (overlap_df.iloc[:, 10] > 0)]
        
        odds1 = filtered.iloc[:, 5]
        odds2 = filtered.iloc[:, 10]
        
        # Calculate statistics
        slope, intercept, r_value, p_value, std_err = linregress(
            np.log(odds1), np.log(odds2)
        )
        line = np.exp(intercept + slope * np.log(odds1))
        r_squared = r_value**2
        spearman_coef, spearman_p = stats.spearmanr(np.log(odds1), np.log(odds2))
        
        # Create plot using sample1's color
        plt.figure()
        plt.scatter(odds1, odds2, color=s2.color, alpha=0.1)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel(f'{sample1_name} odds ratio')
        plt.ylabel(f'{sample2_name} odds ratio')
        plt.plot(odds1, line, color='maroon', linestyle='dashed', 
                linewidth=2, label='Best Fit Line')
        plt.legend()
        
        if output_file:
            plt.savefig(output_file)
        plt.show()
        
        return {
            'r_value': r_value,
            'r_squared': r_squared,
            'spearman_coefficient': spearman_coef,
            'spearman_p_value': spearman_p
        }
    
    def plot_odds_cdf(self, sample_names, feature_type='super_enhancers',
                      output_file=None):
        """
        Plot cumulative distribution function of odds ratios for multiple samples.
        Uses colors assigned to each sample.
        
        Parameters:
        -----------
        sample_names : list
            List of sample names to include
        feature_type : str
            Feature type to plot
        output_file : str, optional
            Path to save the figure
        """
        fig, ax = plt.subplots(figsize=(8, 5))
        
        for sample_name in sample_names:
            sample = self.samples[sample_name]
            odds = sample.get_odds_scores(feature_type)
            sns.ecdfplot(x=odds, ax=ax, label=f'{sample_name} {feature_type}', 
                        color=sample.color, linewidth=2)
        
        ax.set_xscale('log')
        ax.legend()
        ax.set_title('CDF of Odds Ratios')
        ax.set_xlabel('Odds Ratio')
        ax.set_ylabel('Cumulative Probability')
        
        if output_file:
            plt.savefig(output_file)
        plt.show()
    
        
    def plot_nucleosome_cdf(self, sample_names, feature_type='fire', output_file=None):
            
        """
        Plot cumulative distribution function of nucleosome scores.
    
        For 'fire': plots nucleosome density directly from FIRE regions
        For other features: plots nucleosome density within feature regions (enhancers, super_enhancers, typical_enhancers)
    
        Uses colors assigned to each sample.
    
        Parameters:
        -----------
        sample_names : list or str
            Single sample name or list of sample names to include
        feature_type : str
            Feature type to analyze. Options:
                - 'fire': nucleosome scores from FIRE regions (default)
                - 'enhancers': nucleosome scores within all enhancers
                - 'super_enhancers': nucleosome scores within super enhancers
                - 'typical_enhancers': nucleosome scores within typical enhancers
        output_file : str, optional
            Path to save the figure
    
        Returns:
        --------
        dict : Dictionary mapping sample names to their nucleosome score data
            For 'fire': returns Series of nucleosome scores
            For other features: returns DataFrame with full intersection details
        """
            
        # Handle single sample name
        if isinstance(sample_names, str):
            sample_names = [sample_names]
    
        fig, ax = plt.subplots(figsize=(8, 5))
        results = {}
    
        for sample_name in sample_names:
            sample = self.samples[sample_name]
        
            if feature_type == 'fire':
                # Original FIRE-based nucleosome scoring
                if sample.nucs is None:
                    print(f"Warning: No nucleosome data available for {sample_name}, skipping...")
                    continue
            
                nuc_scores = sample.get_nucleosome_scores()
                results[sample_name] = nuc_scores
            
                sns.ecdfplot(x=nuc_scores, ax=ax, label=f'{sample_name} FIRE', color=sample.color, linewidth=2)
        
            else:
                # Feature-based nucleosome scoring (enhancers, super_enhancers, typical_enhancers)
                if sample.nucs is None:
                    print(f"Warning: No nucleosome data available for {sample_name}, skipping...") 
                    continue
            
                # Get the feature and nucleosome data
                feat = getattr(sample, feature_type)
                nuc = sample.nucs
            
                # Intersect feature with nucleosomes
                feat_nucs = feat.intersect(nuc, wa=True, wb=True)
            
                # Parse the intersection results
                column_names = ['Chrom_f', 'Start_f', 'End_f', 'Name', 'OddsScore', 
                           'Chrom_n', 'Start_n', 'End_n', 'TotalFibers', 'NucFibers', 'NucScore']
                df = feat_nucs.to_dataframe(names=column_names)
            
                results[sample_name] = df
            
                # Plot CDF
                sns.ecdfplot(x=df['NucScore'], ax=ax, label=f'{sample_name} {feature_type}', color=sample.color, linewidth=2)
    
        # Set plot labels based on feature type
        if feature_type == 'fire':
            title = 'CDF of Nucleosome Density in FIRE Regions'
        else:
            title = f'CDF of Nucleosome Density in {feature_type.replace("_", " ").title()}'
    
        ax.legend()
        ax.set_title(title)
        ax.set_xlabel('Nucleosome Score')
        ax.set_ylabel('Cumulative Probability')
    
        if output_file:
            plt.savefig(output_file)
            plt.show()
    
        return results
    
    def v2create_venn3(self, sample1_name, sample2_name, sample3_name, 
                 feature_type='nucs',
                 output_file=None,
                 return_subsets=False):
        """
        Create a 3-way Venn diagram comparing three samples.
        Uses colors assigned to each sample.
    
        Parameters:
        -----------
        sample1_name, sample2_name, sample3_name : str
            Names of samples to compare
        feature_type : str
            Feature type to compare
        output_file : str, optional
            Path to save the figure
        return_subsets : bool, optional
            If True, returns a dictionary with BedTool objects for each subset
    
        Returns:
        --------
        dict (if return_subsets=True) : Dictionary with keys:
            'Abc': Regions unique to sample1
            'aBc': Regions unique to sample2
            'abC': Regions unique to sample3
            'ABc': Regions in sample1 and sample2 only
            'AbC': Regions in sample1 and sample3 only
            'aBC': Regions in sample2 and sample3 only
            'ABC': Regions in all three samples
        """
        s1 = self.samples[sample1_name]
        s2 = self.samples[sample2_name]
        s3 = self.samples[sample3_name]
    
        bed1 = getattr(s1, feature_type)
        bed2 = getattr(s2, feature_type)
        bed3 = getattr(s3, feature_type)
    
        # Calculate overlaps for counting (using bed1 as reference)
        bed1_and_bed2 = bed1.intersect(bed2, f=0.5, r=True, u=True)
        bed1_and_bed3 = bed1.intersect(bed3, f=0.5, r=True, u=True)
        bed1_and_bed2_and_bed3 = bed1.intersect(bed2, f=0.5, r=True, u=True)\
                                  .intersect(bed3, f=0.5, r=True, u=True)
    
        # For bed2 perspective
        bed2_and_bed3 = bed2.intersect(bed3, f=0.5, r=True, u=True)
    
        # Calculate counts for Venn diagram
        ABC = bed1_and_bed2_and_bed3.count()
        ABc = bed1_and_bed2.count() - ABC
        AbC = bed1_and_bed3.count() - ABC
        Abc = bed1.count() - (ABc + AbC + ABC)
        aBC = bed2_and_bed3.count() - ABC
        aBc = bed2.count() - (ABc + aBC + ABC)
        abC = bed3.count() - (AbC + aBC + ABC)
    
        # Plot Venn diagram
        plt.figure()
        plt.rcParams.update({'font.size': 6})
        venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC),
              set_labels=(f"{sample1_name} {feature_type}", 
                          f"{sample2_name} {feature_type}",
                          f"{sample3_name} {feature_type}"),
              set_colors=(s1.color, s2.color, s3.color))
    
        if output_file:
            plt.savefig(output_file)
            plt.show()
    
        # Generate BedTool objects for each subset if requested
        if return_subsets:
            subsets = {}
        
            # ABC: All three samples (from bed1's perspective)
            subsets['ABC'] = bed1_and_bed2_and_bed3
        
            # ABc: bed1 and bed2, but not bed3 (from bed1's perspective)
            subsets['ABc'] = bed1_and_bed2.intersect(bed3, f=0.5, r=True, v=True)
        
            # AbC: bed1 and bed3, but not bed2 (from bed1's perspective)
            subsets['AbC'] = bed1_and_bed3.intersect(bed2, f=0.5, r=True, v=True)
        
            # Abc: bed1 only (not in bed2 or bed3)
            subsets['Abc'] = bed1.intersect(bed2, f=0.5, r=True, v=True)\
                                 .intersect(bed3, f=0.5, r=True, v=True)
        
            # aBC: bed2 and bed3, but not bed1 (from bed2's perspective)
            subsets['aBC'] = bed2_and_bed3.intersect(bed1, f=0.5, r=True, v=True)
        
            # aBc: bed2 only (not in bed1 or bed3)
            subsets['aBc'] = bed2.intersect(bed1, f=0.5, r=True, v=True)\
                                 .intersect(bed3, f=0.5, r=True, v=True)
        
            # abC: bed3 only (not in bed1 or bed2)
            subsets['abC'] = bed3.intersect(bed1, f=0.5, r=True, v=True)\
                                 .intersect(bed2, f=0.5, r=True, v=True)
        
            return subsets

    def v3create_venn3(self, sample1_name, sample2_name, sample3_name, 
                 feature_type='nucs',
                 output_file=None,
                 return_subsets=False):
        """
        Create a 3-way Venn diagram comparing three samples.
        Uses colors assigned to each sample.
    
        Parameters:
        -----------
        sample1_name, sample2_name, sample3_name : str
            Names of samples to compare
        feature_type : str
            Feature type to compare
        output_file : str, optional
            Path to save the figure
        return_subsets : bool, optional
            If True, returns a dictionary with BedTool objects for each subset
    
        Returns:
        --------
        dict (if return_subsets=True) : Dictionary with keys:
            'Abc': Regions unique to sample1
            'aBc': Regions unique to sample2
            'abC': Regions unique to sample3
            'ABc': Regions in sample1 and sample2 only
            'AbC': Regions in sample1 and sample3 only
            'aBC': Regions in sample2 and sample3 only
            'ABC': Regions in all three samples
        """
        s1 = self.samples[sample1_name]
        s2 = self.samples[sample2_name]
        s3 = self.samples[sample3_name]
    
        bed1 = getattr(s1, feature_type)
        bed2 = getattr(s2, feature_type)
        bed3 = getattr(s3, feature_type)
    
        # Calculate overlaps using u=True for consistent counting
        bed1_in_bed2 = bed1.intersect(bed2, f=0.5, r=True, u=True)
        bed1_in_bed3 = bed1.intersect(bed3, f=0.5, r=True, u=True)
        bed2_in_bed3 = bed2.intersect(bed3, f=0.5, r=True, u=True)
    
        # Get three-way overlap
        bed1_in_bed2_and_bed3 = bed1_in_bed2.intersect(bed3, f=0.5, r=True, u=True)
    
        # Calculate counts for Venn diagram
        ABC = bed1_in_bed2_and_bed3.count()
        ABc = bed1_in_bed2.count() - ABC
        AbC = bed1_in_bed3.count() - ABC
        Abc = bed1.count() - (ABc + AbC + ABC)
        aBC = bed2_in_bed3.count() - ABC
        aBc = bed2.count() - (ABc + aBC + ABC)
        abC = bed3.count() - (AbC + aBC + ABC)
    
        # Plot Venn diagram
        plt.figure()
        plt.rcParams.update({'font.size': 6})
        venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC),
              set_labels=(f"{sample1_name} {feature_type}", 
                          f"{sample2_name} {feature_type}",
                          f"{sample3_name} {feature_type}"),
              set_colors=(s1.color, s2.color, s3.color))
    
        if output_file:
            plt.savefig(output_file)
            plt.show()
    
        # Generate BedTool objects for each subset if requested
        if return_subsets:
            subsets = {}
        
            # ABC: Regions in all three
            subsets['ABC'] = bed1_in_bed2_and_bed3
        
            # ABc: Regions in bed1 AND bed2 but NOT bed3
            subsets['ABc'] = bed1_in_bed2.intersect(bed3, f=0.5, r=True, v=True)
        
            # AbC: Regions in bed1 AND bed3 but NOT bed2
            subsets['AbC'] = bed1_in_bed3.intersect(bed2, f=0.5, r=True, v=True)
        
            # Abc: Regions ONLY in bed1
            subsets['Abc'] = bed1.intersect(bed2, f=0.5, r=True, v=True)\
                                 .intersect(bed3, f=0.5, r=True, v=True)
        
            # aBC: Regions in bed2 AND bed3 but NOT bed1
            subsets['aBC'] = bed2_in_bed3.intersect(bed1, f=0.5, r=True, v=True)
        
            # aBc: Regions ONLY in bed2
            subsets['aBc'] = bed2.intersect(bed1, f=0.5, r=True, v=True)\
                                 .intersect(bed3, f=0.5, r=True, v=True)
        
            # abC: Regions ONLY in bed3
            subsets['abC'] = bed3.intersect(bed1, f=0.5, r=True, v=True)\
                                 .intersect(bed2, f=0.5, r=True, v=True)
        
            return subsets

    # ------------------------------------------------------------------ #
    #  Hockey stick plot                                                   #
    # ------------------------------------------------------------------ #

    @staticmethod
    def _find_inflection(ranks, odds):
        """
        Find the inflection point of the hockey stick curve using the
        maximum perpendicular distance from the line connecting the first
        and last points (same geometric approach used in ROSE for
        super-enhancer calling).

        Parameters:
        -----------
        ranks : np.ndarray - x values (descending rank)
        odds  : np.ndarray - y values (odds ratio)

        Returns:
        --------
        int : index into ranks/odds where the inflection occurs
        """
        # Normalise both axes to [0,1] so distances are comparable
        x = (ranks - ranks.min()) / (ranks.max() - ranks.min())
        y = (odds  - odds.min())  / (odds.max()  - odds.min())

        dx = x[-1] - x[0]
        dy = y[-1] - y[0]
        line_len = np.sqrt(dx**2 + dy**2)

        distances = np.abs(dy * x - dx * y + x[-1] * y[0] - y[-1] * x[0]) / line_len
        return int(np.argmax(distances))

    def plot_hockey_stick(self, sample_name, gtf_path,
                          tad_path,
                          gene_names=None,
                          feature_type='enhancers',
                          top_only=False,
                          output_file=None):
        """
        Plot a ranked co-accessibility hockey stick graph for one sample.

        Ranks all enhancer pairs in the ce_rank file by odds ratio. The
        curve is split at the automatically detected inflection point:
        the flat portion (low odds) is drawn in black and the steep
        portion (high odds, super-enhancer territory) is drawn in red.

        If gene_names are provided, enhancers are assigned to genes using
        the TAD-based method (sample.assign_genes): each enhancer is mapped
        to the closest protein-coding gene within the same TAD, falling back
        to the globally closest gene if no protein-coding gene shares a TAD.
        Enhancers assigned to the requested gene(s) are then marked as
        labelled dots on the curve.

        Parameters:
        -----------
        sample_name : str
            Key into self.samples for the sample to plot.
        gtf_path : str
            Path to mm10 GTF annotation file (protein-coding genes only used).
        tad_path : str
            Path to mm10 TAD BED file (3-column: chrom, start, end).
            TADs originally in mm9 (e.g. Dowen et al.) should be lifted to
            mm10 with UCSC liftOver before passing here.
        gene_names : list of str, optional
            One or more gene names to highlight (e.g. ['Sox2', 'Mir290a']).
            Each gene gets a distinct colour and labelled dots at all
            enhancers assigned to it via the TAD method.
        feature_type : str, optional
            Which enhancer set to use for gene assignment.
            One of 'enhancers', 'super_enhancers', 'typical_enhancers',
            'intergenic_enhancers'. Default: 'enhancers'.
        top_only : bool, optional
            If True, only the single highest-odds ce pair for each gene
            is marked on the graph. If False, all assigned enhancers are
            marked. Default: False.
        output_file : str, optional
            Path to save the figure (e.g. 'V6.5_hockeystick.svg').

        Returns:
        --------
        pandas.DataFrame : ce_rank data with an added 'rank' column.
        dict : { gene_name: region_id } of the top-ranked enhancer per gene
               (if top_only=True) or all region IDs (if top_only=False).

        Example:
        --------
        analysis.plot_hockey_stick(
            sample_name  = 'V6.5',
            gtf_path     = '/path/to/gencode.vM23.basic.annotation.gtf',
            tad_path     = '/path/to/mESC_TADs_mm10.bed',
            gene_names   = ['Sox2', 'Mir290a'],
            feature_type = 'intergenic_enhancers',
            top_only     = True,
            output_file  = 'V6.5_hockeystick.svg'
        )
        """
        sample = self.samples[sample_name]

        # ---- 1. Load and rank ce_rank file -----------------------------
        col_names = ['Column1', 'Column2', 'Column3', 'Column4', 'Column5', 'Column6']
        df = pd.read_csv(sample.ce_rank_path, sep='\t', header=None, names=col_names, low_memory=False)

        # Sort descending by odds ratio, assign rank (1 = highest odds)
        df = df.sort_values('Column4', ascending=False).reset_index(drop=True)
        df['rank'] = np.arange(1, len(df) + 1)   # rank 1 = highest odds ratio

        ranks = df['rank'].values.astype(float)
        odds  = df['Column4'].values.astype(float)

        # ---- 2. Detect inflection point --------------------------------
        inflection_idx  = self._find_inflection(ranks, odds)
        inflection_rank = ranks[inflection_idx]

        # ---- 3. Find gene-associated enhancers via TAD-based assign_genes()
        gene_enhancer_ranks = {}   # { gene_name: [(rank, odds), ...] }
        gene_region_ids     = {}   # { gene_name: [region_id, ...] }

        if gene_names:
            from pybedtools import BedTool
            import re as _re

            print(f"[{sample_name}] Finding enhancers for requested genes via TAD method...")

            # Build a set of all ce_rank region IDs for fast matching
            ce_rank_ids = set(df['Column1'].tolist() + df['Column2'].tolist())

            # Load TADs and enhancer feature once
            tad_bed      = BedTool(tad_path).sort()
            enhancer_bed = getattr(sample, feature_type)

            # Parse gene coordinates from GTF for requested genes only
            # All biotypes included so miRNAs, lncRNAs etc. are findable
            print(f"[{sample_name}] Parsing GTF for requested gene coordinates...")
            targets     = {g.strip().lower() for g in gene_names}
            gene_coords = {g.strip().lower(): [] for g in gene_names}

            with open(gtf_path, 'r') as _f:
                for _line in _f:
                    if _line.startswith('#'):
                        continue
                    _fields = _line.strip().split('\t')
                    if len(_fields) < 9 or _fields[2] != 'gene':
                        continue
                    _attrs = dict(_re.findall(r'(\w+)\s+"([^"]+)"', _fields[8]))
                    _gname = _attrs.get('gene_name', '').strip().lower()
                    if _gname in targets:
                        gene_coords[_gname].append((
                            _fields[0],
                            int(_fields[3]) - 1,
                            int(_fields[4])
                        ))

            for gene in gene_names:
                key       = gene.strip().lower()
                intervals = gene_coords.get(key, [])

                if not intervals:
                    print(f"  Warning: '{gene}' not found in GTF.")
                    gene_enhancer_ranks[gene] = []
                    gene_region_ids[gene]     = []
                    continue

                # Build BedTool for this gene
                gene_str = '\n'.join(
                    f"{chrom}\t{start}\t{end}"
                    for chrom, start, end in intervals
                )
                gene_bed = BedTool(gene_str, from_string=True).sort()

                # Find TAD(s) containing this gene
                gene_tads = tad_bed.intersect(gene_bed, u=True)

                if gene_tads.count() == 0:
                    # Gene not in any TAD - fall back to 500kb window
                    print(f"  '{gene}': not in any TAD, falling back to 500kb window.")
                    closest = enhancer_bed.closest(gene_bed, d=True, t='all')
                    rows = []
                    for interval in closest:
                        fields = str(interval).strip().split('\t')
                        try:
                            dist = int(fields[-1])
                        except ValueError:
                            continue
                        if 0 <= dist <= 100_000:
                            rows.append(f"{fields[0]}:{fields[1]}-{fields[2]}")
                    matched_ids = set(rows) & ce_rank_ids
                else:
                    # Get all enhancers within the same TAD(s) as the gene
                    enh_in_tad = enhancer_bed.intersect(gene_tads, u=True)
                    if enh_in_tad.count() == 0:
                        print(f"  '{gene}': TAD found but no enhancers inside, "
                              f"falling back to 500kb window.")
                        closest = enhancer_bed.closest(gene_bed, d=True, t='all')
                        rows = []
                        for interval in closest:
                            fields = str(interval).strip().split('\t')
                            try:
                                dist = int(fields[-1])
                            except ValueError:
                                continue
                            if 0 <= dist <= 100_000:
                                rows.append(f"{fields[0]}:{fields[1]}-{fields[2]}")
                        matched_ids = set(rows) & ce_rank_ids
                    else:
                        rows = []
                        for interval in enh_in_tad:
                            fields = str(interval).strip().split('\t')
                            rows.append(f"{fields[0]}:{fields[1]}-{fields[2]}")
                        matched_ids = set(rows) & ce_rank_ids

                if not matched_ids:
                    print(f"  Warning: enhancers for '{gene}' not found in ce_rank.")
                    gene_enhancer_ranks[gene] = []
                    gene_region_ids[gene]     = []
                    continue

                mask       = df['Column1'].isin(matched_ids) | df['Column2'].isin(matched_ids)
                matched_df = df[mask]

                pts = list(zip(matched_df['rank'], matched_df['Column4']))
                gene_enhancer_ranks[gene] = pts
                gene_region_ids[gene]     = sorted(matched_ids)
                print(f"  '{gene}': {len(pts)} enhancer(s) found in TAD")

        # ---- 4. Plot ---------------------------------------------------
        fig, ax = plt.subplots(figsize=(10, 7))

        # Black: flat portion — low odds, high rank numbers (left side after inversion)
        black_mask = ranks > inflection_rank
        ax.plot(ranks[black_mask], odds[black_mask], color='black', linewidth=1.5)

        # Red: steep portion — high odds, low rank numbers (right side after inversion)
        red_mask = ranks <= inflection_rank
        ax.plot(ranks[red_mask], odds[red_mask], color='red', linewidth=1.5)

        # Gene-associated enhancer dots
        dot_colors = plt.cm.tab10.colors
        top_region_ids = {}   # { gene_name: region_id of highest-odds pair }
        for i, (gene, pts) in enumerate(gene_enhancer_ranks.items()):
            if not pts:
                continue
            dot_color = dot_colors[i % len(dot_colors)]

            # If top_only, reduce to just the single highest-odds pair
            pts_to_plot = [max(pts, key=lambda p: p[1])] if top_only else pts

            pt_ranks = [p[0] for p in pts_to_plot]
            pt_odds  = [p[1] for p in pts_to_plot]

            ax.scatter(pt_ranks, pt_odds, color=dot_color, s=70, zorder=5)

            # Label at the plotted point (only one if top_only, else highest)
            best_r, best_o = max(pts_to_plot, key=lambda p: p[1])

            # Use textcoords='offset points' so the label offset is in
            # display units (pixels), not data units — this keeps labels
            # off the curve regardless of where on the axes the dot lands
            ax.annotate(
                gene,
                xy=(best_r, best_o),
                xytext=(40, 40),
                textcoords='offset points',
                fontsize=11,
                fontweight='bold',
                color=dot_color,
                path_effects=[pe.withStroke(linewidth=2, foreground='white')]
            )

            # Store the region ID of the top pair for return
            top_df = df[df['rank'] == best_r]
            if not top_df.empty:
                row = top_df.iloc[0]
                top_region_ids[gene] = {
                    'Column1': row['Column1'],
                    'Column2': row['Column2'],
                    'odds_ratio': best_o
                }

        # Axes formatting — invert x so rank 1 (highest odds) is on the right
        ax.invert_xaxis()
        ax.set_xlabel('Co-Accessible Regions Ranked By Odds Ratio', fontsize=13)
        ax.set_ylabel('Odds Ratio', fontsize=13)
        ax.set_title(
            f'Ranked Co-Accessibility — {sample_name}',
            fontsize=14
        )

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file)
            print(f"Saved: {output_file}")

        plt.show()

        if gene_names:
            return df, top_region_ids if top_only else gene_region_ids
        return df
