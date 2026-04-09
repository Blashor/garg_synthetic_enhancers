#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GenomicSample class for representing a single genomic sample.
@author: michelemeline
"""

import os
import re
import subprocess
import tempfile
import pandas as pd
from pybedtools import BedTool



class GenomicSample:
    """
    Class to represent a single genomic sample with associated files and analysis methods.
    """
    
    def __init__(self, name, ce_rank_path, fire_peaks_intergenic_path=None,
                 fire_bed_path=None, nucs_bed_path=None, color='blue',
                 chain_file=None, liftover_path='liftOver', is_hg38=False):
        """
        Initialize a genomic sample.
        
        Parameters:
        -----------
        name : str
            Sample identifier (e.g., 'Exo+Dox', 'Exo-24h', etc.)
        ce_rank_path : str
            Path to ce_rank.txt file (tab-separated, no header)
        fire_peaks_intergenic_path : str, optional
            Path to FIRE_peaks_intergenic.bed file. When provided, all intergenic
            regions are included in enhancers; those without a scored pair get
            odds ratio of 0.
        fire_bed_path : str, optional
            Path to FIRE peaks BED file
        nucs_bed_path : str, optional
            Path to nucleosome regions BED file
        color : str, optional
            Color to use for this sample in visualizations (default: 'blue')
        chain_file : str, optional
            Path to liftOver chain file (e.g. mm10ToHg38.over.chain.gz).
            Required to use any hg38_* properties on mm10 samples.
        liftover_path : str, optional
            Path to liftOver executable. Defaults to 'liftOver' assuming it
            is on PATH. Provide full path if needed (e.g. '~/ucsc/liftOver').
        is_hg38 : bool, optional
            Set True if sample is already hg38 (e.g. GM12878, K562).
            All hg38_* properties will return the base feature directly.
            (default: False)
        """
        self.name = name
        self.ce_rank_path = ce_rank_path
        self.fire_peaks_intergenic_path = fire_peaks_intergenic_path
        self.fire_bed_path = fire_bed_path
        self.nucs_bed_path = nucs_bed_path
        self.color = color
        self.chain_file = chain_file
        self.liftover_path = liftover_path
        self.is_hg38 = is_hg38

        # Lazy loading - only load when needed
        self._enhancers = None
        self._intergenic_enhancers = None
        self._fire = None
        self._nucs = None
        self._super_enhancers = None
        self._typical_enhancers = None

        # Lazy loading - hg38 lifted versions
        self._hg38_enhancers = None
        self._hg38_intergenic_enhancers = None
        self._hg38_fire = None
        self._hg38_super_enhancers = None
        self._hg38_typical_enhancers = None
    
    @property
    def enhancers(self):
        """Get or create enhancers BedTool."""
        if self._enhancers is None:
            self._enhancers = self._extract_enhancers()
        return self._enhancers
    
    @property
    def fire(self):
        """Get FIRE peaks BedTool."""
        if self._fire is None and self.fire_bed_path:
            self._fire = BedTool(self.fire_bed_path)
        return self._fire
    
    @property
    def nucs(self):
        """Get nucleosome regions BedTool."""
        if self._nucs is None and self.nucs_bed_path:
            self._nucs = BedTool(self.nucs_bed_path)
        return self._nucs
    
    @property
    def super_enhancers(self):
        """Get super enhancers BedTool."""
        if self._super_enhancers is None:
            self._super_enhancers = self._extract_super_enhancers()
        return self._super_enhancers
    
    @property
    def typical_enhancers(self):
        """Get typical enhancers BedTool."""
        if self._typical_enhancers is None:
            self._typical_enhancers = self._extract_typical_enhancers()
        return self._typical_enhancers

    @property
    def intergenic_enhancers(self):
        """Get all intergenic enhancers BedTool, including unscored regions (score=0)."""
        if self._intergenic_enhancers is None:
            self._intergenic_enhancers = self._extract_intergenic_enhancers()
        return self._intergenic_enhancers

    # ------------------------------------------------------------------ #
    #  hg38 properties - pass-through if is_hg38, else run liftOver      #
    # ------------------------------------------------------------------ #

    @property
    def hg38_enhancers(self):
        """Get enhancers in hg38 coordinates. Pass-through if sample is already hg38."""
        if self._hg38_enhancers is None:
            self._hg38_enhancers = self.enhancers if self.is_hg38 else self.liftover(self.enhancers, 'enhancers')
        return self._hg38_enhancers

    @property
    def hg38_intergenic_enhancers(self):
        """Get intergenic enhancers in hg38 coordinates. Pass-through if sample is already hg38."""
        if self._hg38_intergenic_enhancers is None:
            self._hg38_intergenic_enhancers = self.intergenic_enhancers if self.is_hg38 else self.liftover(self.intergenic_enhancers, 'intergenic_enhancers')
        return self._hg38_intergenic_enhancers

    @property
    def hg38_super_enhancers(self):
        """Get super enhancers in hg38 coordinates. Pass-through if sample is already hg38."""
        if self._hg38_super_enhancers is None:
            self._hg38_super_enhancers = self.super_enhancers if self.is_hg38 else self.liftover(self.super_enhancers, 'super_enhancers')
        return self._hg38_super_enhancers

    @property
    def hg38_typical_enhancers(self):
        """Get typical enhancers in hg38 coordinates. Pass-through if sample is already hg38."""
        if self._hg38_typical_enhancers is None:
            self._hg38_typical_enhancers = self.typical_enhancers if self.is_hg38 else self.liftover(self.typical_enhancers, 'typical_enhancers')
        return self._hg38_typical_enhancers

    @property
    def hg38_fire(self):
        """Get FIRE peaks in hg38 coordinates. Pass-through if sample is already hg38."""
        if self._hg38_fire is None:
            self._hg38_fire = self.fire if self.is_hg38 else self.liftover(self.fire, 'fire')
        return self._hg38_fire

    def liftover(self, feature, feature_name):
        """
        Lift over a BedTool feature to hg38 using UCSC liftOver.

        Parameters:
        -----------
        feature : BedTool
            The BedTool object to lift over (e.g. self.enhancers, self.fire).
        feature_name : str
            Short label used for temp/output file naming (e.g. 'enhancers').

        Returns:
        --------
        BedTool : Lifted-over regions in hg38 coordinates. Unmapped regions
                  are written to a temp unlifted file and silently discarded.

        Raises:
        -------
        ValueError : If chain_file is not set on this sample.
        RuntimeError : If liftOver exits with a non-zero return code.
        """
        if self.chain_file is None:
            raise ValueError(
                f"chain_file not set for sample '{self.name}'. "
                "Provide a chain_file path to use liftover()."
            )

        chain = os.path.expanduser(self.chain_file)
        liftover_exe = os.path.expanduser(self.liftover_path)

        # Write BedTool to a named temp file liftOver can read
        with tempfile.NamedTemporaryFile(
                suffix=f'_{self.name}_{feature_name}_input.bed',
                delete=False, mode='w') as tmp_in:
            tmp_in.write(str(feature))
            input_path = tmp_in.name

        output_path   = input_path.replace('_input.bed', '_hg38.bed')
        unlifted_path = input_path.replace('_input.bed', '_unlifted.bed')

        try:
            result = subprocess.run(
                [liftover_exe, input_path, chain, output_path, unlifted_path],
                capture_output=True, text=True
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"liftOver failed for {self.name} {feature_name}:\n{result.stderr}"
                )
            print(f"{self.name} {feature_name}: liftOver complete. "
                  f"Unlifted regions: {unlifted_path}")
            return BedTool(output_path).sort()
        finally:
            os.remove(input_path)

    def _read_ce_rank(self):
        """Read and parse ce_rank.txt file (tab-separated, no header)."""
        column_names = ['Column1', 'Column2', 'Column3', 'Column4', 'Column5', 'Column6']
        return pd.read_csv(self.ce_rank_path, sep='\t', header=None, names=column_names, low_memory=False)

    def _get_max_odds_per_enhancer(self, data):
        """
        Given a ce_rank dataframe, return a dataframe with one row per unique
        enhancer ID and its maximum odds ratio, considering both Column1 and Column2.
        """
        # Max odds for enhancers appearing in Column1
        idx1 = data.groupby('Column1')['Column4'].idxmax()
        df1 = data.loc[idx1][['Column1', 'Column4', 'Column3']].rename(
            columns={'Column1': 'EnhID', 'Column3': 'SuperID'})

        # Max odds for enhancers appearing in Column2
        idx2 = data.groupby('Column2')['Column4'].idxmax()
        df2 = data.loc[idx2][['Column2', 'Column4', 'Column3']].rename(
            columns={'Column2': 'EnhID', 'Column3': 'SuperID'})

        # Combine, reset index so idxmax works on unique indices, then take overall max per enhancer
        combined = pd.concat([df1, df2]).reset_index(drop=True)
        idx_max = combined.groupby('EnhID')['Column4'].idxmax()
        return combined.loc[idx_max].drop_duplicates(subset='EnhID')
    
    def _extract_enhancers(self):
        """Extract enhancers from ce_rank file using max odds across Column1 and Column2."""
        data = self._read_ce_rank()
        df = self._get_max_odds_per_enhancer(data)

        splitcol = df["EnhID"].str.split(pat=":", expand=True)
        splitpos = splitcol[1].str.split(pat="-", expand=True)

        FDRenhID = pd.DataFrame({
            "Chrom": splitcol[0],
            "Start": splitpos[0].astype(int),
            "End": splitpos[1].astype(int),
            "Name": df["SuperID"].values,
            "Score": df["Column4"].values
        })

        FDRenh_string = FDRenhID.to_csv(sep="\t", index=False, header=False)
        return BedTool(FDRenh_string, from_string=True).sort()
    
    def _extract_super_enhancers(self):
        """Extract super enhancers (Column6 == 'Super') using max odds across Column1 and Column2."""
        data = self._read_ce_rank()
        se_data = data[data['Column6'] == 'Super']
        df = self._get_max_odds_per_enhancer(se_data)

        splitcol = df["EnhID"].str.split(pat=":", expand=True)
        splitpos = splitcol[1].str.split(pat="-", expand=True)

        FDRse = pd.DataFrame({
            "Chrom": splitcol[0],
            "Start": splitpos[0].astype(int),
            "End": splitpos[1].astype(int),
            "Name": df["SuperID"].values,
            "Score": df["Column4"].values
        })

        FDRse_string = FDRse.to_csv(sep="\t", index=False, header=False)
        return BedTool(FDRse_string, from_string=True).sort()
    
    def _extract_typical_enhancers(self):
        """Extract typical enhancers (non-super) using max odds across Column1 and Column2."""
        data = self._read_ce_rank()
        te_data = data[data['Column6'] != 'Super']
        df = self._get_max_odds_per_enhancer(te_data)

        splitcol = df["EnhID"].str.split(pat=":", expand=True)
        splitpos = splitcol[1].str.split(pat="-", expand=True)

        FDRte = pd.DataFrame({
            "Chrom": splitcol[0],
            "Start": splitpos[0].astype(int),
            "End": splitpos[1].astype(int),
            "Name": df["SuperID"].values,
            "Score": df["Column4"].values
        })

        FDRte_string = FDRte.to_csv(sep="\t", index=False, header=False)
        return BedTool(FDRte_string, from_string=True).sort()

    def _extract_intergenic_enhancers(self):
        """
        Extract all intergenic enhancers, merging scored pairs with the full
        FIRE_peaks_intergenic.bed. Regions not present in any scored pair
        are assigned an odds ratio of 0. Requires fire_peaks_intergenic_path.
        """
        if self.fire_peaks_intergenic_path is None:
            raise ValueError(
                f"fire_peaks_intergenic_path not set for sample '{self.name}'. "
                "Provide this path to use intergenic_enhancers."
            )

        data = self._read_ce_rank()
        df = self._get_max_odds_per_enhancer(data)

        # Load all intergenic regions as the anchor
        intergenic = pd.read_csv(
            self.fire_peaks_intergenic_path, sep='\t', header=None,
            names=['Chrom', 'Start', 'End']
        )
        intergenic['EnhID'] = (intergenic['Chrom'] + ':'
                                + intergenic['Start'].astype(str) + '-'
                                + intergenic['End'].astype(str))

        # Left merge so all intergenic regions are kept
        merged = intergenic.merge(df[['EnhID', 'Column4', 'SuperID']],
                                  on='EnhID', how='left')
        merged['Column4'] = merged['Column4'].fillna(0)
        merged['SuperID'] = merged['SuperID'].fillna('.')

        FDRinterg = pd.DataFrame({
            "Chrom": merged['Chrom'],
            "Start": merged['Start'],
            "End": merged['End'],
            "Name": merged['SuperID'],
            "Score": merged['Column4']
        })

        FDRinterg_string = FDRinterg.to_csv(sep='\t', index=False, header=False)
        return BedTool(FDRinterg_string, from_string=True).sort()
    
    def get_odds_scores(self, feature_type='enhancers'):
        """
        Get odds ratio scores for a given feature type.
        
        Parameters:
        -----------
        feature_type : str
            One of 'enhancers', 'super_enhancers', 'typical_enhancers'
        
        Returns:
        --------
        pandas.Series : Odds ratio scores
        """
        feature_map = {
            'enhancers': self.enhancers,
            'super_enhancers': self.super_enhancers,
            'typical_enhancers': self.typical_enhancers
        }
        
        bed = feature_map.get(feature_type)
        if bed is None:
            raise ValueError(f"Unknown feature type: {feature_type}")
        
        df = bed.to_dataframe()
        return df['score']
    
    def get_nucleosome_scores(self):
        """Get nucleosome density scores from FIRE regions."""
        if self.nucs is None:
            raise ValueError(f"No nucleosome data available for {self.name}")
        
        column_names = ['Chrom', 'Start', 'End', 'TotalFibers', 'NucFibers', 'NucScore']
        df = self.nucs.to_dataframe(names=column_names)
        return df['NucScore']

    def assign_genes(self, gtf_path, tad_path, feature_type='enhancers'):
        """
        Assign each enhancer to the most proximal protein-coding gene within
        the same topologically associated domain (TAD).

        Strategy (mirrors Dowen et al. / Young lab convention):
          1. Parse protein-coding gene coordinates from the GTF.
          2. Assign enhancers and genes to TADs via bedtools intersect.
          3. For each enhancer, find the closest gene that shares the same TAD
             using bedtools closest restricted to the TAD-gene set.
          4. Enhancers with no gene in their TAD fall back to the globally
             closest protein-coding gene (no distance limit).

        Parameters:
        -----------
        gtf_path : str
            Path to mm10 GTF annotation file (e.g. gencode.vM23.basic.annotation.gtf).
            All gene biotypes are included (protein-coding, miRNA, lncRNA, etc.).
        tad_path : str
            Path to mm10 TAD BED file (3-column: chrom, start, end).
            TADs originally in mm9 should be lifted to mm10 with UCSC liftOver
            (chain: mm9ToMm10.over.chain.gz) before passing here.
        feature_type : str
            Feature to annotate. One of:
                'fire', 'enhancers', 'super_enhancers',
                'typical_enhancers', 'intergenic_enhancers'
            Defaults to 'enhancers'.

        Returns:
        --------
        pandas.DataFrame with columns:
            chrom, start, end, name, score,      <- original enhancer columns
            gene_chrom, gene_start, gene_end,
            gene_id, gene_name, gene_strand,
            distance,                             <- distance to assigned gene (bp)
            assignment_method                     <- 'TAD' or 'nearest_fallback'

        Example:
        --------
        df = samples['V6.5'].assign_genes(
            gtf_path  = '/path/to/gencode.vM23.basic.annotation.gtf',
            tad_path  = '/path/to/mESC_TADs_mm10.bed',
            feature_type = 'intergenic_enhancers'
        )
        """
        # ------------------------------------------------------------------ #
        # 0. Validate inputs                                                  #
        # ------------------------------------------------------------------ #
        valid = ('fire', 'enhancers', 'super_enhancers',
                 'typical_enhancers', 'intergenic_enhancers')
        if feature_type not in valid:
            raise ValueError(
                f"feature_type must be one of {valid}, got '{feature_type}'"
            )
        feature_bed = getattr(self, feature_type)
        if feature_bed is None:
            raise ValueError(
                f"Feature '{feature_type}' is not available for sample '{self.name}'. "
                "Check that the required file paths are set."
            )

        # ------------------------------------------------------------------ #
        # 1. Parse protein-coding genes from GTF                              #
        # ------------------------------------------------------------------ #
        print(f"[{self.name}] Parsing protein-coding genes from GTF: {gtf_path}")
        rows = []
        with open(gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                attrs = _parse_gtf_attributes(fields[8])
                chrom  = fields[0]
                start  = int(fields[3]) - 1   # GTF 1-based -> BED 0-based
                end    = int(fields[4])
                strand = fields[6]
                gene_id   = attrs.get('gene_id', '.')
                gene_name = attrs.get('gene_name', gene_id)
                rows.append(f"{chrom}\t{start}\t{end}\t{gene_id}\t{gene_name}\t{strand}")

        gene_bed = BedTool('\n'.join(rows), from_string=True).sort()
        print(f"[{self.name}] {len(rows)} protein-coding genes loaded.")

        # ------------------------------------------------------------------ #
        # 2. Load TADs                                                        #
        # ------------------------------------------------------------------ #
        tad_bed = BedTool(tad_path).sort()
        print(f"[{self.name}] {tad_bed.count()} TADs loaded from {tad_path}")

        # ------------------------------------------------------------------ #
        # 3. Assign TAD identity to enhancers and genes                       #
        #    bedtools intersect -wa -wb: keep original interval + TAD coords  #
        # ------------------------------------------------------------------ #
        # Enhancers -> TADs: columns [enh(0-4), tad_chrom, tad_start, tad_end]
        enh_in_tad  = feature_bed.intersect(tad_bed, wa=True, wb=True)

        # Genes -> TADs: columns [gene(0-5), tad_chrom, tad_start, tad_end]
        gene_in_tad = gene_bed.intersect(tad_bed, wa=True, wb=True)

        # Parse into DataFrames and create a TAD key (chrom:start-end)
        enh_cols  = ['e_chrom','e_start','e_end','e_name','e_score',
                     'tad_chrom','tad_start','tad_end']
        gene_cols = ['g_chrom','g_start','g_end','gene_id','gene_name','gene_strand',
                     'tad_chrom','tad_start','tad_end']

        enh_df  = enh_in_tad.to_dataframe(names=enh_cols,  comment='#')
        gene_df = gene_in_tad.to_dataframe(names=gene_cols, comment='#')

        enh_df['tad_key']  = (enh_df['tad_chrom']  + ':'
                               + enh_df['tad_start'].astype(str)  + '-'
                               + enh_df['tad_end'].astype(str))
        gene_df['tad_key'] = (gene_df['tad_chrom'] + ':'
                               + gene_df['tad_start'].astype(str) + '-'
                               + gene_df['tad_end'].astype(str))

        print(f"[{self.name}] {len(enh_df)} enhancers fall within a TAD.")
        print(f"[{self.name}] {len(gene_df)} gene-TAD assignments.")

        # ------------------------------------------------------------------ #
        # 4. For each TAD, find the closest gene to each enhancer             #
        # ------------------------------------------------------------------ #
        result_rows = []

        # Group genes by TAD for fast lookup
        genes_by_tad = gene_df.groupby('tad_key')

        # Build an enhancer key for later fallback bookkeeping
        enh_df['enh_key'] = (enh_df['e_chrom'] + ':'
                              + enh_df['e_start'].astype(str) + '-'
                              + enh_df['e_end'].astype(str))

        assigned_enh_keys = set()

        for tad_key, enh_group in enh_df.groupby('tad_key'):
            if tad_key not in genes_by_tad.groups:
                # No protein-coding gene in this TAD — handled in fallback
                continue

            tad_genes = genes_by_tad.get_group(tad_key)

            # Build temporary BedTools for this TAD's enhancers and genes
            enh_str = '\n'.join(
                f"{r.e_chrom}\t{r.e_start}\t{r.e_end}\t{r.e_name}\t{r.e_score}"
                for _, r in enh_group.iterrows()
            )
            gene_str = '\n'.join(
                f"{r.g_chrom}\t{r.g_start}\t{r.g_end}\t{r.gene_id}\t{r.gene_name}\t{r.gene_strand}"
                for _, r in tad_genes.iterrows()
            )

            tad_enh_bed  = BedTool(enh_str,  from_string=True).sort()
            tad_gene_bed = BedTool(gene_str, from_string=True).sort()

            closest = tad_enh_bed.closest(tad_gene_bed, d=True, t='first')

            col_names = ['chrom','start','end','name','score',
                         'gene_chrom','gene_start','gene_end',
                         'gene_id','gene_name','gene_strand','distance']
            df_close = closest.to_dataframe(names=col_names, comment='#')
            df_close['distance'] = pd.to_numeric(df_close['distance'], errors='coerce')
            df_close['assignment_method'] = 'TAD'

            result_rows.append(df_close)

            # Track which enhancers were successfully assigned
            for _, r in df_close.iterrows():
                assigned_enh_keys.add(
                    f"{r['chrom']}:{r['start']}-{r['end']}"
                )

        # ------------------------------------------------------------------ #
        # 5. Fallback: enhancers with no gene in their TAD                    #
        #    (includes enhancers outside all TADs entirely)                   #
        # ------------------------------------------------------------------ #

        # Identify all enhancer keys from the original feature set
        all_enh_df = feature_bed.to_dataframe(
            names=['chrom','start','end','name','score'], comment='#'
        )
        all_enh_df['enh_key'] = (all_enh_df['chrom'] + ':'
                                  + all_enh_df['start'].astype(str) + '-'
                                  + all_enh_df['end'].astype(str))

        fallback_df = all_enh_df[~all_enh_df['enh_key'].isin(assigned_enh_keys)]
        n_fallback  = len(fallback_df)

        if n_fallback > 0:
            print(f"[{self.name}] {n_fallback} enhancers have no gene in their TAD "
                  f"— falling back to globally closest protein-coding gene.")

            fb_str = '\n'.join(
                f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['name']}\t{r['score']}"
                for _, r in fallback_df.iterrows()
            )
            fb_bed  = BedTool(fb_str, from_string=True).sort()
            closest = fb_bed.closest(gene_bed, d=True, t='first')

            col_names = ['chrom','start','end','name','score',
                         'gene_chrom','gene_start','gene_end',
                         'gene_id','gene_name','gene_strand','distance']
            df_fb = closest.to_dataframe(names=col_names, comment='#')
            df_fb['distance'] = pd.to_numeric(df_fb['distance'], errors='coerce')
            df_fb['assignment_method'] = 'nearest_fallback'
            result_rows.append(df_fb)

        # ------------------------------------------------------------------ #
        # 6. Concatenate and report                                           #
        # ------------------------------------------------------------------ #
        final_df = pd.concat(result_rows, ignore_index=True)

        n_tad      = (final_df['assignment_method'] == 'TAD').sum()
        n_fallback = (final_df['assignment_method'] == 'nearest_fallback').sum()
        n_genes    = final_df['gene_name'].nunique()

        print(f"[{self.name}] Assignment complete: "
              f"{n_tad} via TAD, {n_fallback} via nearest fallback, "
              f"{n_genes} unique genes.")

        return final_df
