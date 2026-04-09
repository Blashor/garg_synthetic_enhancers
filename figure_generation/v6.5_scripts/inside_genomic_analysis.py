#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparative analysis of intergenic co-accessible enhancers across cell types.

Produces:
  - Pairwise Venn diagrams comparing F123g intergenic enhancers to V6.5,
    GM12878, and K562
  - Pairwise scatter plots of enhancer odds ratios between cell types
  - Ranked co-accessibility hockey stick plot for V6.5 with Sox2 and
    Mir290a loci annotated using a TAD-constrained gene assignment strategy

Input files required per sample:
  - ce_rank.txt            : tab-separated co-accessibility rankings (no header)
  - FIRE_peaks_intergenic.bed : intergenic FIRE peak calls

Additional inputs:
  - mm10 GTF annotation (GENCODE vM23 basic)
  - mESC TADs in mm10 coordinates (Dixon et al. 2012, lifted from mm9)
  - liftOver chain file (mm10 -> hg38) for cross-species comparisons

Usage:
  Edit the USER CONFIGURATION section below to point to your data files,
  then run:
      python inside_genomic_analysis.py

Dependencies: pybedtools, pandas, scipy, matplotlib, seaborn
See requirements.txt for full version information.

@author: michelemeline
"""

import sys
from pathlib import Path

# -----------------------------------------------------------------------
# USER CONFIGURATION — edit these paths before running
# -----------------------------------------------------------------------

# Directory for output figures (created automatically if absent)
OUTPUT_DIR = Path('results')

# Sequencing data root — update to your local data directory
SEQ_ROOT = Path('/path/to/sequencing/data')

# Sample-specific file paths, colors, and genome build flags
SAMPLE_PATHS = {
    'V6.5': {
        'ce_rank':         SEQ_ROOT / 'mESC_reporter/data/ce_rank.txt',
        'fire_intergenic': SEQ_ROOT / 'mESC_reporter/scored/FIRE_peaks_intergenic.bed',
        'color':           'blue',
        'is_hg38':         False,
    },
    'F123g': {
        'ce_rank':         SEQ_ROOT / 'INSIDE/FDR5pct/data/ce_rank.txt',
        'fire_intergenic': SEQ_ROOT / 'INSIDE/FDR5pct/scored/FIRE_peaks_intergenic.bed',
        'color':           '#C0382F',
        'is_hg38':         False,
    },
    'GM12878': {
        'ce_rank':         SEQ_ROOT / 'GM12878/data/ce_rank.txt',
        'fire_intergenic': SEQ_ROOT / 'GM12878/scored/FIRE_peaks_intergenic.bed',
        'color':           'g',
        'is_hg38':         True,
    },
    'K562': {
        'ce_rank':         SEQ_ROOT / 'K562/data/ce_rank.txt',
        'fire_intergenic': SEQ_ROOT / 'K562/scored/FIRE_peaks_intergenic.bed',
        'color':           '#D4AC0D',
        'is_hg38':         True,
    },
}

# liftOver resources (mm10 -> hg38) — used for V6.5 and F123g cross-species comparisons
MM10_TO_HG38_CHAIN = Path('~/ucsc/mm10ToHg38.over.chain.gz')
LIFTOVER_EXE       = Path('/usr/local/bin/liftOver')

# Annotation files
GTF_PATH = Path('/path/to/gencode.vM23.basic.annotation.gtf')
TAD_PATH = Path('/path/to/mESC_TADs_mm10.bed')  # Dixon et al. 2012, lifted mm9->mm10

# -----------------------------------------------------------------------
# END USER CONFIGURATION
# -----------------------------------------------------------------------

sys.path.insert(0, str(Path(__file__).parent))

from genomic_sample import GenomicSample
from genomic_analysis import GenomicAnalysis

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def build_samples():
    """Instantiate GenomicSample objects from SAMPLE_PATHS config."""
    samples = {}
    for name, cfg in SAMPLE_PATHS.items():
        kwargs = dict(
            name                       = name,
            ce_rank_path               = str(cfg['ce_rank']),
            fire_peaks_intergenic_path = str(cfg['fire_intergenic']),
            color                      = cfg['color'],
            is_hg38                    = cfg.get('is_hg38', False),
        )
        if not cfg.get('is_hg38', False):
            kwargs['chain_file']    = str(MM10_TO_HG38_CHAIN)
            kwargs['liftover_path'] = str(LIFTOVER_EXE)
        samples[name] = GenomicSample(**kwargs)
    return samples


def main():
    samples  = build_samples()
    analysis = GenomicAnalysis(samples)

    # ---- Venn diagrams -------------------------------------------------
    analysis.create_venn2(
        'F123g', 'V6.5',
        feature_type = 'intergenic_enhancers',
        output_file  = str(OUTPUT_DIR / 'V6.5_venn.svg'),
    )
    analysis.create_venn2(
        'F123g', 'GM12878',
        feature_type = 'hg38_intergenic_enhancers',
        output_file  = str(OUTPUT_DIR / 'GM12878_venn.svg'),
    )
    analysis.create_venn2(
        'F123g', 'K562',
        feature_type = 'hg38_intergenic_enhancers',
        output_file  = str(OUTPUT_DIR / 'K562_venn.svg'),
    )

    # ---- Scatter plots -------------------------------------------------
    analysis.create_scatter_comparison(
        'F123g', 'V6.5',
        feature_type = 'intergenic_enhancers',
    )
    analysis.create_scatter_comparison(
        'F123g', 'GM12878',
        feature_type = 'hg38_intergenic_enhancers',
        output_file  = str(OUTPUT_DIR / 'GM12878_scatter.svg'),
    )
    analysis.create_scatter_comparison(
        'F123g', 'K562',
        feature_type = 'hg38_intergenic_enhancers',
        output_file  = str(OUTPUT_DIR / 'K562_scatter.svg'),
    )

    # ---- Hockey stick plot with TAD-based gene annotation --------------
    df, region_ids = analysis.plot_hockey_stick(
        sample_name  = 'V6.5',
        gtf_path     = str(GTF_PATH),
        tad_path     = str(TAD_PATH),
        gene_names   = ['Sox2', 'Mir290a'],
        feature_type = 'intergenic_enhancers',
        top_only     = True,
        output_file  = str(OUTPUT_DIR / 'V6.5_hockey_stick.svg'),
    )


if __name__ == '__main__':
    main()
