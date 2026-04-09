# Repository Overview

This repository contains multiple analysis pipelines used for coaccessibility scoring, model prediction,  sequence design, and downstream visualization. Each pipeline performs a distinct component of the analysis and is associated with a specific Conda environment (see `Requirements` below).

---



### `rank_pipeline`
Performs FIRE peak processing, stitching, and co-accessibility-based ranking. Includes scripts for preprocessing Fiber-seq data and generating ranked regulatory elements used in downstream analyses.

See https://github.com/fiberseq/fibertools-rs and https://github.com/fiberseq/FIRE for further documentation

**Running instructions**  

1. CCS files with HIFI kinetics are run with `run_fiberseq.sh`

2. FIRE peaks are called on fiberseq output with `run_FIRE.sh`

3. FIRE peak co-accessibility is scored and ranked with  `FIRE_to_ranked_coaccessibility.py`
```
usage: FIRE_to_ranked_coaccessibility.py [-h] -p PEAK -c CHROM -a ACC -g GFF
                                         -o OUTPUT

```
---
### `basenji_pipeline`
Handles model-based prediction using Basenji. This includes preparation of genomic tracks (e.g., BigWig files), organization of model inputs, and execution of trained models using provided parameters and weights.

See https://github.com/calico/basenji for further documentation

**Running instructions**  

1. Basenji input tracks are prepared with `make_bigwig_tracks`
   - Basenji input tracks are avaiable at https://zenodo.org/records/19477875
2. Basenji Pipeline is run with `run_basenji.sh`
   - Model weights are available in the subdirectory `basenji_pipeline\model_params`

---
### `INSIDE_pipeline`
Implements the INSIDE framework for sequence design and in silico evolution. This includes sequence generation, genome-context evaluation, scoring of parental sequences, and motif analysis. Parameter configurations are provided as JSON files for reproducible runs.

**Running instructions**  

1. INSIDE run conditions are prepared with `INSIDE_setup.py` and `INSIDE_within_genome_setup.py`
   - INSIDE run conditions are avaiable in the subdirectory `INSIDE_pipeline\json`
2. INSIDE Pipeline is run with `run_INSIDE.sh`


---

### `dafseq_pipeline`
Generates sequence constructs for DAF-seq experiments, including reporter and locus-specific FASTA files (e.g., Cerulean, Meox, Shh). Also includes utilities for summarizing reporter outputs.

See methods in https://www.biorxiv.org/content/10.1101/2024.11.06.622310v1 for further documentation

**Running instructions**  
1. FASTAs are prepared with python scripts in the subdirectory `dafseq_pipeline\make_fastas`
   - Custom FASTAS are avaiable in the subdirectory `INSIDE_pipeline\fastas`
2. dafseq pipeline is run with `dafseq_reporter.py`

---

### `figure_generation`
Generates all figures associated with the study. Scripts are organized such that file names correspond directly to the figure panels they produce (e.g., main figures and supplementary figures). The processed `data` directory is available at https://zenodo.org/records/19477875



---


# Requirements
This repository contains multiple analysis pipelines, each with distinct software and dependency requirements. To ensure reproducibility, Conda environment files (`.yml`) are provided and should be used as specified below.


| Directory              | Dependencies         |
|----------------------|------------------------|
| `INSIDE_pipeline`      | `basenji.yml`          |
| `basenji_pipeline`     | `basenji.yml`          |
| `dafseq_pipeline`      | `dafSeq.yml`           |
| `figure_generation`    | `figure_generation.yml`|
| `rank_pipeline`*       | `PEP_cluster.yml`      |

\* `run_fire.sh` requires `FIRE.yml` instead.

# Data Reproducibility
Processed data for this project is hosted at https://zenodo.org/records/19477875 as well as external public repositories

&emsp;See `data_reproducibility.txt` for download instructions