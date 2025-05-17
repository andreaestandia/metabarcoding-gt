# Metabarcoding of great tit chick faecal samples – Wytham Woods

This repository contains the pipeline used for metabarcoding analyses of great tit chick faecal samples collected in **Wytham Woods**. The goal is to characterise diet composition using high-throughput sequencing data and bioinformatic processing via the DADA2 workflow.

---

## 🔬 Pipeline Overview

The workflow is structured into **7 main steps**, implemented using a mix of `.R` and `.qmd` (Quarto) files:

1. **01_remove_Ns.qmd**  
   Removes reads containing ambiguous bases (Ns) from the raw FASTQ files.

2. **2_cutadapt.R / 2_cutadapt.qmd**  
   Trims primers and adapters from sequences using `cutadapt`.

3. **3_raw_quality_plots.qmd**  
   Generates quality profiles of raw reads to assess sequencing quality.

4. **4_trim_filter.R / 4_trim_filter.qmd**  
   Filters and trims reads based on quality thresholds to prepare for downstream processing.

5. **5_generate_model_error.R / 5_generate_model_error.qmd**  
   Learns and models sequencing error rates, a key step in the DADA2 pipeline.

6. **6_derep_dada2_merge_remove_chimeras.R / .qmd**  
   Performs dereplication, sequence variant inference, merging of read pairs, and removal of chimeras.

7. **7_sequence_tracking.qmd**  
   Tracks the number of sequences remaining at each processing step for quality control.

8. **8_assign_taxonomy.qmd**  
   Assigns taxonomy to amplicon sequence variants (ASVs) using a reference database.

9. **9_stats.qmd**  
   Summarizes the dataset, e.g., number of ASVs per sample, taxonomic composition.

---

## 📁 Structure

Each step of the pipeline is organized as an individual script or Quarto document for transparency and reproducibility.

---

## 👥 Contributors

- Andrea Estandía  (Postdoctoral Researcher, University of Oxford)
- Irene Martínez Baquero (PhD student, University of Oxford)

---

