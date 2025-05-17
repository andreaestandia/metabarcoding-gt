---
  title: "Step 2: Primer Removal with Cutadapt"
author: "Your Name"
date: "`r Sys.Date()`"
format: html
---
  
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
source("common_functions.R")
install_required_packages()
```

## Step 2: Primer Removal with Cutadapt

This step removes PCR primers from the reads using Cutadapt. Primer sequences can heavily influence downstream analyses, so it's important to remove them before proceeding.

### Load Data from Previous Step

```{r load-data}
# Load locus information
LOCUS <- readRDS(file.path(RESULTS_DIR, "locus.rds"))
LOCUS_DIR <- file.path(RESULTS_DIR, LOCUS)

# Load filtered file paths
fnFs.filtN <- readRDS(file.path(LOCUS_DIR, "R_objects", "fnFs.filtN.rds"))
fnRs.filtN <- readRDS(file.path(LOCUS_DIR, "R_objects", "fnRs.filtN.rds"))

# Check which files exist (some may have been filtered out in previous step)
fnFs.filtN.exists <- file.exists(fnFs.filtN)
fnRs.filtN.exists <- file.exists(fnRs.filtN)

# Filter to only existing files
fnFs.filtN <- fnFs.filtN[fnFs.filtN.exists]
fnRs.filtN <- fnRs.filtN[fnRs.filtN.exists]

cat("Number of samples after N filtering:", length(fnFs.filtN), "\n")
```

### Set Primer Sequences

```{r primers}
# Get default primers for the selected locus
primers <- get_default_primers(LOCUS)
FWD <- primers$FWD
REV <- primers$REV

# Optional: Override with custom primers if needed
# FWD <- "CUSTOM_FORWARD_PRIMER"
# REV <- "CUSTOM_REVERSE_PRIMER"

cat("Forward primer:", FWD, "\n")
cat("Reverse primer:", REV, "\n")

# Generate all possible orientations of primers
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Check primer presence in first sample
if (length(fnFs.filtN) > 0) {
  # Count primer occurrences in the first sample
  pre_trim_primer_counts <- rbind(
    FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]])
  )
  
  # Save primer counts
  write.table(pre_trim_primer_counts, 
              file.path(LOCUS_DIR, "outputs", "pre_trim_primer_counts.tsv"), 
              col.names=NA, sep="\t")
  
  # Display counts
  print(pre_trim_primer_counts)
} else {
  cat("No samples available for primer checking\n")
}
```

### Prepare Cutadapt Command

```{r cutadapt-prep}
# Set up cutadapt output files
fnFs.cut <- file.path(LOCUS_DIR, "cutadapt", basename(fnFs.filtN))
fnRs.cut <- file.path(LOCUS_DIR, "cutadapt", basename(fnRs.filtN))

# Save paths for later steps
saveRDS(fnFs.cut, file.path(LOCUS_DIR, "R_objects", "fnFs.cut.rds"))
saveRDS(fnRs.cut, file.path(LOCUS_DIR, "R_objects", "fnRs.cut.rds"))

# Create reverse complements for primer removal
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Set up cutadapt flags
# -g: 5' forward primer
# -a: 3' reverse complement of reverse primer
R1.flags <- paste("-g", FWD, "-a", REV.RC) 

# -G: 5' reverse primer
# -A: 3' reverse complement of forward primer
R2.flags <- paste("-G", REV, "-A", FWD.RC)

cat("Cutadapt R1 flags:", R1.flags, "\n")
cat("Cutadapt R2 flags:", R2.flags, "\n")
```

### Run Cutadapt

```{r run-cutadapt}
# Check if cutadapt is installed
cutadapt_installed <- Sys.which("cutadapt") != ""
if (!cutadapt_installed) {
  cat("WARNING: cutadapt not found in PATH. Please install cutadapt:\n")
  cat("pip install cutadapt\n")
  cat("Then restart the analysis.\n")
} else {
  # Set additional parameters
  min_length <- 50  # Minimum length after trimming
  n_copies <- 2     # Number of copies of primers to search for
  
  # Run cutadapt on each sample
  for (i in seq_along(fnFs.filtN)) {
    # Construct and run cutadapt command
    cmd <- paste("cutadapt", 
                 R1.flags, R2.flags,
                 "-o", fnFs.cut[i], 
                 "-p", fnRs.cut[i],
                 fnFs.filtN[i], fnRs.filtN[i],
                 "-j 0",  # Use all available cores
                 "--discard-untrimmed",  # Discard reads without primers
                 "-m", min_length,  # Minimum length filter
                 "-n", n_copies)  # Number of times to search for adapters
    
    cat("Processing sample", i, "of", length(fnFs.filtN), ":", basename(fnFs.filtN[i]), "\n")
    system(cmd)
  }
  
  # Verify the output files exist
  fnFs.cut.exists <- file.exists(fnFs.cut)
  fnRs.cut.exists <- file.exists(fnRs.cut)
  
  cat("Samples with cutadapt output:", sum(fnFs.cut.exists), "\n")
  
  # Check primer removal in first sample
  if (sum(fnFs.cut.exists) > 0) {
    post_trim_primer_counts <- rbind(
      FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[fnFs.cut.exists][1]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[fnRs.cut.exists][1]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[fnFs.cut.exists][1]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[fnRs.cut.exists][1])
    )
    
    write.table(post_trim_primer_counts, 
                file.path(LOCUS_DIR, "outputs", "post_trim_primer_counts.tsv"), 
                col.names=NA, sep="\t")
    
    print(post_trim_primer_counts)
  } else {
    cat("No samples passed primer trimming\n")
  }
}
```

### Quality Plots After Trimming

```{r post-trim-quality}
# Generate quality plots for trimmed reads
if (sum(fnFs.cut.exists) > 0) {
  # Plot quality profiles for the first few samples
  n_to_plot <- min(3, sum(fnFs.cut.exists))
  
  # Forward reads quality plot
  fwd_plot <- plotQualityProfile(fnFs.cut[fnFs.cut.exists][1:n_to_plot])
  save_quality_plot(fwd_plot, file.path(LOCUS_DIR, "plots", "post_trim_forward_quality_plots.pdf"))
  
  # Reverse reads quality plot
  rev_plot <- plotQualityProfile(fnRs.cut[fnRs.cut.exists][1:n_to_plot])
  save_quality_plot(rev_plot, file.path(LOCUS_DIR, "plots", "post_trim_reverse_quality_plots.pdf"))
  
  # Display plots in the report
  fwd_plot
  rev_plot
  
  # Update sample names based on what passed
  if (sum(fnFs.cut.exists) > 0) {
    sample.names <- sapply(fnFs.cut[fnFs.cut.exists], get_sample_name)
    saveRDS(sample.names, file.path(LOCUS_DIR, "R_objects", "cutadapt_sample_names.rds"))
    cat("Sample names saved for", length(sample.names), "samples after cutadapt\n")
  }
}
```

## Next Steps

After primer removal, you can now proceed to the next step, which is quality filtering and trimming. The quality plots after primer removal will help determine where quality drops and inform the truncation parameters for the next step.