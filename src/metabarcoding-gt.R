# Source file containing common functions, paths, and libraries
# This file should be included at the beginning of each QMD script

# Define project paths
PROJECT_ROOT <- "/media/estandia/LaCie/projects/metabarcoding-gt"
DATA_DIR <- file.path(PROJECT_ROOT, "data")
RAW_DIR <- file.path(DATA_DIR, "raw")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
FIGURES_DIR <- file.path(PROJECT_ROOT, "figures")
TAXONOMY_DB_DIR <- file.path(PROJECT_ROOT, "taxonomy_db")
TEMP_DIR <- file.path(PROJECT_ROOT, "temp")

# Install and load required packages
install_required_packages <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  required_packages <- c(
    # Core packages
    "here", "dada2", "ShortRead", "Biostrings", "ggplot2", "optparse", 
    "gridExtra", "reshape2", "knitr", "rmarkdown", "dplyr",
    # Additional useful packages
    "vegan", "readr", "tidyr", "stringr", "rentrez", "treemapify",
    "patchwork", "purrr"
  )
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% c("dada2", "ShortRead", "Biostrings", "phyloseq")) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
    library(pkg, character.only = TRUE)
  }
}

create_directories <- function(locus_name) {
  # Create directories for storing intermediate and final results
  locus_dir <- file.path(RESULTS_DIR, locus_name)
  
  dirs <- c(
    locus_dir,
    file.path(locus_dir, "filtered_N"),
    file.path(locus_dir, "cutadapt"),
    file.path(locus_dir, "filtered_trimmed"),
    file.path(locus_dir, "outputs"),
    file.path(locus_dir, "plots"),
    file.path(locus_dir, "R_objects")
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  }
  
  return(locus_dir)
}

# Updated get_sample_name function for your specific file structure
get_sample_name <- function(filename) {
  # Extract sample name from the file path
  # Example filename: "/path/Sample_IDT10_UDI_Adapter3553-Irene_A1_COI/IDT10_UDI_Adapter3553-Irene_A1_COI_CGAGAAGATA-CGTCGCCTAT_L001_R1.fastq.gz"
  
  # Get the basename of the file
  basename_file <- basename(filename)
  
  # Extract the sample ID - we'll use the IDT10 part as the sample ID
  # This extracts "IDT10" from the beginning of the filename
  if (grepl("^IDT10_", basename_file)) {
    # Extract the main sample ID
    sample_id <- "IDT10"
    
    # Get the additional sample identifier (e.g., A1, B5, etc.)
    # Looking for pattern like "Irene_A1_COI" and extracting "A1"
    if (grepl("Irene_([A-Z][0-9]+)_", basename_file)) {
      position <- regmatches(basename_file, regexpr("Irene_([A-Z][0-9]+)_", basename_file))
      position <- gsub("Irene_", "", position)
      position <- gsub("_.*$", "", position)
      
      # Combine sample ID with position
      sample_name <- paste0(sample_id, "_", position)
    } else {
      # Fallback if the pattern doesn't match
      sample_name <- sample_id
    }
    
    return(sample_name)
  }
  
  # Fallback: try to extract from parent directory
  parent_dir <- basename(dirname(filename))
  if (grepl("^Sample_IDT10_", parent_dir)) {
    # Extract the position part (A1, B5, etc.)
    if (grepl("Irene_([A-Z][0-9]+)_", parent_dir)) {
      position <- regmatches(parent_dir, regexpr("Irene_([A-Z][0-9]+)_", parent_dir))
      position <- gsub("Irene_", "", position)
      position <- gsub("_.*$", "", position)
      
      # Create sample name
      sample_name <- paste0("IDT10_", position)
      return(sample_name)
    }
  }
  
  # Final fallback: just use the filename without extension and R1/R2 part
  sample_name <- gsub("_L001_R[12]\\.fastq\\.gz$", "", basename_file)
  return(sample_name)
}

# Function to get default primers by locus
get_default_primers <- function(locus) {
  locus <- tolower(locus)  # Convert to lowercase for consistency
  
  if (locus == "16s") {
    return(list(FWD = "GACGAKAAGACCCTA",     # 515F
                REV = "TCTTAATCCAACATCGAGGTC"))   # 806R
  } else if (locus == "rbcl") {
    return(list(FWD = "TTGGCAGCATTYCGAGTAACTCC",  # rbcLa-F
                REV = "AACCYTCTTCAAAAAGGTC"))        # rbcLa-R
  } else if (locus == "coi") {
    return(list(FWD = "ATTCAACCAATCATAAAGATATTGG",     # 515F
                REV = "WACTAATCAATTWCCAAAHCCHCC"))    # jgHCO2198
  } else {
    stop("Unknown locus: ", locus, ". Please provide primer sequences manually.")
  }
}

# Function to save plots
save_quality_plot <- function(plot_data, file_path, width=10, height=8) {
  pdf(file_path, width=width, height=height)
  print(plot_data)
  dev.off()
}

get.sample.name <- function(fname) {
  # Get parts
  parts <- strsplit(basename(fname), "_")[[1]]
  
  # parts look like:
  # [1] "IDT10"       "UDI"     "Adapter3553-Irene" "A1" "COI" "CGAGAAGATA-CGTCGCCTAT" "L001" "R1.fastq.gz"
  
  adapter_part <- parts[3] # e.g. "Adapter3553-Irene"
  well         <- parts[4] # e.g. "A1"
  locus        <- parts[5] # e.g. "COI"
  
  # Extract adapter number
  adapter_num <- as.numeric(gsub("Adapter", "", strsplit(adapter_part, "-")[[1]][1]))
  
  # Determine plate
  plate <- if (!is.na(adapter_num)) {
    if (adapter_num >= 577 & adapter_num <= 671) {
      "P4"
    } else if (adapter_num >= 3553 & adapter_num <= 3648) {
      "P1"
    } else if (adapter_num >= 3649 & adapter_num <= 3744) {
      "P2"
    } else if (adapter_num >= 3745 & adapter_num <= 3839) {
      "P3"
    } else {
      paste0("Adapter", adapter_num)
    }
  } else {
    "UnknownPlate"
  }
  
  paste(well, plate, locus, sep = "_")
}

get.sample.name.16s <- function(fname) {
  parts <- strsplit(basename(fname), "_")[[1]]
  
  # parts:
  # [1] "IDT10" "UDI" "Adapter668-F19" "AAGTTAGGAC-TCACAGCTGC" "L001" "R1.fastq.gz"
  
  adapter_part <- parts[3]
  adapter_split <- strsplit(adapter_part, "-")[[1]]
  
  adapter <- adapter_split[1]
  well    <- adapter_split[2]
  
  adapter_num <- as.numeric(gsub("Adapter", "", adapter))
  
  plate <- if (!is.na(adapter_num)) {
    if (adapter_num >= 579 & adapter_num <= 671) {
      "P4"
    } else if (adapter_num >= 3553 & adapter_num <= 3648) {
      "P1"
    } else if (adapter_num >= 3649 & adapter_num <= 3744) {
      "P2"
    } else if (adapter_num >= 3745 & adapter_num <= 3839) {
      "P3"
    } else {
      paste0("Adapter", adapter_num)
    }
  } else {
    "UnknownPlate"
  }
  
  # Add "_16s" at the end
  paste(well, plate, "16s", sep = "_")
}




# Function to get the number of reads at each step
getN <- function(x) sum(getUniques(x))

# Function to create all primer orientations
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Function to count primer occurrences
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

library(showtext)
library(sysfonts)
library(ggplot2)

# Add Roboto Condensed from Google Fonts
font_add_google("Roboto Condensed", "roboto_condensed")

# Automatically use showtext for new plots
showtext_auto()
