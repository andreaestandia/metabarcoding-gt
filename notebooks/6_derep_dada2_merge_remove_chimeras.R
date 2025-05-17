---
  title: "Step 6: 06_derep_dada2_merge_remove_chimeras.R"
author: "Your Name"
date: "`r Sys.Date()`"
format: html
editor_options: 
  chunk_output_type: console
---
  
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
source(here::here("src", "metabarcoding-gt.R"))
install_required_packages()
```

## read in the filt files
filtFs <- readRDS(file = paste(path, "/R_objects/04_fnFs.filtN.rds", sep=""))
filtRs <- readRDS(file = paste(path, "/R_objects/04_fnRs.filtN.rds", sep=""))

## read in the error models
errF <- readRDS(file = paste(path, "/R_objects/05_errF.rds", sep=""))
errR <- readRDS(file = paste(path, "/R_objects/05_errR.rds", sep=""))

## check files exist
exists <- file.exists(filtFs)

## dereplicate the reads
derepFs <- derepFastq(filtFs[exists])
derepRs <- derepFastq(filtRs[exists])

# get sample names
sample.names <- readRDS(file = paste(path, "/R_objects/03_sample_names.rds", sep=""))

## name the derep class objects
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

## perform dada2 sample inference
dadaFs <- dada(derepFs, err = errF, multithread=TRUE)
dadaRs <- dada(derepRs, err = errR, multithread=TRUE)

## merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

## make ASV table
seqtab <- makeSequenceTable(mergers)

## remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

## write out ASV sequences
sequences <- colnames(seqtab.nochim)
headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  headers[i] <- paste(">ASV", i, sep="_")
}

output_fasta <- c(rbind(headers, sequences))

if (!is.null(opt$marker)){
  write(output_fasta, paste(path, "/working_data/06_", opt$marker, "ASV_seqs.fasta", sep=""))
}
if (is.null(opt$marker)){
  write(output_fasta, paste(path, "/working_data/06_ASV_seqs.fasta", sep=""))
}

# write out ASV counts
tab <- t(seqtab.nochim)

write.table(tab, file = paste(path, "/working_data/06_", opt$marker, "ASV_counts.tsv", sep=""), sep="\t", quote=F, col.names=NA)

saveRDS(dadaFs, file = paste(path, "/R_objects/06_dadaFs.rds", sep=""))
saveRDS(dadaRs, file = paste(path, "/R_objects/06_dadaRs.rds", sep=""))
saveRDS(mergers, file = paste(path, "/R_objects/06_mergers.rds", sep=""))
saveRDS(seqtab.nochim, file = paste(path, "/R_objects/06_seqtab.nochim.rds", sep=""))