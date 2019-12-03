# makes a phylogeny from ASVs
# phyloseq
library(phyloseq)
library(msa)
library(DESeq2)
library(dada2)
library(tidyverse)

# todo: filter fecal sample or let it in!

# input folder
input_folder <- "primer_clipped_reads_22_250250_pool"
# load taxa and RSV table
load(paste0("output/", input_folder, "/taxa.RData"))
load(paste0("output/", input_folder, "/seqtab.RData"))
# sample names
samples_out <- rownames(seqtab_nochim)
nes_df <- data.frame("id" = samples_out)

# construct the phylogenetic tree
seqs <- getSequences(seqtab_nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method = "ClustalW", type = "dna", order = "input")

save(mult, file = "msa_out.RData")
