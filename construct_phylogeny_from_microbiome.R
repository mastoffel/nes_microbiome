# phyloseq
library(phyloseq)
library(ggplot2)
library(dplyr)
library(stringr)
library(msa)
library(inbreedR)
library(readxl)
library(stringr)
library(dplyr)
library(rptR)
library(plyr)
library(reshape2)
library(lme4)
library(DESeq2)
library(dada2)
library(phangorn)
# todo: filter fecal sample or let it in!

# input folder
input_folder <- "primer_clipped_reads_22_220250_pool"
# load taxa and RSV table
load(paste0("output/", input_folder, "/taxa.RData"))
load(paste0("output/", input_folder, "/seqtab.RData"))
# sample names
samples_out <- rownames(seqtab_nochim)
nes_df <- data.frame("id" = samples_out)

# construct the phylogenetic tree

# multiple alignment of the inferred sequences
seqs <- getSequences(seqtab_nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method = "ClustalW", type = "dna", order = "input")
#load("msa_out.RData")

# The phangorn package constructs a phylogenetic tree
phang_align <- as.phyDat(mult, type = "DNA", names=getSequences(seqtab_nochim))
dm <- dist.ml(phang_align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit <- pml(treeNJ, data = phang_align)

## negative edges length changed to 0!
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

save(fitGTR, file = "../nes_microbiome/output/primer_clipped_reads_22_220250_pool/fitGTR.RData")
# detach("package:phangorn", unload=TRUE)