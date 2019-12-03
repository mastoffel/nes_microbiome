# DADA2 pipeline to go from paired-end fastq files 
# to an amplicon sequence variant (ASV) table
# pipeline starts from reads where the primers and adapters have been removed
# from the raw reads

library(dplyr)
library(dada2)
library(magrittr)

# output folder name
output_folder <- "primer_clipped_reads_22_220230_pool"
dir.create(paste0(paste0("output/", output_folder)))

# filtering
ees <- c(2,2)

# trimming
totrim <- c(220,230) #220,230

# path to (nearly) raw reads
# primers and adapters have been clipped, see methods in the paper
path <- "../raw_reads/primer_clipped_reads"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

# sample names
temp_names <- sapply(strsplit(basename(fnFs), "-"), `[`, 7) 
sample_names <- sapply(strsplit(temp_names, "_"), `[`, 1) 

# quality
plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnRs[1:12])

# filtering and trimming
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample_names, "_R_filt.fastq.gz"))

# filter and save filtered files
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = totrim,
                     maxN = 0, maxEE = ees, truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
head(out)

# learn the error rates (computationally expensive)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# check errors, points should be along the line
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# dereplication
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool=TRUE) # pool=TRUE
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool=TRUE) # pool = TRUE

# save.image(file = "dada_adapter_clipped_nopool.RData")
dadaFs[[1]]

# merge the denoised forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
mergers[[1]]

class(mergers)

# sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# inspect distribution of sequence length
table(nchar(getSequences(seqtab)))

# cut sequences that are too short
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(380,450)]

# remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab_nochim)

# how much of the total reads remain?
sum(seqtab_nochim)/sum(seqtab)

# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample_names
head(track)
save(track, file = paste0("output/", output_folder, "/track.RData"))

# Assign Taxonomy
# silva database needs to be downloaded and put into the folder specified in path
taxa <- assignTaxonomy(seqtab_nochim, paste0(path, "/silva_nr_v128_train_set.fa.gz"), multithread=TRUE)

# add species level assignments where possible
taxa <- addSpecies(taxa, paste0(path, "/silva_species_assignment_v128.fa.gz"))

# save ASV table
save(seqtab_nochim, file=paste0("output/", output_folder, "/seqtab.RData"))
save(taxa, file=paste0("output/", output_folder, "/taxa.RData"))








