# from raw reads to ASV table
library(dplyr)
library(dada2)
library(magrittr)
help(package="dada2")
?derepFastq
?dada

# output folder
output_folder <- "combined_reads_1_nopool"
dir.create(paste0(paste0("output/", output_folder)))
# filtering
ees <- 1
# trimming

# path to combined reads
path <- "../data/combined_samples"
list.files(path)

# sort
combn_reads <- sort(list.files(path, pattern="fastq", full.names = TRUE))
plotQualityProfile(combn_reads[1:2])

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
temp_names <- sapply(strsplit(basename(combn_reads), "-"), `[`, 7) 
sample_names <- sapply(strsplit(temp_names, "_"), `[`, 1) 

# assigning file names
filt_path <- file.path(path, "filtered")
filtCs <- file.path(filt_path, paste0(sample_names, "_C_filt.fastq.gz"))

# no trimming here, just filtering
out <- filterAndTrim(combn_reads, filtCs, 
                     maxN = 0, maxEE = 1, truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# error rates
errC <- learnErrors(filtCs, multithread = TRUE)
plotErrors(errC, nominalQ = TRUE)

# dereplication
derepCs <- derepFastq(filtCs, verbose = TRUE)
names(derepCs) <- sample_names 

# sample inference
dadaCs <- dada(derepCs, err = errC, multithread = TRUE)
dadaCs[[1]]

# construct sequence table
seqtab <- makeSequenceTable(dadaCs)
dim(seqtab)

# inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# remove weird sequence lengths
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(400,443)]
table(nchar(getSequences(seqtab2)))

# remove chimaras
seqtab_nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab_nochim)
# fraction of chimeras (total reads)
sum(seqtab_nochim)/sum(seqtab2)

# track reads through the pipeline
# no step should result in the loss of too many reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaCs, getN),
               rowSums(seqtab2), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoised","tabled", "nonchim")
rownames(track) <- sample_names
head(track)
save(track, file = paste0("output/", output_folder, "/track.RData"))

# Assign Taxonomy
taxa <- assignTaxonomy(seqtab_nochim, paste0(path, "/silva_nr_v128_train_set.fa.gz"), multithread=TRUE)

# add species
taxa <- addSpecies(taxa, paste0(path, "/silva_species_assignment_v128.fa.gz"))

save(seqtab_nochim, file= paste0("output/", output_folder, "/seqtab.RData"))
save(taxa, file=paste0("output/", output_folder, "/taxa.RData"))

taxa_print <- taxa # Removing sequence rownames for display only
rownames(taxa_print) <- NULL
head(taxa_print)
table(taxa_print[, 6])

