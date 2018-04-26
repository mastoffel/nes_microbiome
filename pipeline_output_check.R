# check pipeline output

# input folder
input_folder <- "primer_clipped_reads_22_220250_pool"
# load taxa and RSV table
load(paste0("output/", input_folder, "/taxa.RData"))
load(paste0("output/", input_folder, "/seqtab.RData"))
load(paste0("output/", input_folder, "/track.RData"))

head(track)

# check taxa
taxa_print <- taxa # Removing sequence rownames for display only
rownames(taxa_print) <- NULL
head(taxa_print)

# check taxa
taxa_print <- as.data.frame(taxa_print)
table(taxa_print$Kingdom)
# taxa_print2 <- taxa_print #[taxa_print[,1] != "Eukaryota", ]
sum(sort(table(taxa_print[,5]), decreasing = TRUE))
sum(sort(table(taxa_print_nopool[,5]), decreasing = TRUE))

taxa_print_nopool <- taxa_print

taxa_print_22 <- taxa_print

sort(table(taxa_print_combined[,4]), decreasing = TRUE)[1:15]
taxa_print_combined <- taxa_print

