# check pipeline output

# input folder
input_folder <- "primer_clipped_reads_11_220230_pool"
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
taxa_print4 <- as.data.frame(taxa_print)


# 250250
taxa_print1 <- taxa_print
sort(table(taxa_print1$Class), decreasing = TRUE)[1:10]

# 220240
sort(table(taxa_print2$Class), decreasing = TRUE)[1:10]

# 220250
sort(table(taxa_print3$Class), decreasing = TRUE)[1:10]

# 230250
sort(table(taxa_print4$Class), decreasing = TRUE)[1:10]

# taxa_print2 <- taxa_print #[taxa_print[,1] != "Eukaryota", ]
sum(sort(table(taxa_print[,5]), decreasing = TRUE))
sum(sort(table(taxa_print_nopool[,5]), decreasing = TRUE))

taxa_print_nopool <- taxa_print

taxa_print_22 <- taxa_print

sort(table(taxa_print_combined[,4]), decreasing = TRUE)[1:15]
taxa_print_combined <- taxa_print

