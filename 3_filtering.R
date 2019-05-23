# Filter ASV table

# Filtering --------------------------------------------------------------------

# load phyloseq object
load("../data/processed/ps0.RData")

# which samples have low abundance?
#sort(sample_sums(ps0))
#ntaxa(ps0) 

# Remove ASVs that are Chloroplasts, Mitochondria, or have NA on Class level (77 overall)
ps <- ps0 %>%
    subset_taxa(
        ((Family != "Mitochondria") | is.na(Family)) & (Class  != "Chloroplast")
    )
#ntaxa(ps)

# Prevalence of each ASV
# (in how many samples did each ASV appear at least once)
prev0 <- apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts 
prevdf <- data.frame(Prevalence = prev0,
                     TotalAbundance = taxa_sums(ps),
                     tax_table(ps))

# Kingdoms: 2791 Bacteria, 12 Eukaryota, 6 Archaea
#sort(table(prevdf$Family, useNA = "always"), decreasing = TRUE)
#names(sort(table(prevdf$Family, useNA = "always"), decreasing = TRUE))

# checks
prevdf_check <- as_tibble(prevdf)
hist(prevdf_check$Prevalence, breaks = 100)

# check how many phyla and whether there are NAs
sort(table(prevdf$Class, useNA = "always"), decreasing = TRUE) 

#### prepare prevalence filtering
# Spirochaetae appear and three samples and are pathogenic
sort(table(prevdf$Phylum))
keepPhyla <- table(prevdf$Phylum)[(table(prevdf$Phylum) > 2)] # phylum appears minimum 3 samples
prevdf1 <- subset(prevdf, Phylum %in% names(keepPhyla))

# Keep taxa when appearing in minimum 2% samples (3 samples)
prevalenceThreshold <- 0.02 * nsamples(ps)
prevalenceThreshold

# Filter 1: execute prevalence filter minimum 2% samples (3 samples)
ntaxa(ps)
ps1 <- prune_taxa((prev0 > prevalenceThreshold), ps)
ntaxa(ps1)

# FILTER (not used)
# ps2 <- subset_taxa(ps1, Phylum %in% names(keepPhyla))
# ps2
# sort(table(tax_table(ps2)[, "Phylum"], exclude = NULL))

# Filter 2: taxa have to have at least 30 reads
#plot(sort(taxa_sums(ps), TRUE), type="h", ylim=c(0, 100))
min_reads <- 30
ps3 <- filter_taxa(ps1, function(x) sum(x) > min_reads, TRUE)
ntaxa(ps3)

# prevalence-abundance plot to get insights into how to filter
p_prev <- ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = Phylum)) +
    geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
    geom_vline(xintercept = min_reads, alpha = 0.5, linetype = 2) +
    geom_point(size = 2, alpha = 0.6) +
    scale_y_log10() + 
    scale_x_log10(labels = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)), 
                  breaks = c(10, 100, 1000, 10000, 100000)) +
    xlab("Total Abundance") +
    facet_wrap(~Phylum) +
    #ggtitle("Abundance by Phylum") +
    theme_minimal() +
    guides(color=FALSE)
p_prev

#ggsave("../figures/Sup1_PrevVsAbund.jpg", p_prev, width = 7, height = 6)

# save full phyloseq object
save(ps3, file = "../data/processed/ps3.RData")

