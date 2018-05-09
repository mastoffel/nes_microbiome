# script for storing unused code


# gm_mean = function(x, na.rm=TRUE){
#     exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# }
# geoMeans <- apply(counts(GPdds), 1, gm_mean)
# GPdds <- estimateSizeFactors(GPdds, geoMeans=geoMeans)


# core microbiome
# nes_core <- core(ps2, detection = .2/100, prevalence = 80/100)
# plot_core(transform(nes_core, "compositional"),
#             plot.type = "heatmap",
#             colour = gray(seq(0,1,length=5)),
#             prevalences = seq(.05, 1, 05),
#             detections = 10^seq(log10(1e-3), log10(.2), length = 10), 
#             horizontal = TRUE) +
#             xlab("Detection Threshold (Relative Abundance (%))") +
#             theme(axis.text.x = element_blank())



# FILTERS

# Retain taxa with abundance > .1%
# ps_rel <- transform_sample_counts(ps2, function(x) x / sum(x))
# ps_filt <- filter_taxa(ps_rel, function(x) any(x > 0.001),TRUE)
# ps3 <- prune_taxa(taxa_names(ps_filt), ps2)

# alternative filter
# Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
#ps3 <- filter_taxa(ps2, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
#ps3 <- filter_taxa(ps2, function(x) sd(x)/mean(x) > 3.0, TRUE)
# sort(table(tax_table(ps3)[, "Phylum"], exclude = NULL))


# TRANSFORMATIONS

# (1) Simple transformation -----------------------------------
#ps_rel <- transform_sample_counts(ps3, function(x){x / sum(x)})
#ps_rel <- transform_sample_counts(ps3, function(x) log(x + 1))
#ps_rel <- transform_sample_counts(ps3, function(x) log(x + 1) / sum(log(x+1)))

# (2) CLR Centered Log-Ratio transformation -------------------
ps_clr <- microbiome::transform(ps3, 'clr')
# has to be modiefied to avoid negative values for ordination
ps_clr_mod <- ps_clr
otu_table(ps_clr_mod)[otu_table(ps_clr_mod) < 0] <- 0

#define a zero- and NA-tolerant function for calculating CLR
# from https://github.com/joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
# gm_mean = function(x, na.rm=TRUE){
#     # The geometric mean, with some error-protection bits.
#     exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
# }
# clr = function(x, base=2){
#     x <- log((x / gm_mean(x)), base)
#     x[!is.finite(x) | is.na(x)] <- 0.0
#     return(x)
# }
# ps_rel <- transform_sample_counts(ps3, clr)
