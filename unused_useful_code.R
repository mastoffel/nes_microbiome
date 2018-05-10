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


# rlog
# class(GPdds)
# r_log_trans <- rlog(GPdds, blind = FALSE)
# #r_log_trans2 <- rlog(GPdds, blind = TRUE, fast = TRUE)
# rlogMat <- assay(r_log_trans)
# rlogMat[rlogMat<0]<-0
# rownames(rlogMat) <- taxa_names(ps3)
# ps_rlog <- ps3
# otu_table(ps_rlog) <- otu_table(rlogMat, taxa_are_rows = TRUE)
# # plot ordination
# # calculate ordination
# ps_ord <- ordinate(ps_rlog, "MDS", "bray")
# # calculate axis length relationships according to eigenvalues
# evals <- ps_ord$values$Eigenvalues
# # get df
# plot_ordination(ps_rlog, ps_ord, shape = "sex", color = "timepoint")
# p_ord_df <- plot_ordination(ps_rlog, ps_ord, shape = "sex", color = "timepoint", justDF = TRUE)
# 
# p_ord_plot <- ggplot(p_ord_df, aes(Axis.1, Axis.2)) +
#     geom_point(size = 3, alpha = 0.8, aes(shape = sex, color = timepoint)) +
#     #geom_point(size = 3, alpha = 0.8, aes( color = individual)) +
#     theme_martin() +
#     theme(panel.grid = element_blank()) +
#     #scale_shape_manual(values = c(21,2))+
#     scale_color_manual(values = colpal) +
#     coord_fixed(sqrt(evals[2] / evals[1])) +
#     theme(legend.position = "bottom",
#         legend.direction = "horizontal") +
#     xlab("Axis 1 [20,2%]") +
#     ylab("Axis 2 [11,3%]")
# p_ord_plot 




# DESEQ2 ANALYSIS

# calculate deseq output and put in dataframe ======================================
calc_deseq_table_timepoint <- function(not_timepoint, sex = NULL, phseq_obj){
    ps_mod <- phseq_obj
    # prune
    if (!is.null(sex)){
        nes_sub_sex <- sample_data(ps_mod)$sex == sex
        ps_mod <- prune_samples(nes_sub_sex, ps_mod)
    }
    nes_sub_time <- !(sample_data(ps_mod)$timepoint == not_timepoint)
    ps_temp <- prune_samples(nes_sub_time, ps_mod)
    # analysis
    timedds <- phyloseq_to_deseq2(ps_temp, ~ timepoint)
    timedds$timepoint
    timedds <- DESeq(timedds, test="Wald", fitType="parametric")
    # formatting
    res <- results(timedds, cooksCutoff = FALSE)
    alpha <- 0.01
    sigtab <- res[which(res$padj < alpha), ]
    sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps3)[rownames(sigtab), ], "matrix"))
    sigtab <- as_tibble(sigtab)
    sigtab$comparison <- stringr::str_remove("T1T2T3", not_timepoint)
    if (!is.null(sex)){
        sigtab$sex <- sex
    }
    sigtab
}

ggplot(sigtab, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
    theme_martin() +
    facet_wrap(~comparison) +
    coord_flip()
ps_mod <- ps3
nes_sub_sex <- sample_data(ps3)$sex == "M"
ps_mod <- prune_samples(nes_sub_sex, ps3)

ps_mod <- subset_samples(ps3, timepoint != "T3")