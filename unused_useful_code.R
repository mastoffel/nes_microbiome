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














# overview over phyloseq options -------------------------------------------------------------------
#ps
#ntaxa(ps)
#nsamples(ps)
#sample_names(ps)[1:5]
#rank_names(ps)
#sample_variables(ps)
# #otu_table(ps)[1:3, 1:3]
#tax_table(ps)[1:5, 1:5]
#taxa_names(ps)[1]


# old microbiome overview plot
# summarize taxa
RAs <- summarize_taxa(ps3, "Phylum") %>% 
    arrange(desc(meanRA)) %>% 
    mutate(Phylum = fct_inorder(Phylum)) %>% 
    filter(count > 1500)

ggplot(RAs, aes(Phylum, meanRA)) +
    geom_point(alpha = 1, size = 2) +
    geom_errorbar(aes(ymin = meanRA-sdRA, ymax = meanRA+sdRA), width = 0.1) +
    theme_martin() +
    theme(axis.text = element_text(angle = 90)) +
    ylab("Relative abundance")


calc_RAs_across_time <- function (timepoint, physeq, Rank) {
    if (timepoint != "all") {
        ps <- prune_samples(sample_data(physeq)$timepoint == timepoint, physeq)
    } else {
        ps <- physeq
    }
    
    # Rank <- enquo(Rank)
    Rank_quo <- enquo(Rank)
    out <- ps %>% 
        summarize_taxa(Rank) %>% 
        as_tibble() %>% 
        arrange(desc(meanRA)) %>%
        filter(count > 200) %>%
        mutate(timepoint = timepoint)
}

Rank <- "Phylum"
all_RAs <- bind_rows(lapply(c("all", "T1", "T2", "T3"), calc_RAs_across_time, ps3, Rank = Rank)) %>% 
    dplyr::mutate(Phylum = fct_inorder(factor(Phylum)))

colpal <- c("grey", wes_palette("Moonrise2", 4, type = "discrete"))

p_comp <- ggplot(all_RAs, aes(Phylum, meanRA, fill = timepoint, col = timepoint, alpha = timepoint)) +
    geom_errorbar(aes(ymin = meanRA-sdRA, ymax = meanRA+sdRA), width = 0.0,
                  position=position_dodge(width=0.5), alpha = 1, lwd = 0.3) +
    geom_point(size = 2.5, shape = 21, colour = "black", lwd = 0.1, position=position_dodge(width=0.5)) +
    theme_martin() +
    scale_alpha_manual(values = c(1, 0.6, 0.6, 0.6), name = "Timepoint") +
    scale_fill_manual(values = colpal, name = "Timepoint") +
    scale_color_manual(values = colpal, name = "Timepoint") +
    scale_x_discrete(limits = rev(levels(all_RAs$Phylum))) +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          axis.title = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20), size = 14),
          axis.text = element_text(size = 11)) +
    coord_flip() +
    ylab("Relative abundance")
p_comp


# plot broad Abundances --------------------------------------------------------
# Taxonomic agglomeration
# How many genera are present after filtering? 135

rank_names(ps3)
# Abundance value plotting
plot_abundance = function(physeq, ylabn = "",
                          Facet = "Genus",
                          Color = "Phylum"){
    mphyseq = psmelt(physeq)
    mphyseq <- subset(mphyseq, Abundance > 0)
    ggplot(data = mphyseq,
           mapping = aes_string(x = "timepoint", y = "Abundance",
                                color = Color, fill = Color)) +
        geom_violin(fill = NA) +
        geom_point(size = 1, alpha = 0.3,
                   position = position_jitter(width = 0.3)) +
        facet_wrap(facets = Facet) + ylab(ylabn) +
        scale_y_log10()
}
plot_abundance(ps3)



# Create df with samples in rows and exact ASV sequence names in columns -------
# full data and each core-microbiota (TODO)
# Extract abundance matrix from the phyloseq object
ASV1 <- as(otu_table(ps3), "matrix")
# transpose if necessary
if(taxa_are_rows(ps3)){ASV1 <- t(ASV1)}
# Coerce to data.frame
OTUdf <- as.data.frame(ASV1)
names(OTUdf)



# split ASV sequences into multiple lines for table plotting (not used ATM) ------------------------
# slice <- function(input, by= 4) { 
#     strlen <- str_length(input) 
#     split_end <- seq(from= 0, to = strlen, by = strlen %/% by)[c(-1)] +1
# 
#     input_split <- unlist(strsplit(input, split=""))
#     split_end[length(split_end)] <- length(input_split)
#     split_start <- c(1, split_end[-length(split_end)] + 1)
#     
#     # split input, add new line and concatenate again
#     input_split2 <- purrr::map2(split_start, split_end, function(x,y) {
#                                             out <- c(input_split[x:y], "\n") # \n
#     })
#     input_split3 <- unlist(input_split2)
#     input_split4 <- stringr::str_c(input_split3[-length(input_split3)], collapse = "")
# }
# tp_asv_table$Sequence <- sapply(tp_asv_table$Sequence, slice, by = 4)
# tp_asv_table <- tp_asv_table %>% 
#     mutate_all(linebreak)

# linesep = c("\\addlinespace") # "\\hline"
# kable(tp_asv_table, format = "latex", booktabs = TRUE, escape = F, align = "l",
#     linesep = linesep) %>% 
#     column_spec(2, width = "30cm") %>% 
#     kable_as_image("core_tp_asv", keep_pdf = TRUE, file_format = "jpeg")




# Ordination individual similarity plot ----------------------------------------

make_subset_plots <- function(iter, p_ord_df, plotcols) {
    to_sample <- sample(all_samples, 10)
    all_samples <<- all_samples[!(all_samples%in%to_sample)]
    p_ord_df_sub <- filter(p_ord_df, individual %in% to_sample)
    
    p_ord_plot_ind <- ggplot(p_ord_df_sub, aes(Axis.1, Axis.2, color = individual)) +
        geom_point(size = 1, alpha = 1) +
        #geom_point(size = 3, alpha = 0.8, aes( color = individual)) +
        # theme_martin() +
        geom_polygon(aes(fill=individual), alpha = 0.4, size = 0.1) +
        theme(panel.grid = element_blank(),
              axis.line = element_line(colour = "black", size = 0.5),
              legend.text=element_text(size=7),
              legend.key.size = unit(0.6,"line"),
              legend.title = element_blank(),
              plot.margin = unit(c(-3,0,0,0), "cm"),
              legend.position = "none",
              axis.title = element_blank()) + #unit(c(-1.1,0.3,-1.1,1), "cm")
        scale_x_continuous(breaks = seq(from = -0.25, to = 0.25, by = 0.25)) + # limits = c(-0.45, 0.45))
        scale_y_continuous(breaks = seq(from = -0.2, to = 0.2, by = 0.2)) + #, limits = c(-0.4, 0.4)
        #scale_shape_manual(values = c(21,2))+
        scale_color_manual(values = plotcols) +
        scale_fill_manual(values = plotcols) +
        coord_fixed(sqrt(evals[2] / evals[1])) +
        theme(legend.position = "bottom",
              legend.direction = "horizontal") +
        # xlab("\n Axis 1 [28,5%]") +
        # ylab("Axis 2 [13,4%] \n ") +
        theme(panel.grid = element_blank(),
              axis.line.x = element_line(colour = "black", size = 0.3, linetype = 1),
              axis.line.y = element_line(colour = "black", size = 0.3, linetype = 1),
              axis.ticks = element_line(colour = "black", size = 0.3)) +
        guides(fill=guide_legend(override.aes=list(shape=21))) +
        # theme(legend.position = "bottom",
        #    legend.direction = "horizontal") +
        xlab("Axis 1 [28,5%]") +
        ylab("Axis 2 [13,4%]") 
    p_ord_plot_ind
}
all_samples <- unique(p_ord_df$individual)
set.seed(24)
plotcols <- sample(c(get_colors("Paired")))

all_plots <- map(1:4, make_subset_plots, p_ord_df, plotcols)

p_all <- all_plots[[1]] + all_plots[[2]] + all_plots[[3]] + all_plots[[4]]
p_all
# ggsave("../figures/Sup2_mds_by_int_2.jpg", p_all, width = 7, height = 5.5)