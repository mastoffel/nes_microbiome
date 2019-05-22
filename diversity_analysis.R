# Full analysis of ASVs for the paper

# load packages ------------------------------------------------------------------------------------
library(pacman)
p_load(phyloseq, tidyverse, msa, inbreedR, rptR, lme4, DESeq2, 
       dada2, phangorn, wesanderson, grid, cowplot, readxl, RColorBrewer, 
       blogdown, patchwork, Demerelate, vegan, ecodist, reshape2, microbiome,
       kableExtra)
source("martin.R")
#library(microbiomeSeq)

# Prepare and load data ----------------------------------------------------------------------------

# folder with ASV table 
input_folder <- "primer_clipped_reads_22_220230_pool"

# load ASV table and ...
load(paste0("output/", input_folder, "/taxa.RData"))
load(paste0("output/", input_folder, "/seqtab.RData"))
load(paste0("output/", input_folder, "/fitGTR.RData"))

# sample names
samples_out <- rownames(seqtab_nochim)
# make data.frame
nes_df <- data.frame("id" = samples_out)

# read and process other northern elephant seal data
nes_sampling <- read_xlsx("../data/processed/sampling_data_processed.xlsx") %>% 
  dplyr::rename(id = ID,
                date = DATE,
                territory = TERRITORY, 
                sex = SEX) %>% 
  mutate(id = str_replace(id, "17BEMA0", "17BEMa")) %>% 
  mutate(id = str_replace(id, "17BEMA", "17BEMa")) %>% 
  mutate(id = str_c(id, timepoint)) %>% 
  mutate(id = str_replace(id, "T1", "")) %>% 
  mutate(sex = ifelse(sex == "M", "F", "M")) %>% 
  mutate(individual = str_replace(id, "T[2,3]", "")) %>% 
  mutate(died = ifelse(individual %in% names(which(table(individual) < 3)), "died", "survived")) %>% 
  dplyr::select(-date) %>% 
  as.data.frame()

# join to nes_df to get the right sample sequence
nes <- left_join(nes_df, nes_sampling, by = "id") %>% 
          arrange(timepoint) %>% 
          mutate(id = fct_inorder(id))

# rownames have to match ASV table
rownames(nes) <- nes$id 

# add birth territory variable
impute_birth_territory <- function(row_num, df){
  sample_row <- df[row_num, ]
  if (str_detect(as.character(sample_row$id), "T")){
    sample_row$birth_territory <- as.factor(df[stringr::str_replace(sample_row$id,"T[2,3]", ""), "territory"])
  } else {
    sample_row$birth_territory <- as.factor(sample_row$territory)
  }
  sample_row
}
nes2 <- bind_rows(lapply(1:nrow(nes), impute_birth_territory, nes))
rownames(nes2) <- nes2$id 

# create phyloseq object
ps0 <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE), 
               sample_data(nes2), 
               tax_table(taxa),
               phy_tree(fitGTR$tree))

# filter out fecal sample
ps0 <- subset_samples(ps0, id != "17BEMa11Fec")


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

# filtering ----------------------------------------------------------------------------------------

# which samples have low abundance?
sort(sample_sums(ps0))

ntaxa(ps0)
# 8 Chloroplast, 3 Mitochondria, NAs on Class level removed too (77)
ps <- ps0 %>%
    subset_taxa(
        ((Family != "Mitochondria") | is.na(Family)) & (Class  != "Chloroplast")
    )
ntaxa(ps)


# Define prevalence of each taxon
# (in how many samples did each taxa appear at least once)
prev0 <- apply(X = otu_table(ps),
                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prev0,
                     TotalAbundance = taxa_sums(ps),
                     tax_table(ps))

# Kingdoms: 2791 Bacteria, 12 Eukaryota, 6 Archaea
sort(table(prevdf$Family, useNA = "always"), decreasing = TRUE)
names(sort(table(prevdf$Family, useNA = "always"), decreasing = TRUE))
# checks
prevdf_check <- as_tibble(prevdf)
hist(prevdf_check$Prevalence, breaks = 100)

# check out how many phyla and whether there are NAs
sort(table(prevdf$Class, useNA = "always"), decreasing = TRUE) # 11 unidentified Phyla --> filter


# prepare prevalence filtering

# Spirochaetae appear and three samples and are pathogenic
sort(table(prevdf$Phylum))
keepPhyla <- table(prevdf$Phylum)[(table(prevdf$Phylum) > 2)] # phylum appears minimum 3 samples
prevdf1 <- subset(prevdf, Phylum %in% names(keepPhyla))

# Keep taxa when appearing in minimum 2% samples (3 samples)
prevalenceThreshold <- 0.02 * nsamples(ps)
prevalenceThreshold

# FILTER 1
# execute prevalence filter minimum 2% samples (3 samples)
ntaxa(ps)
ps1 <- prune_taxa((prev0 > prevalenceThreshold), ps)
ps1
ntaxa(ps1)

# FILTER 2
# ps2 <- subset_taxa(ps1, Phylum %in% names(keepPhyla))
# ps2
# sort(table(tax_table(ps2)[, "Phylum"], exclude = NULL))

# FILTER 3, either of:

# taxa have to have at least 30 reads
plot(sort(taxa_sums(ps), TRUE), type="h", ylim=c(0, 100))

min_reads <- 30
ps3 <- filter_taxa(ps1, function(x) sum(x) > min_reads, TRUE)
ntaxa(ps3)


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
  theme_martin() +
  guides(color=FALSE)
# p_prev

#ggsave("../figures/Sup1_PrevVsAbund.jpg", p_prev, width = 7, height = 6)

# Microbiome composition ---------------------------------------------------------------------------

# plotting abundances -----------------------

# plot relative abundances for all samples
ps_rel <- transform_sample_counts(ps3, function(x) x / sum(x) )
#plot_bar(ps_rel, x = "id", fill="Phylum") +
#    theme(axis.text = element_text(size = 3))

# summarize taxa
source("microbiome_composition_funs.R")
library(forcats)
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

# some other plotting options
ggplot(all_RAs, aes(timepoint, meanRA)) + 
  geom_col(aes(fill = Phylum)) +
  ylab("Relative abundance") +
  xlab("Timepoint") +
  scale_fill_manual(values = colpal, name = "Phylum") +
  coord_flip() 

# 
nes_phylum <- ps3 %>% 
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>% 
  psmelt() %>% 
  filter(Abundance > 0.01) %>% 
  arrange(desc(Abundance))

p_bar <- ggplot(nes_phylum, aes(individual, Abundance, fill = Phylum)) +
  facet_grid(timepoint~.)+
  geom_bar(stat = "identity") +
  theme_martin(base_family = "Helvetica") +
  #scale_fill_manual(values = phylum_colors) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14), 
        axis.text.y = element_text(size=8))+
  scale_fill_manual(values = c(wes_palette("Chevalier1")[1],wes_palette("FantasticFox1")[3:1],  
                               wes_palette("IsleofDogs1")[c(2,1,3)] ), name = "Phylum") +
  #scale_fill_brewer(type = "qual", palette = "Set3") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  theme(legend.text = element_text(size = 11),legend.title = element_text(face="bold")) +
  ylab("Relative abundance \n") +
  xlab("Individuals") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
        #axis.line = element_line(size = 0.5, colour = "black")) +
  theme(panel.spacing.x = unit(10, "lines"))+
  theme(strip.text.y = element_text(size=12)) 
p_bar
ggsave(filename = "../figures/composition_phylum_bar.jpg", p_bar, width = 8, height = 4.5)

#plot_bar(ps3, "Phylum", facet_grid=~timepoint)

#ggsave(filename = "../figures/composition_phylum.jpg", p_comp, width = 5.8, height = 4)

# number of ASVs per sample ----------------------------------------------------
ps_otu <- as.data.frame(otu_table(ps3))
asv_per_sample <- rowSums(ps_otu > 0)
mean(asv_per_sample) # 286
sd(asv_per_sample) # 67
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



# define some extra variables ----------------------------------------------------------------------
# mark low abundance samples and impute birth territory
sample_data(ps)$abundance <- sample_sums(ps) < 10000
sample_data(ps3)$abundance <- sample_sums(ps3) < 10000

# Create data frame with samples in rows and exact ASV sequence names in columns for----------------
# full data and each core-microbiota (TODO)

# Extract abundance matrix from the phyloseq object
ASV1 <- as(otu_table(ps3), "matrix")
# transpose if necessary
if(taxa_are_rows(ps3)){ASV1 <- t(ASV1)}
# Coerce to data.frame
OTUdf <- as.data.frame(ASV1)
names(OTUdf)

# core microbiome across all three timepoints ------------------------------------------------------

# Calculate relative abundances
ps_rel <- microbiome::transform(ps3, "compositional")

# Core Microbiome for each timepoint
core_microbiome <- function(timepoint, ps) {
    # ps_rel_t1 <- subset_samples(ps_rel, timepoint == "T1")
    ps_tp <- prune_samples(sample_data(ps)$timepoint %in% timepoint, ps)
    # get core microbiome (prevalence = 90%)
    core_taxa_tp <- core(ps_tp, detection = 0, prevalence = 90/100)
    # function for assessing means and medians of relative abundances
    x <- otu_table(core_taxa_tp)
    # taxa_are_rows(x)
    taxa_means <- round(colMeans(x), 4)
    # taxa_medians <- apply(x, 2, median)
    # get tax table
    tp_tax_table <- as.data.frame(tax_table(core_taxa_tp)) %>% 
        tibble::rownames_to_column("ASV") %>% 
        mutate(`Mean rel. abundance %` = 100 * taxa_means) %>% 
        arrange(desc(`Mean rel. abundance %`)) 
}

# three timepoints
all_timepoints <- c("T1", "T2", "T3")
core_across_time <- purrr::map(all_timepoints, core_microbiome, ps_rel)
names(core_across_time) <- all_timepoints
# save tax tables with exact sequences
sapply(all_timepoints, function(x) write_excel_csv(core_across_time[[x]], 
                                   path = paste0("../data/processed/core_microbiome_", x, ".txt")))

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



# Core microbiome tables ---------------------------------------------------------------------------
linesep <- c("") 

# core microbiome table
for (i in all_timepoints) {
    core_timepoint <- dplyr::select(core_across_time[[i]], -ASV)
    kable(core_timepoint, format = "latex", booktabs = TRUE, escape = T, linesep = linesep,
        align = "l") %>% 
        kable_styling(latex_options = c( "scale_down")) %>% 
        row_spec(0, bold = T) %>% 
        kable_as_image(paste0("../tables/core_microbiome_", i), keep_pdf = TRUE, file_format = "jpeg")
    
}

library(forcats)
source("martin.R")
# core microbiome plots
# T1
core_T1 <- core_across_time[[1]] %>% 
    dplyr::rename(Abundance = `Mean rel. abundance %`) %>% 
    dplyr::mutate(ASV_id = 1:nrow(.)) %>% 
    mutate(Genus = fct_inorder(Genus))

p1 <- ggplot(core_T1, aes(ASV_id, Abundance)) +
        geom_col(fill = "ghostwhite", col = "#333333") +
        theme_martin() +
        scale_x_reverse(breaks = 1:nrow(core_T1),
                            labels = as.character(core_T1$Genus)) +
        ylab("Mean relative abundance %") +
        xlab("ASV (Genus)") +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        coord_flip() 
ggsave(filename = "../figures/core_t1.jpg", width = 4, height = 4)


core_T2 <- core_across_time[[2]] %>% 
    dplyr::rename(Abundance = `Mean rel. abundance %`) %>% 
    dplyr::mutate(ASV_id = 1:nrow(.)) %>% 
    mutate(Genus = fct_inorder(Genus))

p2 <- ggplot(core_T2, aes(ASV_id, Abundance)) +
    geom_col(fill = "ghostwhite", col = "#333333") +
    theme_martin() +
    scale_x_reverse(breaks = 1:nrow(core_T2),
        labels = as.character(core_T2$Genus)) +
    ylab("Mean relative abundance %") +
    xlab("ASV (Genus)") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    coord_flip() 
ggsave(filename = "../figures/core_t2.jpg", width = 4, height = 4)


core_T3 <- core_across_time[[3]] %>% 
    dplyr::rename(Abundance = `Mean rel. abundance %`) %>% 
    dplyr::mutate(ASV_id = 1:nrow(.)) %>% 
    mutate(Genus = fct_inorder(Genus))

p3 <- ggplot(core_T3, aes(ASV_id, Abundance)) +
    geom_col(fill = "ghostwhite", col = "#333333") +
    theme_martin() +
    scale_x_reverse(breaks = 1:nrow(core_T3),
        labels = as.character(core_T3$Genus)) +
    ylab("Mean relative abundance %") +
    xlab("ASV (Genus)") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    coord_flip() 
ggsave(filename = "../figures/core_t3.png", width = 5, height = 5)

# transformation and ordination --------------------------------------------------------------------

# which transformation is appropriate?
qplot(log(rowSums(otu_table(ps3)))) +
    xlab("Logged counts-per-sample")

# variance stabilizing transformation =========================

# convert to deseq
nes_dds <- phyloseq_to_deseq2(ps3, ~1)
# estimate size factors with not including 0 in geometric mean calc
nes_dds  <- estimateSizeFactors(nes_dds, type = "poscounts") %>% 
                estimateDispersions(fitType = "local")
# create new phyloseq object with variance stabilised ASV table
ps_vst <- ps3
otu_table(ps_vst) <- otu_table(getVarianceStabilizedData(nes_dds), taxa_are_rows = TRUE)

# Ordination sex / time plot ==================================

# for plotting
get_colors <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)

colpal <- c(  "#0000ff", "#ffb14e", "#ea5f94"  )
colpal <- wes_palette("Moonrise2", 3, type = "discrete")   
# calculate ordination
ps_ord <- ordinate(ps_vst, "MDS", "bray")
# calculate axis length relationships according to eigenvalues
evals <- ps_ord$values$Eigenvalues
# get df
#ps_ord <- ordinate(ps_vst, method = "MDS",distance = "wunifrac")
plot_ordination(ps_vst, ps_ord, shape = "sex", color = "timepoint")

p_ord_df <- plot_ordination(ps_vst, ps_ord, shape = "sex", color = "timepoint", justDF = TRUE,
    axes = 1:4)

#which((p_ord_df$Axis.1 > 0) & (p_ord_df$Axis.2 < 0) & (p_ord_df$sex == "F"))
#p_ord_df[32, ]
p_ord_plot <- ggplot(p_ord_df, aes(Axis.1, Axis.2)) +
    geom_point(size = 3.5, alpha = 0.8, aes(shape = sex, fill = timepoint)) +
    #geom_point(size = 3, alpha = 0.8, aes( color = individual)) +
    #theme_martin() +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = colpal, name = "Timepoint") +
    coord_fixed(sqrt(evals[2] / evals[1])) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal") +
    xlab("Axis 1 [28,5%]") +
    ylab("Axis 2 [13,4%]") +
    theme(panel.grid = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3, linetype = 1),
        axis.line.y = element_line(colour = "black", size = 0.3, linetype = 1),
        axis.ticks = element_line(colour = "black", size = 0.3)) +
    guides(fill=guide_legend(override.aes=list(shape=21)))
    
p_ord_plot
ggsave("../figures/Fig1_sex_time_MDS_new.jpg", p_ord_plot, width = 6, height = 4)

# Ordination outlier plot ----------------------------------
p_ord_plot_outlier <- ggplot(p_ord_df, aes(Axis.1, Axis.2)) +
    geom_point(size = 3, alpha = 0.8, aes(color = abundance)) +
    #geom_point(size = 3, alpha = 0.8, aes( color = individual)) +
    theme_martin() +
    theme(panel.grid = element_blank()) +
    scale_x_discrete(labels = c("")) +
    #scale_shape_manual(values = c(21,2))+
    #scale_color_manual(values = colpal) +
    coord_fixed(sqrt(evals[2] / evals[1])) +
    theme(legend.position = "bottom",
        legend.direction = "horizontal") +
    xlab("Axis 1 [28,4%]") +
    ylab("Axis 2 [13,3%]")

p_ord_plot_outlier

# Within individual similarity --------------------------------------------------------------------

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
      xlab("\n Axis 1 [28,5%]") +
      ylab("Axis 2 [13,4%] \n ") +
      theme(panel.grid = element_blank(),
            axis.line.x = element_line(colour = "black", size = 0.3, linetype = 1),
            axis.line.y = element_line(colour = "black", size = 0.3, linetype = 1),
            axis.ticks = element_line(colour = "black", size = 0.3)) +
      guides(fill=guide_legend(override.aes=list(shape=21)))
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
ggsave("../figures/Sup2_mds_by_int_2.jpg", p_all, width = 7, height = 5.5)

# Within individual similarity2 --------------------------------------------------------------------

# p_ord_df %>% 
#   filter(individual %in% sample(unique(p_ord_df$individual), 13)) %>% 
# ggplot(aes(Axis.1, Axis.2, shape = sex, color = individual)) +
#   geom_point() +
#   geom_polygon(alpha = 0.02)
set.seed(25)
plotcols <- sample(c(wes_palette("FantasticFox1")[-4], wes_palette("GrandBudapest2")[-1], wes_palette("Moonrise2")), 6)
# good inds to show
# 17BEMa32, 17BEMa27 or 24, 17BEMa38, 17BEMa3
set.seed(113)
data_sub <- filter(p_ord_df, individual %in% sample(unique(p_ord_df$individual), 8)) %>% 
              filter(!(individual %in% c("17BEMa4", "17BEMa13")))

p_host <- ggplot(p_ord_df, aes(Axis.1, Axis.2, shape = sex)) +
  geom_point(size = 3.5, alpha = 0.6, fill = "grey90", stroke = 0.2) +
  geom_point(data = data_sub, 
             aes(col=individual, fill=individual),
             size = 3.5, alpha = 0.9, stroke = 0.5) +
  geom_polygon(data = data_sub, 
               aes(col=individual, fill=individual), alpha = 0.5, size = 0.5) + # , fill=NA
  scale_shape_manual(values = c(21,24), name = "Sex") +
 scale_x_continuous(breaks = seq(from = -0.25, to = 0.25, by = 0.25)) + # limits = c(-0.45, 0.45))
  scale_y_continuous(breaks = seq(from = -0.2, to = 0.2, by = 0.2)) + #, limits = c(-0.4, 0.4)
  #scale_shape_manual(values = c(21,2))+
  #scale_shape_manual(values = c(21,24), name = "Sex") +
  scale_color_manual(values = c(plotcols), labels = rep("", 6)) +
  scale_fill_manual(values = c(plotcols), labels = rep("", 6)) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  xlab("Axis 1 [28,5%]") +
  ylab("Axis 2 [13,4%]") +
  guides(color= FALSE, shape = FALSE,
         fill = guide_legend("Host", nrow = 1,
                             override.aes=list(fill=plotcols, col = plotcols, shape=c(17,16,16,17,17,16))),
         keywidth = 1) +
  theme(panel.grid = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3, linetype = 1),
        axis.line.y = element_line(colour = "black", size = 0.3, linetype = 1),
        axis.ticks = element_line(colour = "black", size = 0.3),
        legend.position = "bottom") 
 
library(cowplot)
p_full <- plot_grid(p_ord_plot, p_host, labels = c("A", "B"))
p_full 
p_full_vert <-  p_ord_plot / p_host
ggsave("../figures/beta_div.jpg",plot = p_full, width = 10, height = 3.8)
ggsave("../figures/beta_div_vert.jpg",plot = p_full_vert, width = 6, height = 7)



# plotting time trends -----------------------------------------------------------------------------
ps_rel <- transform_sample_counts(ps3, function(x) x/sum(x))
ps_df <- psmelt(ps_rel)
ps_df <- subset(ps_df, Abundance > 0)
#ps_df <- as_tibble(ps_df)

length(unique(ps_df$OTU))


# time trends across Classes ----------------------------------------------------------
plot_df <- ps_df %>% 
    group_by(Sample, Class) %>% 
    dplyr::summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else dplyr::first(.))) %>% 
    dplyr::filter(!is.na(timepoint)) %>% 
    dplyr::filter(Abundance > 0.0001) %>% 
   dplyr::filter(Class %in% c("Actinobacteria", "Bacilli","Bacteroidia", "Betaproteobacteria",
                       "Clostridia", "Coriobacteriia", "Deferribacteres", "Deltaproteobacteria",
                       "Erysipelotrichia", "Flavobacteriia", "Fusobacteriia", "Gammaproteobacteria",
                       "Mollicutes",  "Negativicutes", "Spirochaetes", "Epsilonproteobacteria")) %>%
   #mutate(Abundance = log(Abundance + 0.001))
    mutate(Abundance = log(Abundance)) 

p_class <- ggplot(plot_df , aes(x = timepoint, y = Abundance, by = sex, shape = sex, fill = sex)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    #geom_jitter(size = 2.7, alpha = 0.6,  col = "black", aes(shape = sex, fill = "grey"), width = 0.3, stroke =0.7)
    geom_point(position=position_jitterdodge(jitter.width = 0.1), size=1, alpha=0.6, color = "black") +
   # geom_jitter(alpha=0.3, width = 0.1, aes(by = sex)) +
    facet_wrap(~Class) +
    scale_y_continuous(breaks = log(c(0.001, 0.010, 0.100, 1.000)), labels = c(c("0.001", "0.010", "0.100", "1.000"))) +
    theme_martin() +
    ylab("Relative abundance") +
    xlab("Timepoint") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    guides(fill = guide_legend(),
           shape = guide_legend())

p_class

ggsave("../figures/Fig2_classes_time_trends.jpg", p_class, width = 7.5, height = 6.5)
#ggsave("../figures/Sup3_classes_time_trends_full.jpg", p_class, width = 9.7, height = 7.5)


# time trends across Phyla ------------------------------------
plot_df <- ps_df %>% 
    group_by(Sample, Phylum) %>% 
    summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else dplyr::first(.))) %>% 
    filter(!is.na(timepoint)) %>% 
    dplyr::filter(Abundance > 0.0001) %>% 
    mutate(Abundance = log(Abundance)) 
    # filter(Phylum %in% c("Actinobacteria", "Bacteriodetes", "Deferribacteres", "Firmicutes", "Fusobacteria",
    #                      "Proteobacteria", "Spirochaetae", "Tenericutes")) # %>% 
   # mutate(Abundance = log(Abundance + 0.001)) 

p_phylum <- ggplot(plot_df , aes(x = timepoint, y = Abundance, by = sex, shape = sex, fill = sex)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    #geom_jitter(size = 2.7, alpha = 0.6,  col = "black", aes(shape = sex, fill = "grey"), width = 0.3, stroke =0.7)
    geom_point(position=position_jitterdodge(jitter.width = 0.1), size=1, alpha=0.6, color = "black") +
    # geom_jitter(alpha=0.3, width = 0.1, aes(by = sex)) +
    facet_wrap(~Phylum) +
    scale_y_continuous(breaks = log(c(0.001, 0.010, 0.100, 1.000)), labels = c(c("0.001", "0.010", "0.100", "1.000"))) +
    theme_martin() +
    ylab("Relative abundance") +
    xlab("Timepoint") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    guides(fill = guide_legend(),
        shape = guide_legend())

p_phylum

ggsave("../figures/Sup4_phylum_time_trends.jpg", p_phylum, width = 7.5, height = 6.5)

# time trends across Order--------------------------------------------------------------------------
plot_df <- ps_df %>% 
    group_by(Sample, Order) %>% 
    summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else dplyr::first(.))) %>% 
    filter(!is.na(timepoint)) %>% 
    dplyr::filter(Abundance > 0.0001) %>% 
    filter(Order %in% c("Actinomycetales", "Aeromonodales", "Anaeroplasmatales",
        "Bacillales", "Bacteroidales", "Burkholderiales", "Campylobacterales",
        "Clostridiales", "Coriobacteriales", "Corynebacteriales", "Desulfovibrionales",
        "Enterobacteriales", "Erysipelotrichales", "Flavobacteriales", "Fusobacteriales",
        "Lactobacillales", "Mycoplasmatales", "Neisseriales", "Pasteurellales",
        "Pseudomonadales", "Rhodospirillaleles", "Selenomonadales", "Spirochaetales")) %>%
    mutate(Abundance = log(Abundance)) 

p_order <- ggplot(plot_df , aes(x = timepoint, y = Abundance, by = sex, shape = sex, fill = sex)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    #geom_jitter(size = 2.7, alpha = 0.6,  col = "black", aes(shape = sex, fill = "grey"), width = 0.3, stroke =0.7)
    geom_point(position=position_jitterdodge(jitter.width = 0.1), size=1, alpha=0.6, color = "black") +
    # geom_jitter(alpha=0.3, width = 0.1, aes(by = sex)) +
    facet_wrap(~Order) +
    scale_y_continuous(breaks = log(c(0.001, 0.010, 0.100, 1.000)), labels = c(c("0.001", "0.010", "0.100", "1.000"))) +
    theme_martin() +
    ylab("Relative abundance") +
    xlab("Timepoint") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    guides(fill = guide_legend(),
        shape = guide_legend()) +
    theme(strip.text = element_text(size=9),
          axis.text = element_text(size = 9))

p_order

ggsave("../figures/Sup5_order_time_trends.jpg",p_order, width = 8.5, height = 7.5)


# time trends across Genera ------------------------------
# get top 10 core microbiota across all three timepoints
top10 <- as.character(unlist(map(core_across_time, function(x) x$Genus[1:10])))
top10 <- unique(top10)

plot_df <- ps_df %>% 
    group_by(Sample, Genus) %>% 
    dplyr::summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else dplyr::first(.))) %>% 
    dplyr::filter(!is.na(timepoint)) %>% 
    dplyr::filter(Abundance > 0.0001) %>% 
    dplyr::filter(Genus %in% top10) %>%
    mutate(Abundance = log(Abundance + 0.00001)) %>% 
    dplyr::filter(!is.na(Genus))
    #mutate(Abundance = sqrt(Abundance)) 

p_genus <- ggplot(plot_df , aes(x = timepoint, y = Abundance, by = sex, shape = sex, fill = sex)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    #geom_jitter(size = 2.7, alpha = 0.6,  col = "black", aes(shape = sex, fill = "grey"), width = 0.3, stroke =0.7)
    geom_point(position=position_jitterdodge(jitter.width = 0.1), size=1, alpha=0.6, color = "black") +
    # geom_jitter(alpha=0.3, width = 0.1, aes(by = sex)) +
    facet_wrap(~Genus) +
    scale_y_continuous(breaks = log(c(0.001, 0.010, 0.100, 1.000)), labels = c(c("0.001", "0.010", "0.100", "1.000"))) +
    theme_martin() +
    ylab("Relative abundance") +
    xlab("Timepoint") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    guides(fill = guide_legend(),
        shape = guide_legend())

p_genus

ggsave("../figures/Fig2_core_microbiota_genus_time_trends.jpg", p_genus, width = 7.5, height = 5.5)


# alpha diversity ----------------------------------------------------------------------------------
diversity_df <- estimate_richness(ps, measures = c("Shannon", "Simpson", "InvSimpson", "Observed", "Fisher")) %>% 
                    tibble::rownames_to_column("id") %>% 
                    mutate(id = str_replace(id, "X", "")) %>% 
                    left_join(as_tibble(sample_data(ps)), by = "id")

colpal_cavalanti <- wes_palette("Cavalcanti1", 2, type = "discrete")
as.character(wes_palette("Darjeeling2"))
colpal_moonrise <- c("#899DA4",  "#79402E")

#library(tidyr)
#div_lf <- gather(diversity_df, div_meas, div, c(Observed, Shannon, Simpson, InvSimpson))
set.seed(12)
p_div <- ggplot(diversity_df, aes(timepoint, Shannon, by = sex)) + #colour = sex
    geom_boxplot(alpha = 0.4, outlier.shape = NA, aes(fill = sex)) + #, aes(color = sex)
    geom_point(position=position_jitterdodge(jitter.width = 0.3), size = 2.3, alpha = 0.8,  
        col = "black", aes(shape = sex, fill = sex), stroke =0.7) +
   # geom_jitter(size = 2.3, alpha = 0.8,  col = "black", aes(shape = sex, fill = sex), width = 0.3, stroke =0.7) + #shape = 21,
    #scale_shape_manual(values = c(21,24))+
    #facet_wrap(~timepoint) +
    theme_martin(base_family = "Helvetica", highlight_family = "Helvetica") +
   # ggtitle("Timepoint") +
    #scale_fill_manual(values =  "grey") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    #scale_x_discrete(labels = c(stri_unescape_unicode("\\u2640"), stri_unescape_unicode("\\u2642"))) +
    #scale_fill_manual(values = colpal_moonrise, name = "Sex") +
    #scale_color_manual(values = colpal_moonrise, name = "Sex") +
    ylab("Shannon diversity\n")+
    xlab("\nTimepoint")+
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line = element_line(colour = "black", size = 0.3, linetype = 1),
          axis.ticks = element_line(colour = "black", size = 0.3, linetype = 1))
    #scale_x_discrete(labels = stri_unescape_unicode("a\\u0105!\\u0032\\n")) 
p_div
ggsave("../figures/Fig4_diversity.jpg", p_div, width = 4.5, height = 2.9)

# modeling 
library(lme4)
library(partR2)
div_mod <- lmer(Shannon ~ sex + timepoint + (1|individual), data = diversity_df)
summary(div_mod)
tidy(div_mod)
VarCorr(div_mod)
set.seed(17) 
boot_div_mod <- confint(div_mod, method = c("boot"), nsim = 1000)
boot_div_mod
R2_div <- partGaussian(div_mod, partvars = c("sex", "timepoint"), nboot = 1000) #
rpt_div <- rptGaussian(Shannon ~ (1|individual),grname = "individual", data = diversity_df)#"timepoint"


# permanova analysis (check pseudoreplication problem) --------------------------------------------
library(vegan)
# overall effects 
metadata <- as(sample_data(ps_vst), "data.frame")
mod_full <- adonis(phyloseq::distance(ps_vst, method="bray") ~ timepoint + sex + individual,
    data = metadata, by = "terms")
mod_full
# t1t2
ps_vst_t1t2 <- subset_samples(ps_vst, timepoint %in% c("T1", "T2"))
metadata_t1t2 <- as(sample_data(ps_vst_t1t2), "data.frame")
mod_t1t2 <- adonis(phyloseq::distance(ps_vst_t1t2, method="bray") ~ timepoint + sex + individual,
    data = metadata_t1t2, by = "terms")
# t2t3
ps_vst_t2t3 <- subset_samples(ps_vst, timepoint %in% c("T2", "T3"))
metadata_t2t3 <- as(sample_data(ps_vst_t2t3), "data.frame")
mod_t2t3 <- adonis(phyloseq::distance(ps_vst_t2t3, method="bray") ~ timepoint + sex + individual,
    data = metadata_t2t3, by = "terms")

# individuals, repeatability?
metadata <- as(sample_data(ps_vst), "data.frame")
mod_rpt <- adonis(phyloseq::distance(ps_vst, method="bray") ~  sex + individual, 
    data = metadata)

# group dispersion assumption testing
mod <- betadisper(d = phyloseq::distance(ps_vst, method="bray"), group = metadata$sex)
anova(mod) # parametric anova
TukeyHSD(mod) # Tukey Honest Significant Difference Method
permutest(mod) # Permutation test
mod <- betadisper(d = phyloseq::distance(ps_vst, method="bray"), group = metadata$timepoint)
anova(mod)
TukeyHSD(mod)
permutest(mod)

# 
sample_data(ps_vst)$individual <- as.factor(sample_data(ps_vst)$individual)
p_ord_rpt <- plot_ordination(ps_vst, ps_ord, color = "individual", type = "samples")
p_ord_rpt  + geom_polygon(aes(fill=individual), alpha = 0.05) + geom_point(size=1) + ggtitle("samples") +
    theme_martin()


# experimental repeatability for distance matrices -------------------------------------------------
library(reshape2)
library(vegan)
metadata <- as(sample_data(ps_vst), "data.frame")
mod_ind <- adonis(phyloseq::distance(ps_vst, method="bray") ~    individual, #timepoint + sex +
    data = metadata) # ,  strata = metadata$sex
mod_df <- as_tibble(mod_ind$aov.tab, rownames="varcomp")

# aov method from the repeatability paper
groups <- factor(metadata$individual)
k  <- length(levels(groups))  # number of groups
N  <- ncol(otu_table(ps_vst)) # number of samples
# Anova repeatability
n0 <- 1/(k-1)*(N - sum(table(metadata$individual)^2)/N)
# mean between ind squares
MSa <- as.numeric(mod_df[mod_df$varcomp == "individual", "MeanSqs"])
MSw <- as.numeric(mod_df[mod_df$varcomp == "Residuals", "MeanSqs"])
# point estimate 0.3516666
rpt_res <- (MSa - MSw)/ (MSa + (n0 - 1) * MSw)
# CI      0.1365467 0.5667866
se  <- sqrt((2*(N-1)*(1-rpt_res)^2 * (1 + (n0-1)*rpt_res)^2) / (n0^2 * (N-k)*(k-1)))
CI <- 0.95
rpt_CI <- rpt_res + c(1,-1) * qt((1-CI)/2,k-1)*se 



# deseq2 modeling and plotting (1) -----------------------------------------------------------------

run_deseq_per_sex <- function(sex){
    # run two models for females and males seperately to account for repeated measures
    ps_m <- prune_samples(sample_data(ps3)$sex %in% sex, ps3)
    # create deseq object
    nes_mod_dds_m <- phyloseq_to_deseq2(ps_m, ~ individual + timepoint)
    # estimate size factos
    vst_dds_m <- estimateSizeFactors(nes_mod_dds_m , type = "poscounts")
    # colData(nes_mod_dds)
    # see again if groups differ =========
    #vst_dds_temp <- varianceStabilizingTransformation(vst_dds_m, fitType = "local")
    #plotPCA(vst_dds_temp, intgroup="group")
    # deseq analysis    ===========
    vst_dds <- DESeq(vst_dds_m, test="Wald", fitType="local")
    #colData(vst_dds)
    #resultsNames(vst_dds)
    # Investigate results table
    get_dds_results_tables <- function(contrast, deseq_analysis_object, alpha) {
        res <- results(deseq_analysis_object, cooksCutoff = FALSE, contrast)
        res <- res[order(res$padj, na.last=NA), ]
        sigtab <- res[which(res$padj < alpha), ]
        sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps3)[rownames(sigtab), ], "matrix"))
        sigtab <- as_tibble(sigtab)
        sigtab$comparison <- paste0(contrast[2], "vs" , contrast[3])
        sigtab
    }
    all_contrasts <- list(c("timepoint", "T2", "T1"), c("timepoint", "T3", "T2")) 
    all_sigtabs <- map_df(all_contrasts,get_dds_results_tables, vst_dds, 0.01)
    all_sigtabs$sex <- sex
    all_sigtabs
}

# run models with data subsetted per sex to be able to control for individuals
sigtabs_per_sex <- map_df(c("F", "M"), run_deseq_per_sex)
sigtabs_per_sex$combination <- paste0(sigtabs_per_sex$sex, sigtabs_per_sex$comparison)


# releveling and changing NA to Not assigned for plotting
all_sigtabs_plot <- sigtabs_per_sex %>% 
    mutate(combination = combination %>% as_factor() %>% 
            fct_relevel("FT2vsT1", "FT3vsT2", "MT2vsT1", "MT3vsT2")) %>% 
    tidyr::replace_na(replace = list(Class = "Not assigned", Family = "Not assigned", 
        Genus = "Not assigned", Order = "Not assigned")) %>% 
    mutate_at(c("Family", "Genus", "Order"), function(x) {
        as_factor(x) %>% fct_relevel("Not assigned", after = Inf) %>% fct_rev()
    }) %>% 
    mutate_at(c("Class"), function(x) {
        as_factor(x) %>% fct_relevel("Not assigned", after = 0) %>% fct_rev()
    })

# save Class levels for later plotting
class_levels <- levels(all_sigtabs_plot$Class)
library(RColorBrewer)
get_colors <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)

set.seed(2049) # 5 / 559 / 1053
plotcols <- sample(c(get_colors("Paired"), "grey", "black", "white"))
#plotcols <- c(sample(c(get_colors("Paired"), "grey", "black")),  "white")

# sort colors
# "#A6CEE3" #lightblue
    # "#FB9A99" #lightred
    # "#B2DF8A" #lightgreen
    # "#B15928" #brown
    # "#6A3D9A" #purple
    # "#33A02C" #green
    # "#FF7F00" #orange
    # "#CAB2D6" #lightpurple
# "#E31A1C" #red
    # "#FFFF99" #yellow
    # "#FDBF6F" #lightorange
    # "#1F78B4" #blue
    # "black"
# "grey"
    # "white"

# make a color table

taxa_classes <- as.character(rev(levels(all_sigtabs_plot$Class)))
taxa_colors <- c("#1F78B4", "#FF7F00", "lightgrey", "#B15928", "#CAB2D6", "#6A3D9A",
                "black", "#33A02C", "white", "#FB9A99", "#FDBF6F", "#B2DF8A", "#FFFF99", 
                "#E31A1C", "#A6CEE3")
names(taxa_colors) <- taxa_classes 


dat_text <- data.frame(
    label = c("","\u2640","","\u2642"),
    comparison  = factor(c("FT2vsT1", "FT3vsT2", "MT2vsT1", "MT3vsT2")),
    y = c(Inf, Inf, Inf, Inf),
    x = c(15,15,15,15)
)

# create labeller function
wrap_labels <- c(FT2vsT1 = "T1 => T2", FT3vsT2 = "T2 => T3", MT2vsT1 = "T1 => T2", MT3vsT2 = "T2 => T3")

p_diff_time <- ggplot(all_sigtabs_plot, aes(x=Family, y=log2FoldChange)) + 
    #geom_point(size = 3) + 
    geom_point(size=3, alpha = 0.7, shape = 21, colour = "black", lwd = 0.1, aes(fill=Class)) + 
    theme_martin() +
    facet_wrap(~combination, labeller = labeller(combination = wrap_labels)) +
    scale_fill_manual(values = taxa_colors) +
    # limits = unique(rev(all_sigtabs$Family))
    scale_y_continuous(breaks = seq(from = -30, to = 30, by = 10), limits = c(-35, 35)) +
    geom_hline(yintercept = 0.0, alpha = 0.5, linetype = 2) +
    #annotate("text", x = 10, y = 40, label = "\u2642", size = 20, lwd = 3, color = "darkgrey") +
    coord_flip() +
    ylab(expression(log[2]~fold~change)) +
    theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 15),
       # legend.position = c(0.37, -0.18),
        plot.margin = margin(b = 5, r = 10, l = 5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        strip.text = element_text(size=15),
        panel.spacing = unit(1, "lines")
    ) +
    guides(fill=guide_legend(
        reverse = TRUE,
        #title.position = "left",
        keywidth=0.1,
        keyheight=0.1,
       # nrow = 5,
        default.unit="inch")) 

    #annotate(geom = "text", x = 15, y = Inf, label = "Some text")
    # geom_text(
    #     data    = dat_text,
    #     mapping = aes(x = x, y = y, label = label),
    #     size = 20,
    #     color = "darkgrey"
    #     #hjust   = -0.1,
    #     #vjust   = -1
    # )

p_diff_time_final <- ggdraw(p_diff_time) + 
    draw_label("\u2642",size = 35, colour = "black", x = 0.476, y = 0.3,lineheight = 5) +
    draw_label("\u2640",size = 35, colour = "black", x = 0.476, y = 0.75,lineheight = 5) 
p_diff_time_final 
ggsave("../figures/timepoint_diff_abund_family.jpg", p_diff_time_final, width = 11, height = 8)


# Summary statistics of differential abundance for paper -------------------------------------------
# Overall number of differentially abundant taxa
sigtabs_per_sex %>% 
    group_by(sex, combination) %>% 
    tally()

sigtabs_per_sex %>% 
    group_by(sex, combination) %>% 
    filter(log2FoldChange < 0) %>% 
    tally()

sigtabs_per_sex %>% 
    group_by(sex, combination, Class) %>% 
    tally() %>% 
    arrange(sex, combination, desc(n)) %>% 
    mutate(n_prop = n / sum(n)) %>% 
    print(n = Inf)

sigtabs_per_sex %>% 
    group_by(sex, combination, Family) %>% 
    tally() %>% 
    arrange(sex, combination, desc(n)) %>% 
    mutate(n_prop = n / sum(n)) %>% 
    print(n = Inf)
#


    #summarise(sum_fold_change = sum(log2FoldChange))

# deseq2 modeling and plotting (2) -------------------------------------------------------------------------------

# create factor representing the combination of time and sex
sample_data(ps3)$group <- factor(paste0(sample_data(ps3)$timepoint, sample_data(ps3)$sex))

nes_mod_dds_sex <- phyloseq_to_deseq2(ps3, ~ group)
# estimate size factos
vst_dds_sex <- estimateSizeFactors(nes_mod_dds_sex, type = "poscounts")
colData(nes_mod_dds_sex)

# see again if groups differ 
#vst_dds_temp <- varianceStabilizingTransformation(vst_dds_sex , fitType = "local")
#plotPCA(vst_dds_temp, intgroup="group")

# deseq analysis  
vst_dds_sex_mod <- DESeq(vst_dds_sex, test="Wald", fitType="local")
colData(vst_dds_sex_mod)
resultsNames(vst_dds_sex_mod)
# Investigate results table
get_dds_results_tables <- function(contrast, deseq_analysis_object, alpha) {
    res <- results(deseq_analysis_object, cooksCutoff = FALSE, contrast)
    res <- res[order(res$padj, na.last=NA), ]
    sigtab <- res[which(res$padj < alpha), ]
    sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps3)[rownames(sigtab), ], "matrix"))
    sigtab <- as_tibble(sigtab)
    sigtab$comparison <- paste0(contrast[2], "vs" , contrast[3])
    sigtab
}

all_contrasts <- list(c("group", "T1F", "T1M"),c("group", "T2F", "T2M"), c("group", "T3F", "T3M")) 
all_sigtabs <- map_df(all_contrasts,get_dds_results_tables, vst_dds_sex_mod, 0.01)

# releveling and changing NA to Not assigned for plotting
all_sigtabs_plot <- all_sigtabs %>% 
    mutate(combination = comparison %>% as_factor()) %>% 
            #fct_relevel("FT2vsT1", "FT3vsT2", "MT2vsT1", "MT3vsT2")) %>% 
    tidyr::replace_na(replace = list(Class = "Not assigned", Family = "Not assigned", 
        Genus = "Not assigned", Order = "Not assigned")) %>% 
    mutate_at(c("Family", "Genus", "Order"), function(x) {
        as_factor(x) %>% fct_relevel("Not assigned", after = Inf) %>% fct_rev()
    }) %>% 
    mutate_at(c("Class"), function(x) {
        as_factor(x) %>% fct_relevel("Not assigned", after = 0) %>% fct_rev()
    })
get_colors <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)

#set.seed(2049) # 5 / 559 / 1053
#plotcols <- sample(c(get_colors("Paired"), "grey", "black", "white"))
#plotcols <- c(sample(c(get_colors("Paired"), "grey", "black")),  "white")

# dat_text <- data.frame(
#     label = c("","\u2640","","\u2642"),
#     comparison  = factor(c("FT2vsT1", "FT3vsT2", "MT2vsT1", "MT3vsT2")),
#     y = c(Inf, Inf, Inf, Inf),
#     x = c(15,15,15,15)
# )

# create labeller function
wrap_labels <- c(T1FvsT1M = "T1: F vs. M", T2FvsT2M = "T2: F vs. M", T3FvsT3M = "T3: F vs. M")
#levels(all_sigtabs_plot$Class)

p_diff_sex <- ggplot(all_sigtabs_plot, aes(x=Family, y=log2FoldChange)) + 
    #geom_point(size = 3) + 
    geom_point(size=3, alpha = 0.7, shape = 21, colour = "black", lwd = 0.1, aes(fill=Class)) + 
    theme_martin() +
    facet_wrap(~combination, labeller = labeller(combination = wrap_labels)) + #, 
    scale_fill_manual(values = taxa_colors) +
    # limits = unique(rev(all_sigtabs$Family))
    scale_y_continuous(breaks = seq(from = -30, to = 30, by = 10), limits = c(-35, 35)) +
    geom_hline(yintercept = 0.0, alpha = 0.5, linetype = 2) +
    #annotate("text", x = 10, y = 40, label = "\u2642", size = 20, lwd = 3, color = "darkgrey") +
    coord_flip() +
    ylab(expression(log[2]~fold~change)) +
    theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 15),
        # legend.position = c(0.37, -0.18),
        plot.margin = margin(b = 5, r = 10, l = 5),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15),
        strip.text = element_text(size=15),
        panel.spacing = unit(1, "lines")
    ) +
    guides(fill=guide_legend(
        #title.position = "left",
        reverse = TRUE,
        keywidth=0.1,
        keyheight=0.1,
        # nrow = 5,
        default.unit="inch")) 
p_diff_sex
#annotate(geom = "text", x = 15, y = Inf, label = "Some text")
# geom_text(
#     data    = dat_text,
#     mapping = aes(x = x, y = y, label = label),
#     size = 20,
#     color = "darkgrey"
#     #hjust   = -0.1,
#     #vjust   = -1
# )
ggsave("../figures/timepoint_diff_abund_sex.jpg", p_diff_sex, width = 11, height = 5)



# Summary statistics of differential abundance for paper -------------------------------------------
# Overall number of differentially abundant taxa
all_sigtabs %>% 
    group_by(comparison) %>% 
    tally()

all_sigtabs  %>% 
    group_by(comparison) %>% 
    filter(log2FoldChange > 0) %>% 
    tally()

all_sigtabs %>% 
    group_by(comparison, Class) %>% 
    tally() %>% 
    arrange(comparison, desc(n)) %>% 
    mutate(n_prop = n / sum(n)) %>% 
    print(n = Inf)







# gut microbiota and relatedness -------------------------------------------------------------------

# microsatellites preparation
nes_msats <- readxl::read_xls("../data/raw/nes_msats_cleaned.xls") %>% 
    add_column(factor(rep("SB", nrow(.))), .after = 1)
names(nes_msats)[1:2] <- c("Sample-ID", "Population")
nes_msats <- data.frame(nes_msats)
str(nes_msats)

# test whether number of loci is enough -----------------------------------
rel_test <- Loci.test(nes_msats, bt=1000, object=TRUE, value= "loiselle", file.output=TRUE)
rel_means <- unlist(lapply(rel_test, mean))
yplus <- rel_means + unlist(lapply(rel_test, sd))
yminus <- rel_means - unlist(lapply(rel_test, sd))
rel_test_df <- data.frame(rel_mean = rel_means, lower_sd = yminus, upper_sd = yplus, loci = c(1:21))

p_rel_test <- ggplot(rel_test_df, aes(loci, rel_mean)) +
    geom_pointrange(aes(ymin = yminus, ymax = yplus), fatten = 8, alpha = 1, fill = "lightgrey", shape=21) +
    theme_martin() +
    ylab("Mean difference in relatedness") +
    xlab("Number of loci") +
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20))
p_rel_test
#ggsave("../figures/Sup6_rel_test.jpg",p_rel_test, width = 5, height = 3.5)

# (1) relatedness estimation ----------------------------------------------
nes_rel <- Demerelate(as.data.frame(nes_msats), value= "loiselle", file.output=FALSE,  #loiselle
                      object = TRUE, pairs = 100)

hist(unlist(nes_rel$Empirical_Relatedness))

# Format 
nes_rel_pair_names <- names(unlist(nes_rel[[2]])) %>% 
                        str_replace("SB.", "") %>% 
                        str_replace_all("E", "17BEMa") %>% 
                        str_split("_") %>% 
                        map(function(x) as_tibble(t(x))) %>% 
                        bind_rows()
nes_rel_df <- data.frame("ind1" = nes_rel_pair_names[[1]], "ind2" = nes_rel_pair_names[[2]], 
                        "rel" = as.numeric(unlist(nes_rel$Empirical_Relatedness)))

# WriteXLS(nes_rel_df, "../../../../../Desktop/nes_relatedness.xls")
# microbial distances estimation

# merge samples coming from one individual
#ps_merged <- subset_samples(ps3, sex == "M")
ps_merged <- merge_samples(ps3, "individual", fun = mean)

# ps_merged <- subset_taxa(ps_merged, Order %in% c("Bacteroidales", "Clostridiales"))
# convert to deseq
ps_merged_dds <- phyloseq_to_deseq2(ps_merged, ~ 1) # here ps_merged
# estimate size factors with not including 0 in geometric mean calc
ps_merged_dds <- estimateSizeFactors(ps_merged_dds, type = "poscounts") %>% 
    estimateDispersions(fitType = "local")
# create new phyloseq object with variance stabilised ASV table
ps_merged_vst <- ps3
otu_table(ps_merged_vst) <- otu_table(getVarianceStabilizedData(ps_merged_dds), taxa_are_rows = TRUE)

# create dataframe for each sex and combine for plotting
get_dist_df <- function(sub_sex, ps){
    # subset Males
  #  ps_merged_vst_sub  <- subset_samples(ps, (sex == sub_sex))
    ps_merged_vst_sub <- prune_samples(sample_data(ps)$sex == sub_sex, ps)
    # calculate bray curtis dissimilarity
    nes_mic_dist <- phyloseq::distance(ps_merged_vst_sub, method = "bray")
    labels(nes_mic_dist)
    # create df
    nes_mic_dist_df <- melt(as.matrix(nes_mic_dist), varnames = c("ind1", "ind2")) %>% 
        mutate(ind1 = str_replace(ind1, "T3", "")) %>% 
        mutate(ind2 = str_replace(ind2, "T3", "")) 
    # create joint df with microbial and genetic distance
    both_dist_df <- inner_join(nes_mic_dist_df, nes_rel_df, by = c("ind1" = "ind1", "ind2" = "ind2")) %>%
        mutate(sex = sub_sex)
    # mutate(rel_classes = cut(rel, breaks = quantile(rel, probs = seq(0, 1, 0.15))))
}

sex_distances <- map_df(c("F", "M"), get_dist_df, ps_merged_vst)

p_rel <- ggplot(sex_distances , aes(rel, 1-value)) +
           geom_point(size = 2, alpha = 0.5, aes(shape = sex, fill = sex)) +
           geom_smooth(method = "lm", se = FALSE, aes(color = sex)) +
           facet_wrap(~sex) +
           theme_martin() +
           scale_shape_manual(values = c(21,24), name = "Sex") +
           scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
           scale_color_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
           xlab("\nGenetic relatedness") +
           ylab("Microbiome similarity\n") +
           theme( strip.background = element_blank(),
                panel.grid = element_blank(),
               strip.text.x = element_blank(), 
               axis.text = element_text(color = "black"),
               axis.line = element_line(colour = "black", size = 0.3),
               axis.ticks = element_line(colour = "black", size = 0.3),
               panel.spacing = unit(2, "lines")) 
           
p_rel
ggsave("../figures/Fig5_relatedness.jpg", p_rel, width = 7, height = 3)

# nes_males <- unique(unlist(both_dist_df[c(1,2)]))
# dist_mat <- matrix(nrow = length(nes_males), ncol = length(nes_males), dimnames = list(nes_males, nes_males))

create_distmat <- function(dist_df, sample_var1, sample_var2, value_var) {
    unique_ids <- unique(unlist(dist_df[c(sample_var1, sample_var2)]))
    dist_mat <- matrix(nrow = length(unique_ids), ncol = length(unique_ids), 
        dimnames = list(unique_ids, unique_ids))
    
    for (i in 1:nrow(dist_df)) {
        x <- dist_df[i, ]
        dist_mat[x[[sample_var1]], x[[sample_var2]]]  <- x[[value_var]]
    }
    
    t(dist_mat)
}

# calculate mantel test seperately for the two sexes
both_dist_df <- filter(sex_distances, sex == "M")

distmat_microbes <- as.dist(1 - create_distmat(both_dist_df, "ind1", "ind2", "value"))
distmat_msats <- as.dist(create_distmat(both_dist_df, "ind1", "ind2", "rel"))

mantel_test <- ecodist::mantel(formula =  distmat_microbes ~ distmat_msats, 
    mrank = FALSE, pboot = 0.9, nperm = 10000, nboot = 10000)

# M: mantelr = 0.2623, CI 0.05882922 0.41531428, p two tailed = 0.00090000
# F: mantelr = 0.0588, CI -0.0676119  0.2006059 , p two tailed = 0.4209000
mantel_test 

# check whether slopes are different to avoid some sort of fallacy that Holger mentioned
# idea: permutation test. Permute pairwise relatedness and recalculate interaction slope
library(simpleboot)
mod_rel <- lm(rel ~ value*sex, data = sex_distances)
boot_rel <- lm.boot(mod_rel, 1000, rows = FALSE)
perc(boot_rel , p = c(0.025, 0.975)) # interaction slope -0.11 CI: -[0.23, -0.006]
org_slope <- summary(lm(rel ~ value*sex, data = sex_distances))$coefficients[4, 1]

permed_slope <- function(iter, sex_distances) {
  data_perm <- sex_distances
  data_perm[, "rel"] <- sex_distances$rel[sample(1:nrow(sex_distances))]
  out <- summary(lm(value ~ rel*sex, data = data_perm))$coefficients[4, 1]
  out
}

perm_slopes <- sapply(1:1000, permed_slope, sex_distances)
hist(perm_slopes)
p_val <- 1 - sum(perm_slopes > org_slope) / length(perm_slopes)
org_slope

# how many substances explain relatedness ?

# for later use also
create_distmat <- function(dist_df, sample_var1, sample_var2, value_var) {
    unique_ids <- unique(unlist(dist_df[c(sample_var1, sample_var2)]))
    dist_mat <- matrix(nrow = length(unique_ids), ncol = length(unique_ids), 
        dimnames = list(unique_ids, unique_ids))
    
    for (i in 1:nrow(dist_df)) {
        x <- dist_df[i, ]
        dist_mat[x[[sample_var1]], x[[sample_var2]]]  <- x[[value_var]]
    }
    
    t(dist_mat)
}

calc_cor_rel_mic <- function(topx, ps3, sex) {
    
    # find most abundant taxa
    # ps_rel <- transform_sample_counts(ps3, function(x) log(x))
    # topXX <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:topx]
    # ps_topXXX <- prune_taxa(topXX, ps3)
    
    # merge samples coming from one individual
    #ps_merged <- subset_samples(ps3, sex == "M")
    ps_merged <- merge_samples(ps3, "individual", fun = mean)
    # ps_merged <- subset_taxa(ps_merged, Order %in% c("Bacteroidales", "Clostridiales"))
    # convert to deseq
    ps_merged_dds <- phyloseq_to_deseq2(ps_merged , ~1) # here ps_merged
    # estimate size factors with not including 0 in geometric mean calc
    ps_merged_dds <- estimateSizeFactors(ps_merged_dds, type = "poscounts") %>% 
        estimateDispersions(fitType = "local")
    # create new phyloseq object with variance stabilised ASV table
    ps_merged_vst <- ps3
    otu_table(ps_merged_vst) <- otu_table(getVarianceStabilizedData(ps_merged_dds), taxa_are_rows = TRUE)
    
    # subset Males
    logi_sex <- sample_data(ps_merged_vst)$sex == sex
    ps_merged_vst_sub <- prune_samples(logi_sex, ps_merged_vst)
    #ps_merged_vst_sub2  <- subset_samples(ps_merged_vst, (sex == "F" ))
    
    # subset most abundant
    topXX <- names(sort(taxa_sums(ps_merged_vst_sub), decreasing = TRUE))[1:topx]
    ps_merged_vst_sub <- prune_taxa(topXX, ps_merged_vst_sub)
    
    # calculate bray curtis dissimilarity
    nes_mic_dist <- phyloseq::distance(ps_merged_vst_sub, method = "bray")
    labels(nes_mic_dist)
    
    
    nes_mic_dist_df <- melt(as.matrix(nes_mic_dist), varnames = c("ind1", "ind2")) %>% 
                            mutate(ind1 = str_replace(ind1, "T3", "")) %>% 
                            mutate(ind2 = str_replace(ind2, "T3", "")) 
    
    both_dist_df <- inner_join(nes_mic_dist_df, nes_rel_df, by = c("ind1" = "ind1", "ind2" = "ind2")) %>% 
                        mutate(rel_classes = cut(rel, breaks = quantile(rel, probs = seq(0, 1, 0.15))))
    
    
    distmat_microbes <- as.dist(1 - create_distmat(both_dist_df, "ind1", "ind2", "value"))
    distmat_msats <- as.dist(create_distmat(both_dist_df, "ind1", "ind2", "rel"))
    
    mantel_test <- ecodist::mantel(formula =  distmat_microbes ~ distmat_msats, 
                                   mrank = FALSE, nperm = 1000, nboot = 1000)
    out <- data.frame(t(mantel_test), "topX" = topx)
}

## takes some time
#out_f <- map_df(seq(from = 4, to = 1063, by = 2), calc_cor_rel_mic, ps3, "F")
#out_m <- map_df(seq(from = 4, to = 1063, by = 2), calc_cor_rel_mic, ps3, "M")
#mantel_subset_df <- rbind(out_f, out_m) %>% 
#                mutate(sex = c(rep("F", nrow(.)/2 ), rep("M", nrow(.)/2)))

load("output/mantel_subset_test_final.RData")
#save(mantel_subset_df, file = "output/mantel_subset_test_final.RData")
p_rel_sub <- ggplot(mantel_subset_df, aes(topX, mantelr, fill = sex, shape = sex, color = sex, by = sex)) +
    geom_pointrange(aes(ymin = llim.2.5., ymax = ulim.97.5.), fatten = 15, size = 0.1, alpha = 0.8) +
    theme_martin() +
    ylab("Mantel r\n") +
    xlab("\nNumber of most abundant bacterial taxa") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_color_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    scale_x_continuous(breaks = c(2, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1064)) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(colour = "black", size = 0.3),
          axis.text = element_text(colour = 'black'))
p_rel_sub
ggsave("../figures/Fig6_relatedness_subset.jpg", p_rel_sub, width = 4.8, height = 3)



# create dataframe for each timepoint per sex and combine for plotting

get_dist_df2 <- function(sub_timepoint, sub_sex, ps){
    # subset Males
    #  ps_merged_vst_sub  <- subset_samples(ps, (sex == sub_sex))
    ps_merged_vst_sub <- prune_samples((sample_data(ps)$sex == sub_sex) & 
                                       (sample_data(ps)$timepoint == sub_timepoint) , ps)
    # calculate bray curtis dissimilarity
    nes_mic_dist <- phyloseq::distance(ps_merged_vst_sub, method = "bray")
    labels(nes_mic_dist)
    # create df
    nes_mic_dist_df <- melt(as.matrix(nes_mic_dist), varnames = c("ind1", "ind2")) %>% 
        mutate(ind1 = str_replace(ind1, sub_timepoint, "")) %>% 
        mutate(ind2 = str_replace(ind2, sub_timepoint, "")) 
    # create joint df with microbial and genetic distance
    both_dist_df <- inner_join(nes_mic_dist_df, nes_rel_df, by = c("ind1" = "ind1", "ind2" = "ind2")) %>%
        mutate(sex = sub_sex) %>% 
        mutate(timepoint = sub_timepoint)
    # mutate(rel_classes = cut(rel, breaks = quantile(rel, probs = seq(0, 1, 0.15))))
}

# create data.frame
sex_distances <- map_df(c("F", "M"), function(sex) {
    map_df(c("T1", "T2", "T3"), get_dist_df2, sex, ps_vst)
    })
                   

p_rel_by_time <- ggplot(sex_distances, aes(rel, value, fill = sex, shape = sex)) +
    geom_point(size = 1.7, alpha = 0.3) +
    geom_smooth(se = FALSE, method = "lm", aes(color = sex)) +
    facet_grid(sex ~ timepoint) +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    scale_color_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    scale_x_continuous(breaks = c(-0.2, 0, 0.2)) +
    xlab("Genetic relatedness") +
    ylab("Bray-Curtis dissimilarity") +
    theme_martin() +
    theme( strip.background = element_blank(),
        strip.text.y = element_blank(), 
        panel.spacing = unit(1, "lines")) 

#ggsave("../figures/Sup7_relatedness_by_time.jpg", p_rel_by_time , width = 6.2, height = 4)


# mantel tests

all_mantels <- function(timepoint, sex, dist_df, nes_rel_df) {
    
    both_dist_df <- dplyr::filter(dist_df, timepoint == timepoint, sex == sex)
    
    #both_dist_df <- inner_join(nes_mic_dist_df, nes_rel_df, by = c("ind1" = "ind1", "ind2" = "ind2"))
    
    distmat_microbes <- as.dist(1 - create_distmat(both_dist_df, "ind1", "ind2", "value"))
    distmat_msats <- as.dist(create_distmat(both_dist_df, "ind1", "ind2", "rel"))
    
    mantel_test <- ecodist::mantel(formula =  distmat_microbes ~ distmat_msats, 
        mrank = FALSE, nperm = 10000, nboot = 10000)
    out <- data.frame(t(mantel_test), "timepoint" = timepoint, "sex" = sex)
}


all_mantels("T1", "M", sex_distances, nes_rel_df)




# inbreedR analyses
# calc het and g2
nes_geno <- convert_raw(nes_msats[3:ncol(nes_msats)])
g2_microsats(nes_geno, nboot = 1000, nperm = 1000)
# het df
nes_het <- data.frame(id = nes_msats$id, het = sMLH(nes_geno)) %>% 
  arrange(by = het)

# combine data
nes_full <- nes_divs %>% 
  mutate(id = str_replace(id, "17BEMa", "E")) %>% 
  right_join(nes_het)

ggplot(nes_full, aes(het, divs, col = sex)) + geom_point() + 
  facet_wrap(~timepoint) +
  geom_smooth(method = lm, se =FALSE)












≈