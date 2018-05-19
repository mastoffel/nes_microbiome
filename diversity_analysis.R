# Full analysis of ASVs for the paper

# load packages ------------------------------------------------------------------------------------
library(pacman)
p_load(phyloseq, tidyverse, msa, inbreedR, rptR, lme4, DESeq2, 
       dada2, phangorn, wesanderson, grid, cowplot, readxl, RColorBrewer, blogdown)
library(patchwork)
source("martin.R")
#library(microbiomeSeq)

# Prepare and load data ----------------------------------------------------------------------------

# input folder
input_folder <- "primer_clipped_reads_22_220230_pool"

# load taxa and RSV table
load(paste0("output/", input_folder, "/taxa.RData"))
load(paste0("output/", input_folder, "/seqtab.RData"))
load(paste0("output/", input_folder, "/fitGTR.RData"))

# sample names
samples_out <- rownames(seqtab_nochim)
nes_df <- data.frame("id" = samples_out)

# nes data
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

# rownames have to match SVG table
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


# plotting abundances ------------------------------------------------------------------------------

# plot relative abundances for all samples
ps_abund <- transform_sample_counts(ps0, function(x) x / sum(x))
plot_bar(ps_abund, fill="Phylum", x = "id")

top20 <- names(sort(taxa_sums(ps_abund), decreasing = TRUE))[1:50]
# ps_top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps_top20 <- prune_taxa(top20, ps)
plot_bar(ps_abund, fill="Class", x = "id")

ps_abund <- transform_sample_counts(ps, function(x) x / sum(x))
topXX <- names(sort(taxa_sums(ps_abund), decreasing = TRUE))[1:50]
ps_top <- prune_taxa(topXX, ps_abund)
plot_bar(ps_abund, fill="Class", x = "id", facet_grid = "sex") 


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


# Define prevalence of each taxa
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
ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = min_reads, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.6) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum) +
  ggtitle("Abundance by Phylum")


# plot Abundances --------------------------------------------------------
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

p_ord_df <- plot_ordination(ps_vst, ps_ord, shape = "sex", color = "timepoint", justDF = TRUE)

#which((p_ord_df$Axis.1 > 0) & (p_ord_df$Axis.2 < 0) & (p_ord_df$sex == "F"))
#p_ord_df[32, ]
p_ord_plot <- ggplot(p_ord_df, aes(Axis.1, Axis.2)) +
    geom_point(size = 3, alpha = 0.8, aes(shape = sex, fill = timepoint)) +
    #geom_point(size = 3, alpha = 0.8, aes( color = individual)) +
    theme_martin() +
    theme(panel.grid = element_blank()) +
    scale_shape_manual(values = c(21,24), name = "Sex")+
    scale_fill_manual(values = colpal, name = "Timepoint") +
    coord_fixed(sqrt(evals[2] / evals[1])) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal") +
    xlab("Axis 1 [28,5%]") +
    ylab("Axis 2 [13,4%]") +
    guides(fill=guide_legend(override.aes=list(shape=21)))
    
p_ord_plot
#ggsave("../figures/sex_time_MDS.jpg", p_ord_plot, width = 6, height = 4)

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
        geom_point(size = 1, alpha = 0.8) +
        #geom_point(size = 3, alpha = 0.8, aes( color = individual)) +
        theme_martin() +
        geom_polygon(aes(fill=individual), alpha = 0.05) +
        theme(panel.grid = element_blank(),
            legend.text=element_text(size=7),
            legend.key.size = unit(0.6,"line"),
            legend.title = element_blank(),
            plot.margin = unit(c(-1.1,0.3,-1.1,1), "cm")) +
        scale_x_continuous(breaks = seq(from = -0.45, to = 0.45, by = 0.2), limits = c(-0.45, 0.45)) +
        scale_y_continuous(breaks = seq(from = -0.4, to = 0.4, by = 0.2), limits = c(-0.4, 0.4)) +
        #scale_shape_manual(values = c(21,2))+
        scale_color_manual(values = plotcols) +
        coord_fixed(sqrt(evals[2] / evals[1])) +
       # theme(legend.position = "bottom",
        #    legend.direction = "horizontal") +
        xlab("Axis 1 [28,4%]") +
        ylab("Axis 2 [13,3%]") 
    p_ord_plot_ind
}
all_samples <- unique(p_ord_df$individual)
set.seed(22)
plotcols <- sample(c(get_colors("Paired")))
all_plots <- map(1:4, make_subset_plots, p_ord_df, plotcols)

p_all <- all_plots[[1]] + all_plots[[2]] + all_plots[[3]]  + all_plots[[4]]
p_all
#ggsave("../figures/Sup1_mds_by_int.jpg", p_all, width = 10, height = 5)


# plotting time trends -----------------------------------------------------------------------------
ps_rel <- transform_sample_counts(ps3, function(x) x/sum(x))
ps_df <- psmelt(ps_rel)
ps_df <- subset(ps_df, Abundance > 0)
#ps_df <- as_tibble(ps_df)

length(unique(ps_df$OTU))


# time trends across Classes ----------------------------------------------------------
plot_df <- ps_df %>% 
    group_by(Sample, Class) %>% 
    summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% 
    filter(!is.na(timepoint)) %>% 
    filter(Class %in% c("Actinobacteria", "Bacilli","Bacteroidia", "Betaproteobacteria",
                        "Clostridia", "Coriobacteriia", "Deferribacteres", "Deltaproteobacteria",
                        "Erysipelotrichia", "Flavobacteriia", "Fusobacteriia", "Gammaproteobacteria",
                        "Mollicutes",  "Negativicutes", "Spirochaetes", "Epsilonproteobacteria")) %>% 
    mutate(Abundance = log(Abundance + 0.001)) 

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

ggsave("../figures/classes_time_trends.jpg", p_class, width = 7.5, height = 6.5)



# time trends across Phyla ------------------------------------
plot_df <- ps_df %>% 
    group_by(Sample, Phylum) %>% 
    summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% 
    filter(!is.na(timepoint)) %>% 
    filter(Phylum %in% c("Actinobacteria", "Bacteriodetes", "Deferribacteres", "Firmicutes", "Fusobacteria",
                        "Proteobacteria", "Spirochaetae", "Tenericutes")) %>% 
    mutate(Abundance = log(Abundance + 0.001)) 

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

ggsave("../figures/phylum_time_trends.jpg", p_phylum, width = 7.5, height = 6.5)

# time trends across Order--------------------------------------------------------------------
plot_df <- ps_df %>% 
    group_by(Sample, Order) %>% 
    summarise_all(funs(if(is.numeric(.)) sum(., na.rm = TRUE) else first(.))) %>% 
    filter(!is.na(timepoint)) %>% 
    filter(Order %in% c("Actinomycetales", "Aeromonodales", "Anaeroplasmatales",
        "Bacillales", "Bacteroidales", "Burkholderiales", "Campylobacterales",
        "Clostridiales", "Coriobacteriales", "Corynebacteriales", "Desulfovibrionales",
        "Enterobacteriales", "Erysipelotrichales", "Flavobacteriales", "Fusobacteriales",
        "Lactobacillales", "Mycoplasmatales", "Neisseriales", "Pasteurellales", 
        "Pseudomonadales", "Rhodospirillaleles", "Selenomonadales", "Spirochaetales")) %>% 
    mutate(Abundance = log(Abundance + 0.001)) 

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
        shape = guide_legend())

p_order

ggsave("../figures/order_time_trends.jpg",p_order, width = 7.5, height = 6.5)


# alpha diversity ----------------------------------------------------------------------------------
diversity_df <- estimate_richness(ps, measures = c("Shannon", "Simpson", "InvSimpson", "Observed", "Fisher")) %>% 
                    tibble::rownames_to_column("id") %>% 
                    mutate(id = str_replace(id, "X", "")) %>% 
                    left_join(sample_data(ps), by = "id")

colpal_cavalanti <- wes_palette("Cavalcanti1", 2, type = "discrete")
as.character(wes_palette("Darjeeling2"))
colpal_moonrise <- c("#899DA4",  "#79402E")

#library(tidyr)
#div_lf <- gather(diversity_df, div_meas, div, c(Observed, Shannon, Simpson, InvSimpson))
set.seed(12)
p_div <- ggplot(diversity_df, aes(sex, Shannon)) + #colour = sex
    geom_boxplot(alpha = 1, outlier.shape = NA) + #, aes(color = sex)
    geom_jitter(size = 2.7, alpha = 0.6,  col = "black", aes(shape = sex, fill = "grey"), width = 0.3, stroke =0.7) + #shape = 21,
    #scale_shape_manual(values = c(21,24))+
    facet_wrap(~timepoint) +
    theme_martin() +
    ggtitle("Timepoint") +
    #scale_fill_manual(values =  "grey") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    #scale_x_discrete(labels = c(stri_unescape_unicode("\\u2640"), stri_unescape_unicode("\\u2642"))) +
    #scale_fill_manual(values = colpal_moonrise, name = "Sex") +
    #scale_color_manual(values = colpal_moonrise, name = "Sex") +
    ylab("Shannon diversity")+
    xlab("Sex")+
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
          #axis.text.x = element_text(family = "sans serif")) +
    guides(fill=FALSE, shape = FALSE) 
    #scale_x_discrete(labels = stri_unescape_unicode("a\\u0105!\\u0032\\n")) 
p_div
ggsave("../figures/diversity.jpg", p_div, width = 5.5, height = 3.4)

# modeling 
library(lme4)
library(partR2)
div_mod <- lmer(Shannon ~ sex + timepoint   + (1|individual), data = diversity_df)
summary(div_mod)
tidy(div_mod)
VarCorr(div_mod)
set.seed(17) 
boot_div_mod <- confint(div_mod, method = c("boot"), nsim = 1000)
boot_div_mod
R2_div <- partGaussian(div_mod, partvars = c("sex", "timepoint"), nboot = 1000)
rpt_div <- rptGaussian(Shannon ~ sex + timepoint + (1|individual),grname = "individual", data = diversity_df)


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
mod_ind <- adonis(phyloseq::distance(ps_vst, method="bray") ~ timepoint + sex + individual,
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
all_sigtabs <- map_df(all_contrasts,get_dds_results_tables, vst_dds, 0.01)

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










# gut microbiota and relatedness -------------------------------------------------------------------
library(Demerelate)
library(vegan)
library(ecodist)
library(reshape2)
nes_msats <- readxl::read_xls("../data/raw/nes_msats_cleaned.xls") %>% 
    add_column(factor(rep("SB", nrow(.))), .after = 1)
names(nes_msats)[1:2] <- c("Sample-ID", "Population")

# wang estimator performs well for small sample sizes and polymorphic loci (see paper)
# test whether number of loci is enough
#Loci.test(as.data.frame(nes), bt=1000, ref.pop="NA", object=TRUE, value= "wang", file.output=TRUE)
nes_rel <- Demerelate(as.data.frame(nes_msats), value= "morans", file.output=FALSE,  #loiselle
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

#as.dist(nes_rel_df)

# Calculate microbial distances

calc_cor_rel_mic <- function(topx, ps3) {
    
# find most abundant taxa
ps_rel <- transform_sample_counts(ps3, function(x) x/sum(x))
topXX <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:topx]
ps_topXXX <- prune_taxa(topXX, ps3)


# merge samples coming from one individual
#ps_merged <- subset_samples(ps3, sex == "M")
ps_merged <- merge_samples(ps_topXXX, "individual", fun = mean)
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
ps_merged_vst_sub  <- subset_samples(ps_merged_vst, (sex == "M" ))
# calculate bray curtis dissimilarity
nes_mic_dist <- phyloseq::distance(ps_merged_vst_sub, method = "bray")
labels(nes_mic_dist)


nes_mic_dist_df <- melt(as.matrix(nes_mic_dist), varnames = c("ind1", "ind2")) %>% 
                        mutate(ind1 = str_replace(ind1, "T3", "")) %>% 
                        mutate(ind2 = str_replace(ind2, "T3", "")) 


both_dist_df <- inner_join(nes_mic_dist_df, nes_rel_df, by = c("ind1" = "ind1", "ind2" = "ind2")) %>% 
                    mutate(rel_classes = cut(rel, breaks = quantile(rel, probs = seq(0, 1, 0.15))))

#ggplot(both_dist_df, aes(rel, value)) + 
#            geom_point(size = 3, alpha = 0.5) +
#            geom_smooth(method = "lm")

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

distmat_microbes <- as.dist(1 - create_distmat(both_dist_df, "ind1", "ind2", "value"))
distmat_msats <- as.dist(create_distmat(both_dist_df, "ind1", "ind2", "rel"))

mantel_test <- ecodist::mantel(formula =  distmat_microbes ~ distmat_msats, 
                               mrank = FALSE, pboot = 0.7, nperm = 10000, nboot = 10000)
out <- data.frame(t(mantel_test), "topX" = topx)
}


out <- map_df(seq(from = 5, to = 100, by = 5), calc_cor_rel_mic, ps3)


plot(out$mantelr, out$topx)




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