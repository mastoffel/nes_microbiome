# Complete data analysis and plotting for the paper

# Packages ----------------------------------------------------------------
library(pacman)
p_load(phyloseq, tidyverse, msa, inbreedR, rptR, lme4, DESeq2, 
       dada2, phangorn, wesanderson, grid, cowplot, readxl, RColorBrewer, 
       blogdown, patchwork, Demerelate, vegan, ecodist, reshape2, microbiome,
       kableExtra, broom, simpleboot, egg, insight, patchwork)
# package to calculate R2 
# devtools::install_github("mastoffel/partR2")
library(partR2)
# theme for plotting
source("theme_simple.R")
# useful function for microbiome data
source("microbiome_composition_funs.R")

# data
# load (mostly) unfiltered phyloseq object ps0
load("data/ps0.RData")
# load abundance and prevalence filtered phyloseq object ps3
load("data/ps3.RData")

# Variance stabilizing transformation ==========================================

# convert to deseq
nes_dds <- phyloseq_to_deseq2(ps3, ~1)
# estimate size factors with not including 0 in geometric mean calc
nes_dds  <- estimateSizeFactors(nes_dds, type = "poscounts") %>% 
  estimateDispersions(fitType = "local")
# create new phyloseq object with variance stabilised ASV table
ps_vst <- ps3
otu_table(ps_vst) <- otu_table(getVarianceStabilizedData(nes_dds), 
                               taxa_are_rows = TRUE)


# Microbiome composition =======================================================

# select phylum, transform to relativ abundancens and filter above 1% for plotting
nes_phylum <- ps3 %>% 
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>% 
  psmelt() %>% 
  filter(Abundance > 0.01) %>% 
  arrange(desc(Abundance))

# Figure 1, barplot
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
#ggsave(filename = "../figures/composition_phylum_bar.jpg", p_bar, width = 8, height = 4.5)


# Calculate summary statistics for ASVs ========================================

ps_otu <- as.data.frame(otu_table(ps3))
asv_per_sample <- rowSums(ps_otu > 0)
mean(asv_per_sample) # 286
sd(asv_per_sample) # 67


# Define new variable detect samples with low coverage =========================

# mark low abundance samples
sample_data(ps3)$abundance <- sample_sums(ps3) < 10000

# Core microbiome across individuals within each timepoint =====================
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
#sapply(all_timepoints, function(x) write_excel_csv(core_across_time[[x]], 
#                                   path = paste0("../data/processed/core_microbiome_", x, ".txt")))


# Core microbiome tables =======================================================
linesep <- c("") 

# core microbiome table in Supplementary
# uncomment lines for actually saving tables
for (i in all_timepoints) {
    core_timepoint <- dplyr::select(core_across_time[[i]], -ASV)
    kable(core_timepoint, format = "latex", booktabs = TRUE, escape = T, 
          linesep = linesep, align = "l") %>% 
          kable_styling(latex_options = c( "scale_down")) %>% 
          row_spec(0, bold = T) # %>% 
          # this saves the tables
          #kable_as_image(paste0("../tables/core_microbiome_", i), keep_pdf = TRUE, 
          #             file_format = "jpeg")
}

# Core microbiome plots (not in paper) -----------------------------------------

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
        coord_flip() +
        ggtitle("Core ASVs at T1")
#ggsave(filename = "../figures/core_t1.jpg", width = 4, height = 4)
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
    coord_flip() +
    ggtitle("Core ASVs at T2")
#ggsave(filename = "../figures/core_t2.jpg", width = 4, height = 4)

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
    coord_flip() +
    ggtitle("Core ASVs at T3")
#ggsave(filename = "../figures/core_t3.png", width = 5, height = 5)
# plot
p1 + p2 + p3

# Transformation and ordination ================================================

# which transformation is appropriate?
qplot(log(rowSums(otu_table(ps3)))) +
    xlab("Logged counts-per-sample")

# Figure 2: MDS plots ==========================================================

# Figure 2: MDS sex / age ======
# for plotting
get_colors <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
colpal <- c(  "#0000ff", "#ffb14e", "#ea5f94"  )
colpal <- wes_palette("Moonrise2", 3, type = "discrete")   

# multidimensional scaling based on bray-curtis
ps_ord <- ordinate(ps_vst, "MDS", "bray")
# calculate axis length relationships according to eigenvalues
evals <- ps_ord$values$Eigenvalues
# extract data to plot
p_ord_df <- plot_ordination(ps_vst, ps_ord, shape = "sex", color = "timepoint", 
                            justDF = TRUE, axes = 1:5)


# plot
p_ord_plot <- ggplot(p_ord_df, aes(Axis.1, Axis.2)) +
    geom_point(size = 3.5, alpha = 0.8, aes(shape = sex, fill = timepoint)) +
    #geom_point(size = 3, alpha = 0.8, aes( color = individual)) +
    #theme_martin() +
    theme_classic(base_size = 14) +
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
        axis.ticks = element_line(colour = "black", size = 0.3),
        legend.key.width = unit(0,  unit = "cm"),
        legend.spacing.x = unit(0.17, unit = "cm"),
        axis.text = element_text(color = "black")) +
    guides(fill=guide_legend(override.aes=list(shape=21)))
    
p_ord_plot
# ggsave("../figures/Fig1_sex_time_MDS_new.jpg", p_ord_plot, width = 6, height = 4)


# Figure 2: MDS health ======
# multidimensional scaling based on bray-curtis
ps_ord <- ordinate(ps_vst, "MDS", "bray")
evals <- ps_ord$values$Eigenvalues
p_ord_df <- plot_ordination(ps_vst, ps_ord, justDF = TRUE, axes = 1:10) %>% 
            mutate(health_status = as.numeric(health_status)) %>% 
            mutate(health_status = ifelse(is.na(health_status), "not_sampled", health_status))
          
# p_ord_T13 <- plot_ordination(ps_vst_T13, ps_ord, shape = "sex", color = "health_status",
#                             justDF = TRUE, axes = 1:4) %>% 
#   mutate(health_status = as.factor(as.numeric(health_status)))
data_sub <- filter(p_ord_df, timepoint %in% c("T1", "T3")) %>% 
            filter(category != "NA")

#plotcols <- c("#8da0cb", "#fc8d62", "grey90")
plotcols <- c(wes_palette("Royal1")[2],"#39312F", "grey90")
p_health_plot <- ggplot(p_ord_df, aes(Axis.1, Axis.2, shape = sex)) +
  geom_point(size = 3.5, alpha = 0.6, fill = "grey90", stroke = 0.2) +
  geom_point(data = data_sub, 
             aes(fill=health_status),
             size = 3.5, alpha = 0.8, stroke = 0.5) +
  scale_shape_manual(values = c(21,24), name = "Sex") +
  scale_x_continuous(breaks = seq(from = -0.25, to = 0.25, by = 0.25)) + # limits = c(-0.45, 0.45))
  scale_y_continuous(breaks = seq(from = -0.2, to = 0.2, by = 0.2)) + #, limits = c(-0.4, 0.4)
  scale_fill_manual("Health status", values = c(plotcols), labels = c("Clinically\nabnormal", "Clinically\nhealthy")) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  xlab("Axis 1 [28,5%]") +
  ylab("Axis 2 [13,4%]") +
  guides(shape = FALSE,
         fill = guide_legend(override.aes = list(shape = c(22,22))),
         keywidth = 1) +
  theme_classic(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3, linetype = 1),
        axis.line.y = element_line(colour = "black", size = 0.3, linetype = 1),
        axis.ticks = element_line(colour = "black", size = 0.3),
        axis.text = element_text(color = "black"),
        legend.position = "bottom") 
p_health_plot

# Figure 2: MDS host ID ======
set.seed(25)
plotcols <- sample(c(wes_palette("FantasticFox1")[-4], 
                     wes_palette("GrandBudapest2")[-1], wes_palette("Moonrise2")), 6)
# good inds to show
# 17BEMa32, 17BEMa27 or 24, 17BEMa38, 17BEMa3
# random number sampler changed across R versions, this keeps it the old way
RNGkind(sample.kind = "Rounding")
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
  theme_classic(base_size = 14) +
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
        axis.text = element_text(color = "black"),
        legend.position = "bottom") 
 
#p_full <- plot_grid(p_ord_plot, p_host, p_health_plot, labels = c("A", "B", "C"), ncol = 2)

p_full <- (p_ord_plot + p_host) / (p_health_plot + plot_spacer()) +
  plot_annotation(tag_levels = 'A') 
ggsave("../figures/beta_div3.jpg",plot = p_full, width = 9, height = 7)
ggsave("../figures/beta_div3.pdf",plot = p_full, width = 9, height = 7)
# Ordination outlier plot ======================================================

# This plot is to check whether samples with low read abundances
# are large outliers.
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

# plotting time trends as supplementary figures ================================
# relative abundances
ps_rel <- transform_sample_counts(ps3, function(x) x/sum(x))
ps_df <- psmelt(ps_rel)
ps_df <- subset(ps_df, Abundance > 0)
length(unique(ps_df$OTU))

# Time trends across classes (Supplementary Figure 2) ==========================
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
    theme_minimal() +
    ylab("Relative abundance") +
    xlab("Timepoint") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    guides(fill = guide_legend(),
           shape = guide_legend())

p_class

#ggsave("../figures/Fig2_classes_time_trends.jpg", p_class, width = 7.5, height = 6.5)
#ggsave("../figures/Sup3_classes_time_trends_full.jpg", p_class, width = 9.7, height = 7.5)


# Time trends across phyla (Supplementary Figure 1) ============================

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
    theme_minimal() +
    ylab("Relative abundance") +
    xlab("Timepoint") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    guides(fill = guide_legend(),
        shape = guide_legend())

p_phylum
# ggsave("../figures/Sup4_phylum_time_trends.jpg", p_phylum, width = 7.5, height = 6.5)

# Time trends across orders (Supplementary Figure 3) ===========================
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
    theme_minimal() +
    ylab("Relative abundance") +
    xlab("Timepoint") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    guides(fill = guide_legend(),
        shape = guide_legend()) +
    theme(strip.text = element_text(size=9),
          axis.text = element_text(size = 9))

p_order
# ggsave("../figures/Sup5_order_time_trends.jpg",p_order, width = 8.5, height = 7.5)


# Time trends across genera (not in supplementary)   ===========================

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
    theme_minimal() +
    ylab("Relative abundance") +
    xlab("Timepoint") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    guides(fill = guide_legend(),
        shape = guide_legend())

p_genus
# ggsave("../figures/Fig2_core_microbiota_genus_time_trends.jpg", p_genus, width = 7.5, height = 5.5)

# Alpha diversity across time and sex   =============================
diversity_df <- estimate_richness(ps0, measures = c("Shannon", "Simpson", "InvSimpson", "Observed", "Fisher")) %>% 
                    tibble::rownames_to_column("id") %>% 
                    mutate(id = str_replace(id, "X", "")) %>% 
                    left_join(as_tibble(sample_data(ps0)), by = "id")

colpal_cavalanti <- wes_palette("Cavalcanti1", 2, type = "discrete")
as.character(wes_palette("Darjeeling2"))
colpal_moonrise <- c("#899DA4",  "#79402E")

set.seed(12)
p_div <- ggplot(diversity_df, aes(timepoint, Shannon, by = sex)) + #colour = sex
    geom_boxplot(alpha = 0.4, outlier.shape = NA, aes(fill = sex)) + #, aes(color = sex)
    geom_point(position=position_jitterdodge(jitter.width = 0.3), size = 2.3, alpha = 0.8,  
        col = "black", aes(shape = sex, fill = sex), stroke =0.7) +
    theme_martin(base_family = "Helvetica", highlight_family = "Helvetica") +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    ylab("Shannon diversity\n")+
    xlab("\nTimepoint")+
    ylim(1.8, 5) +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.line = element_line(colour = "black", size = 0.3, linetype = 1),
          axis.ticks = element_line(colour = "black", size = 0.3, linetype = 1),
          plot.margin = margin(t = 25)) # plot.margin = margin(t = 25)
    #scale_x_discrete(labels = stri_unescape_unicode("a\\u0105!\\u0032\\n")) 
p_div
# ggsave("../figures/Fig4_diversity.jpg", p_div, width = 4.5, height = 2.9)

# Alpha diversity for health and non-healthy individuals (Figure 1B)  ==========
plotcols <- c("#39312F", wes_palette("Royal1")[2])
diversity_df2 <- diversity_df %>% filter(!is.na(health_status)) %>%  # "Shannon", "Simpson", "InvSimpson", "Observed", "Fisher"
                  mutate(health_status = ifelse(health_status == 1, 0, 1)) %>% 
                  mutate(health_status = as.factor(health_status))


my_tag <- c("T1", "T3")
p_div2 <- ggplot(diversity_df2, aes(sex, Shannon, by = health_status)) + #colour = sex
  geom_boxplot(alpha = 0.4, outlier.shape = NA, aes(fill = health_status)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.1), size = 2.3, alpha = 0.8,  
             col = "black", aes(shape = sex, fill = health_status), stroke =0.7) +
  facet_wrap(~timepoint) +
  scale_fill_manual("Health status", values = c(plotcols), labels = c( "Clinically healthy", "Clinically abnormal")) +
  scale_shape_manual(values = c(21,24), name = "Sex", guide = FALSE) +
  theme_martin(base_family = "Helvetica", highlight_family = "Helvetica") +
  ylab("Shannon diversity\n")+
  xlab("\nSex")+
  ylim(1.8, 5) +
  scale_x_discrete(labels = c("Female", "Male")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        panel.grid = element_blank(),
       # strip.text = element_text(size = 11), #face="bold", 
        strip.text = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black", size = 0.3, linetype = 1),
        axis.ticks = element_line(colour = "black", size = 0.3, linetype = 1),
        plot.margin = margin(t = 25)) 
p_div2_final <- tag_facet(p_div2, tag_pool = my_tag, open="",close="", x = 1.5, y=5, fontface = 1)

# Alpha diversity Figure 4 =====================================================
p_div_full <- plot_grid(p_div, p_div2_final, labels = c("A", "B"), ncol = 2, 
                        rel_widths = c(0.42,0.58), scale = 0.95)
p_div_full
ggsave("../figures/Fig4_diversity2.jpg", p_div_full, width = 8.5, height = 3.1)


# Alpha diversity stats ========================================================

####### check partR2 package here ###########
# sex, time point and host ID (individual)
div_mod <- lmer(Shannon ~ sex + timepoint + (1|individual), data = diversity_df)
summary(div_mod)
tidy(div_mod)
#VarCorr(div_mod)
set.seed(17) 
# parametric bootstraps
boot_div_mod <- confint(div_mod, method = c("boot"), nsim = 1000)
boot_div_mod
# calculate R2 with convidence intervals
R2_div <- partR2(div_mod, partvars = c("sex", "timepoint"), nboot = 1000) 
# calculate repeatability
rpt_div <- rptR::rptGaussian(Shannon ~ (1|individual), grname = "individual", 
                       data = diversity_df)

# health analysis
diversity_df %>% 
  filter(!is.na(health_status)) %>% 
  group_by(sex, timepoint) %>% 
  summarise(prop_healthy = mean(as.numeric(as.character(health_status))))
# subset timepoints 1 and 3 for which we measured health status
div_t1t3 <- diversity_df %>% 
  filter(timepoint %in% c("T1", "T3"))
# model
div_mod2 <- lmer(Shannon ~ sex + health_status + timepoint + (1|individual), data = div_t1t3)
summary(div_mod2)
plot(div_mod2)
boot_div_mod2 <- confint(div_mod2, method = c("boot"), nsim = 1000)

R2 <- function(mod) {
  var_comps <- insight::get_variance(mod)
  R2 <- var_comps$var.fixed / (var_comps$var.fixed + var_comps$var.residual + 0) # 0 variance in random effect
}

R2s <- bootMer(div_mod2, R2, nsim = 1000, type = "parametric")
# CI
quantile(R2s$t, probs = c(0.025, 0.975))
#partR2(div_mod2, partvars = c("sex", "timepoint", "health_status"))


# Beta diversity stats =========================================================

# all permanova models based on variance-stabilised data

# overall effects 
metadata <- as(sample_data(ps_vst), "data.frame")
mod_full <- adonis2(phyloseq::distance(ps_vst, method="bray") ~ timepoint + sex + individual,
    data = metadata, by = "terms", strata = "individual") # "terms"
mod_full

# sex effect at T1, T2 and T3 (check that repeated samples didn't effect sex effect) ======
#
sex_eff <- function(tpoint) {
  ps_vst_t <- prune_samples( sample_data(ps_vst)$timepoint == tpoint, ps_vst)
  metadata_t <- as(sample_data(ps_vst_t), "data.frame")
  mod_t <- adonis(phyloseq::distance(ps_vst_t, method="bray") ~ sex, data = metadata_t)
}
all_sex_eff <- purrr::map(c("T1", "T2", "T3"), sex_eff)

# t1 vs t2
ps_vst_t1t2 <- subset_samples(ps_vst, timepoint %in% c("T1", "T2"))
metadata_t1t2 <- as(sample_data(ps_vst_t1t2), "data.frame")
mod_t1t2 <- adonis(phyloseq::distance(ps_vst_t1t2, method="bray") ~ timepoint + sex + individual,
    data = metadata_t1t2, by = "terms")

# t2 vs t3
ps_vst_t2t3 <- subset_samples(ps_vst, timepoint %in% c("T2", "T3"))
metadata_t2t3 <- as(sample_data(ps_vst_t2t3), "data.frame")
mod_t2t3 <- adonis(phyloseq::distance(ps_vst_t2t3, method="bray") ~ timepoint + sex + individual,
    data = metadata_t2t3, by = "terms")

# individuals, repeatability 
metadata <- as(sample_data(ps_vst), "data.frame")
mod_rpt <- adonis(phyloseq::distance(ps_vst, method="bray") ~  sex + individual, 
    data = metadata)

# Are we seeing differences in group means of their dispersion? To check this,
# we use betadisper (looks like there is no difference in dispersion)
mod <- betadisper(d = phyloseq::distance(ps_vst, method="bray"), group = metadata$sex)
anova(mod) # parametric anova
TukeyHSD(mod) # Tukey Honest Significant Difference Method
permutest(mod) # Permutation test
mod <- betadisper(d = phyloseq::distance(ps_vst, method="bray"), group = metadata$timepoint)
anova(mod)
TukeyHSD(mod)
permutest(mod)

# plotting all individuals
sample_data(ps_vst)$individual <- as.factor(sample_data(ps_vst)$individual)
p_ord_rpt <- plot_ordination(ps_vst, ps_ord, color = "individual", type = "samples")
p_ord_rpt  + geom_polygon(aes(fill=individual), alpha = 0.05) + geom_point(size=1) + ggtitle("samples") +
    theme_martin()

# health analysis (only T1 and T3)
ps_vst_health <- subset_samples(ps_vst, !is.na(health_status))
metadata <- as(sample_data(ps_vst_health), "data.frame") 
mod_full <- adonis2(phyloseq::distance(ps_vst_health, method="bray") ~ health_status + sex + timepoint + individual,
                   data = metadata, by = "terms", strata = "individual")
mod_full

# health effect at T1, T3 (analysis without repeated measures) =================
ps_vst_health <- subset_samples(ps_vst, !is.na(health_status) & timepoint == "T1")
metadata <- as(sample_data(ps_vst_health), "data.frame") 
mod_health_t1 <- adonis2(phyloseq::distance(ps_vst_health, method="bray") ~ health_status + sex,
                    data = metadata, by = "terms")
mod_health_t1

ps_vst_health <- subset_samples(ps_vst, !is.na(health_status) & timepoint == "T3")
metadata <- as(sample_data(ps_vst_health), "data.frame") 
mod_health_t3 <- adonis2(phyloseq::distance(ps_vst_health, method="bray") ~ health_status + sex,
                         data = metadata, by = "terms")
mod_health_t3


# check again whether it's differences in means or dispersion
mod <- betadisper(d = phyloseq::distance(ps_vst_health, method="bray"), group = metadata$health_status)
anova(mod) # parametric anova
TukeyHSD(mod) # Tukey Honest Significant Difference Method
permutest(mod) # Permutation test


# ///experimental/// repeatability for distance matrices -----------------------------
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



# Which ASVs are differentially abundant among sexes and across time? ==========
# We are using deseq2 to test that.

# Differential abundance modeling (1): Across time, per sex ====================

# function runs deseq seperately for males and females
# to explore differential ASVs across time
run_deseq_per_sex <- function(sex){
    # run two models for females and males seperately to account for repeated measures
    ps_m <- prune_samples(sample_data(ps3)$sex %in% sex, ps3)
    # create deseq object
    nes_mod_dds_m <- phyloseq_to_deseq2(ps_m, ~ individual + timepoint)
    # estimate size factos
    vst_dds_m <- estimateSizeFactors(nes_mod_dds_m , type = "poscounts")
    # deseq analysis    ===========
    vst_dds <- DESeq(vst_dds_m, test="Wald", fitType="local")
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
get_colors <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)

set.seed(2049) # 5 / 559 / 1053
plotcols <- sample(c(get_colors("Paired"), "grey", "black", "white"))

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
    theme_minimal() +
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

p_diff_time_final <- ggdraw(p_diff_time) + 
    draw_label("\u2642",size = 35, colour = "black", x = 0.476, y = 0.3,lineheight = 5) +
    draw_label("\u2640",size = 35, colour = "black", x = 0.476, y = 0.75,lineheight = 5) 
p_diff_time_final 
# ggsave("../figures/timepoint_diff_abund_family.jpg", p_diff_time_final, width = 11, height = 8)


# Summary statistics of differential abundance for paper =======================

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


# Differential abundance modeling (2): Across sexes, within timepoints =========

# create factor representing the combination of time and sex
sample_data(ps3)$group <- factor(paste0(sample_data(ps3)$timepoint, sample_data(ps3)$sex))
nes_mod_dds_sex <- phyloseq_to_deseq2(ps3, ~ group)

# estimate size factos
vst_dds_sex <- estimateSizeFactors(nes_mod_dds_sex, type = "poscounts")
colData(nes_mod_dds_sex)

# deseq analysis  
vst_dds_sex_mod <- DESeq(vst_dds_sex, test="Wald", fitType="local")
#colData(vst_dds_sex_mod)
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

# create labeller function
wrap_labels <- c(T1FvsT1M = "T1: F vs. M", T2FvsT2M = "T2: F vs. M", T3FvsT3M = "T3: F vs. M")
#levels(all_sigtabs_plot$Class)

# plot
p_diff_sex <- ggplot(all_sigtabs_plot, aes(x=Family, y=log2FoldChange)) + 
    #geom_point(size = 3) + 
    geom_point(size=3, alpha = 0.7, shape = 21, colour = "black", lwd = 0.1, aes(fill=Class)) + 
    theme_minimal() +
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
# ggsave("../figures/timepoint_diff_abund_sex.jpg", p_diff_sex, width = 11, height = 5)

# Summary statistics of differential abundance for paper -----------------------
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



# Microbiome and genetic relatedness ===========================================

# microsatellites data preparation
nes_msats <- readxl::read_xls("data/nes_msats_cleaned.xls") %>% 
    add_column(factor(rep("SB", nrow(.))), .after = 1)
names(nes_msats)[1:2] <- c("Sample-ID", "Population")
nes_msats <- data.frame(nes_msats)
str(nes_msats)

# test whether number of microsat loci is enough to estimate relatedness =======
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

# (1) relatedness estimation ===================================================
nes_rel <- Demerelate(as.data.frame(nes_msats), value= "loiselle", 
                      file.output=FALSE,  #loiselle
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

# convert to deseq
ps_merged_dds <- phyloseq_to_deseq2(ps_merged, ~ 1) # here ps_merged

# estimate size factors with not including 0 in geometric mean calc
ps_merged_dds <- estimateSizeFactors(ps_merged_dds, type = "poscounts") %>% 
    estimateDispersions(fitType = "local")

# create new phyloseq object with variance stabilised ASV table
ps_merged_vst <- ps3
otu_table(ps_merged_vst) <- otu_table(getVarianceStabilizedData(ps_merged_dds), 
                                      taxa_are_rows = TRUE)

# create dataframe for each sex and combine for plotting
get_dist_df <- function(sub_sex, ps){
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
# ggsave("../figures/Fig5_relatedness.jpg", p_rel, width = 7, height = 3)

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

# calculate mantel test seperately for the two sexes ---------------------------
# "M" for males, "F" for females
both_dist_df <- filter(sex_distances, sex == "M")

distmat_microbes <- as.dist(1 - create_distmat(both_dist_df, "ind1", "ind2", "value"))
distmat_msats <- as.dist(create_distmat(both_dist_df, "ind1", "ind2", "rel"))

mantel_test <- ecodist::mantel(formula =  distmat_microbes ~ distmat_msats, 
    mrank = FALSE, pboot = 0.9, nperm = 10000, nboot = 10000)

mantel_test 
# results
# M: mantelr = 0.2623, CI 0.05882922 0.41531428, p two tailed = 0.00090000
# F: mantelr = 0.0588, CI -0.0676119  0.2006059 , p two tailed = 0.4209000

# differences in significance are not necessarily significant differences!
# check whether slopes are different (not possible with mantel tests)

# idea: permutation test. Permute pairwise relatedness and recalculate interaction slope
set.seed(3110)
mod_rel <- lm(rel ~ value*sex, data = sex_distances)
boot_rel <- lm.boot(mod_rel, 10000, rows = FALSE)
perc(boot_rel , p = c(0.025, 0.975)) # interaction slope -0.11 CI: -[0.23, -0.006]
org_slope <- summary(lm(rel ~ value*sex, data = sex_distances))$coefficients[4, 1]

# function to recalculate slope using permuted data
permed_slope <- function(iter, sex_distances) {
  data_perm <- sex_distances
  data_perm[, "rel"] <- sex_distances$rel[sample(1:nrow(sex_distances))]
  out <- summary(lm(value ~ rel*sex, data = data_perm))$coefficients[4, 1]
  out
}

perm_slopes <- sapply(1:10000, permed_slope, sex_distances)
hist(perm_slopes)
# permutation based p-value, are the slopes different?
p_val <- 1 - sum(perm_slopes > org_slope) / length(perm_slopes)
p_val

# How many ASVs carry an association with genetic relatedness? =================

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

# this analysis takes time, so load pre-saved data
load("output/mantel_subset_test_final.RData")
#save(mantel_subset_df, file = "output/mantel_subset_test_final.RData")
p_rel_sub <- ggplot(mantel_subset_df, aes(topX, mantelr, fill = sex, shape = sex, color = sex, by = sex)) +
    geom_pointrange(aes(ymin = llim.2.5., ymax = ulim.97.5.), fatten = 15, size = 0.1, alpha = 0.8) +
    theme_minimal() +
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
# ggsave("../figures/Fig6_relatedness_subset.jpg", p_rel_sub, width = 4.8, height = 3)



# create dataframe for each timepoint per sex and combine for plotting ---------

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
                   
p_rel_by_time <- ggplot(sex_distances, aes(rel, 1-value, fill = sex, shape = sex)) +
    geom_point(size = 1.7, alpha = 0.3) +
    geom_smooth(se = FALSE, method = "lm", aes(color = sex)) +
    facet_grid(sex ~ timepoint) +
    scale_shape_manual(values = c(21,24), name = "Sex") +
    scale_fill_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    scale_color_manual(values = c("#046C9A", "#D69C4E"), name = "Sex") +
    scale_x_continuous(breaks = c(-0.2, 0, 0.2)) +
    xlab("Genetic relatedness") +
    ylab("Bray-Curtis similarity") +
    theme_martin2() +
    theme( strip.background = element_blank(),
        strip.text.y = element_blank(), 
        panel.spacing = unit(1, "lines")) 

ggsave("../figures/Sup7_relatedness_by_time.jpg", p_rel_by_time , width = 6.2, height = 4)

