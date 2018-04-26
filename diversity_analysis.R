# phyloseq
library(phyloseq)
library(ggplot2)
library(dplyr)
library(stringr)
library(msa)
library(inbreedR)
library(readxl)
library(stringr)
library(dplyr)
library(rptR)
library(plyr)
library(reshape2)
library(lme4)
library(DESeq2)
library(dada2)
library(phangorn)
library(purrr)
library(related)
library(forcats)
library(wesanderson)
source("martin.R")
# todo: filter fecal sample or let it in!


# Prepare and load data --------------------------------------------------------

# input folder
# input_folder <- "combined_reads_1_nopool"
input_folder <- "primer_clipped_reads_22_250250_pool"
# load taxa and RSV table
load(paste0("output/", input_folder, "/taxa.RData"))
load(paste0("output/", input_folder, "/seqtab.RData"))
load(paste0("output/", input_folder, "/fitGTR.RData"))
# sample names
samples_out <- rownames(seqtab_nochim)
nes_df <- data.frame("id" = samples_out)

# nes data
nes_sampling <- read_xlsx("../data/sampling_data_processed.xlsx") %>% 
  dplyr::rename(id = ID,
                date = DATE,
                territory = TERRITORY, 
                sex = SEX) %>% 
  mutate(id = str_replace(id, "17BEMA0", "17BEMa")) %>% 
  mutate(id = str_replace(id, "17BEMA", "17BEMa")) %>% 
  mutate(id = str_c(id, timepoint)) %>% 
  mutate(id = str_replace(id, "T1", "")) %>% 
  mutate(sex = ifelse(sex == "M", "F", "M")) %>% 
  dplyr::select(-date) %>% 
  as.data.frame()


# join to nes_df to get the right sample sequence
nes <- left_join(nes_df, nes_sampling, by = "id") %>% 
          arrange(timepoint) %>% 
          mutate(id = fct_inorder(id))
          
          
#which(nes$id == "17BEMa11Fec")
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
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE), 
               sample_data(nes2), 
               tax_table(taxa),
               phy_tree(fitGTR$tree))

# filter out fecal sample
ps <- subset_samples(ps, id != "17BEMa11Fec")

# plotting abundances
plot_bar(ps, fill = "Phylum", x = "id")
# plot top twenty
top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:50]
# ps_top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps_top20 <- prune_taxa(top20, ps)
plot_bar(ps_top20, fill="Phylum", x = "id") +
  scale_y_sqrt()

top20 <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:20]
# ps_top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps_top20 <- prune_taxa(top20, ps_rel)
plot_bar(ps_top20, fill="Order", x = "id") +
  ggtitle("20 most abundant ASVs across samples")
# overview over phyloseq options -----------------------------------------------
#ps
#ntaxa(ps)
#nsamples(ps)
#sample_names(ps)[1:5]
#rank_names(ps)
#sample_variables(ps)
# #otu_table(ps)[1:3, 1:3]
#tax_table(ps)[1:5, 1:5]
#taxa_names(ps)[1]

# prevalence filtering ---------------------------------------------------------
# which samples have low abundance?
sort(sample_sums(ps))
# ps <- prune_samples(sample_sums(ps) >= 600, ps)

# in the dataset, which we will define here as the number of samples 
# in which a taxon appears at least once.
# Compute prevalence of each feature, store as data.frame
prev0 <- apply(X = otu_table(ps),
                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prev0,
                     TotalAbundance = taxa_sums(ps),
                     tax_table(ps))

# checks
prevdf_check <- as_tibble(prevdf)
hist(prevdf_check$Prevalence, breaks = 100)

# check out how many phyla and whether there are NAs)
sort(table(prevdf$Phylum, useNA = "always"), decreasing = TRUE) # 11 unidentified Phyla --> filter

# checked that the abundant phyla were found in other species too
# cut off everything <11
sort(table(prevdf$Phylum))
keepPhyla <- table(prevdf$Phylum)[(table(prevdf$Phylum) > 5)]
prevdf1 <- subset(prevdf, Phylum %in% names(keepPhyla))

# Keep taxa when appearing in minimum 2% samples
prevalenceThreshold <- 0.02 * nsamples(ps)
prevalenceThreshold

# execute prevalence filter 
ps1 <- prune_taxa((prev0 > prevalenceThreshold), ps)
ps1

# filter entries with unidentified phylum
ps2 <- subset_taxa(ps1, Phylum %in% names(keepPhyla))

ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = Class)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.6) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum) +
  ggtitle("Abundance by Phylum")


# plot Abundances and potentially transform ------------------------------------
# Taxonomic agglomeration
# How many genera are present after filtering? 135
taxGlomRank <- "Genus"
length(get_taxa_unique(ps2, taxonomic.rank = taxGlomRank)) 

rank_names(ps2)
# Abundance value plotting
plot_abundance = function(physeq, ylabn = "",
                          Facet = "Order",
                          Color = "Phylum"){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq,
         mapping = aes_string(x = "sex", y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + ylab(ylabn) +
    scale_y_log10()
}
plot_abundance(ps2)

# define some extra variables --------------------------------------------------
# mark low abundance samples and impute birth territory
sample_data(ps2) <- sample_data(ps2) %>% 
  mutate(abundance = sample_sums(ps2) < 5000) 


# ordination -------------------------------------------------------------------
# Transform to relative abundance. Save as new object.
ps_rel <- transform_sample_counts(ps2, function(x){x / sum(x)})
ps_rel <- transform_sample_counts(ps2, function(x) log(x + 1))
ps_ord <- ordinate(ps_rel, "NMDS", "bray")
# plot taxa
# plot_ordination(ps_rel, ps_ord, type="taxa", color="Phylum", title="taxa")
# plot samples
colpal <- c(  "#0000ff", "#ffb14e", "#ea5f94"  )
colpal <- wes_palette("Moonrise2", 3, type = "discrete")      

p1 <- plot_ordination(ps_rel, ps_ord, type="samples", color="abundance", shape = "sex") 
p1 +
  theme_martin() +
  theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c(19,17))+
  geom_point(size = 3, alpha = 1) +
  scale_color_manual(values = colpal) +
  ggtitle("abundance < 5000 reads")

#p1 + geom_polygon(aes(fill=timepoint)) + geom_point(size=5) + ggtitle("samples")

plot_richness_estimates(ps, x = "timepoint", measure = c("Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot() +
  theme_martin() +
  theme(panel.grid = element_blank()) 

# try variance-stabilising transformation Now perform the variance-stabilizing transformation and replace the original OTU abundances in a copy of GP with them. We'll call the copy GPvst.
# Use `~1` as the experimental design so that the actual design doesn't
# influence your tranformation.
ps_sd <- phyloseq_to_deseq2(ps2, ~1)
ps_sd <- estimateSizeFactors(ps_sd)
ps_sd <- estimateDispersions(ps_sd, fitType = "local")

ps_vst <- ps2
otu_table(ps_vst) <- otu_table(getVarianceStabilizedData(ps_sd), taxa_are_rows = TRUE)
plot_ordination(ps_vst, ordinate(ps_vst, "NMDS", "bray"), color="timepoint", shape = "sex")+
  theme_martin() +
  theme(panel.grid = element_blank()) +
  scale_shape_manual(values = c(19,17))+
  geom_point(size = 3, alpha = 1) +
  scale_color_manual(values = colpal)

ps_rl <- ps2
rld <- rlog(ps_sd, blind = FALSE)
rownames(rld) <- taxa_names(ps2)
otu_table(ps_rl) <- otu_table(rld, taxa_are_rows = TRUE)
plot_ordination(rld, ordinate(Grld, "NMDS", "bray"), color = "SampleType")


# deseq2 analysis --------------------------------------------------------------
timedds <- phyloseq_to_deseq2(ps2, ~ timepoint)
timedds <- DESeq(timedds, test="Wald", fitType="parametric")
res <- results(timedds, cooksCutoff = FALSE)
alpha <- 0.01
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps2)[rownames(sigtab), ], "matrix"))
sigtab_check <- as_tibble(sigtab)

# which ASV were different?
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x <- sort(x, TRUE)
sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab$Genus <- factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  ggtitle("Sex differences")

# preprocessing 

# cumulative sum scaling -------------------------------------------------------
library("metagenomeSeq")
## Convert the phyloseq object to a metagenomeSeq object (MRexperiment).
## The Phyloseq_to_metagenomeSeq function is included in the phyloseq package.
metagenome.obj <- phyloseq_to_metagenomeSeq(ps2)
## Calculate the proper percentile by which to normalize counts
cNstat <- metagenomeSeq::cumNormStatFast(metagenome.obj) 
# cNstat  #0.5
## Normalise counts
metagenome.obj <- metagenomeSeq::cumNorm(metagenome.obj, p = cNstat)
## Export the normalised count table
metag.norm.counts <- metagenomeSeq::MRcounts(metagenome.obj, norm = TRUE)
## Add a pseudocount of 0.0001 to the table and log transform
metag.norm.counts_log <- log(metag.norm.counts+0.0001)
## Substract the value from log of pseudocount to preserve zeros of the original counts
metag.norm.counts_log2 <- metag.norm.counts_log-(log(0.0001))

pstest <- phyloseq(otu_table(metag.norm.counts_log2, taxa_are_rows = TRUE), 
               sample_data(nes), 
               tax_table(taxa),
               phy_tree(fitGTR$tree))

##  Calulate PCoA for CSS normalised OTU table with Bray-Curtis.
nes_ord <- ordinate(physeq = pstest, method = "PCoA", distance = "bray")
plot_ordination(pstest, nes_ord, color = "timepoint")

# other transformation ---------------------------------------------------------
# log 10 transformation looks good
qplot(log10(rowSums(otu_table(ps2)))) +
  xlab("Logged counts-per-sample")

# log transform
pslog <- transform_sample_counts(ps2, function(x) log(1 + x))
# Transform to relative abundance. Save as new object.
ps3 <- transform_sample_counts(ps2, function(x){x / sum(x)})

# unifrac MDS
out_wuf_log <- ordinate(pslog, method = "MDS", distance = "bray")
# 
evals <- out_wuf_log$values$Eigenvalues
# plot the ordination
plot_ordination(ps3, out_wuf_log, shape = "sex", color = "timepoint") +
  #labs(col = "sex") +
  geom_point(size = 3) +
  coord_fixed(sqrt(evals[2] / evals[1]))

# rel_abund <- t(apply(otu_table(ps2), 1, function(x) x / sum(x)))
# qplot(rel_abund[, 14], geom = "histogram") +
#   xlab("Relative abundance")

# bray curtis with MDS
out_bc_log <- ordinate(ps3, method = "MDS", distance = "bray")
evals <- out_bc_log$values$Eigenvalues
plot_ordination(pslog, out_bc_log, color = "timepoint") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "timepoint")

# PCoA with weighted unifrac
out_wuf_log <- ordinate(pslog, method = "PCoA", distance = "wunifrac")
evals <- out_wuf_log$values$Eigenvalues
# plot the ordination
plot_ordination(pslog, out_wuf_log, color = "timepoint") +
  labs(col = "timepoint") +
  coord_fixed(sqrt(evals[2] / evals[1]))


# PCA on ranks
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 15)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = .7) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")

# all those microbes with rank below some threshold are set to be tied
# at 1. The ranks for the other microbes are shifted down, so there is no large gap between ranks.
abund_ranks <- abund_ranks - 400
abund_ranks[abund_ranks < 1] <- 1

# statistical tests
library(phyloseqGraphTest)
# (1) Minimum spanning tree with Jaccard dissimilarity
gt <- graph_perm_test(ps3, sampletype = "timepoint",
                      #grouping = "timepoint",
                      distance = "bray", type = "knn", knn = 3)
plot_test_network(gt) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9)) 
plot_permutations(gt) 


# try rarefaction
hist(sample_sums(ps2), breaks = 112)
sort(sample_sums(ps2))

resample_ASV <- function(num_boot, physeq, sample_size = 1000) {
  ps_rf <- rarefy_even_depth(ps2, sample.size = sample_size, replace = TRUE, trimOTUs = TRUE)
  nes_richness <- estimate_richness(ps_rf, split = TRUE, measures = c("Shannon"))
}
# subsampling 1000 times
all_div <- map(1:1000, resample_ASV, ps2, sample_size = 3000)
all_div_df <- bind_cols(all_div)

all_div_means <- data.frame("id" = rownames(all_div[[1]]), divs = rowMeans(all_div_df))
nes_divs <- left_join(all_div_means, nes) %>% 
              mutate(id = str_replace(id, "T[23]", ""))

hist(all_div_means)

# repeatability of diversity
library(rptR)
nes_rpt <- rptGaussian(divs ~ (1|id), grname = "id", data = nes_divs)
nes_rpt

# load microsats ---------------------------------------------------------------
nes_msats <- read_xls("../data/nes_msats_cleaned.xls")
# calc het and g2
nes_geno <- convert_raw(nes_msats[2:ncol(nes_msats)])
g2_microsats(nes_geno)
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



# transform and filter
# ps_trans <- transform_sample_counts(ps, function(x) x / sum(x) )
# ps_filt <- filter_taxa(ps_trans, function(x) mean(x) > 1e-5, TRUE)
# others: 
# subset_taxa(ps, Phylum=="Chlamydiae")
# prune_samples(sample_sums(ps)>=20, ps)

# standardize samples to median sequencing depth
total <- median(sample_sums(ps))
standf <- function(x, t=total) round(t * (x / sum(x)))
stand_ps <- transform_sample_counts(ps, standf)

# Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
stand_ps_filt <- filter_taxa(stand_ps, function(x) sd(x)/mean(x) > 3.0, TRUE)

# Subset the data to Bacteroidetes, used in some plots
#gpsfb <- subset_taxa(stand_ps_filt, Phylum=="Bacteroidetes")
#title <- "plot_bar; Bacteroidetes-only"
#plot_bar(gpsfb, "timepoint", "Abundance", title=title)
#plot_bar(gpsfb, "timepoint", "Abundance", "Family", title=title)

# estimate alpha diversity and clean data.frame
nes_richness <- estimate_richness(ps, split = TRUE, measures = c("Shannon", "Simpson")) %>% 
  mutate(id = rownames(.)) %>% 
  mutate(timepoint = case_when(
    str_detect(id, "T2") ~ "T2",
    str_detect(id, "T3") ~ "T3",
    !((str_detect(id, "T2") | (str_detect(id, "T3")))) ~ "T1"
  )) %>% 
  mutate(id = str_replace(id, "T[23]", "")) %>% 
  filter(!(id == "17BEMa11Fec")) %>% 
  mutate(id = str_replace(id, "17BEMa", "E"))
plot(nes_richness$Simpson)
#nes_rpt <- rptGaussian(Shannon ~ timepoint + (1|id), npermut = 1000, grname = "id", data = nes_richness)
#summary(nes_rpt)

# Show available ranks in the dataset
rank_names(ps)
# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)
# some phyla could not be assigned --> artifacts? --> filter
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))




# Are there phyla that are comprised of mostly low-prevalence features? 
# Compute the total and average prevalences of the features in each phylum.
sum_tab <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# aim: filter anything less abundant than 1%: 
less_than_1perc <- 0.001*(sum(sum_tab[,3]))
# phyla less than 1perc
tofilter_phyla <- as.character(sum_tab[sum_tab[, 3] < less_than_1perc, "Phylum"])
# potential filtering
# Filter entries with unidentified Phylum.
ps1 <- subset_taxa(ps, !Phylum %in% tofilter_phyla)
ps1

# prevalence filtering
# every point is a taxon
# Subset to the remaining phyla
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps1),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Define prevalence threshold as 5% of total samples
prevalenceThreshold <- 0.05 * nsamples(ps1)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 <- prune_taxa(keepTaxa, ps)

# abundance value transformation
plot_abundance = function(physeq, title = "",
Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# Transform to relative abundance. Save as new object.
ps2ra = transform_sample_counts(ps2, function(x){x / sum(x)})
plotBefore = plot_abundance(ps2,"")
plotAfter = plot_abundance(ps2ra,"")

# The histogram suggests that a log(1+x) transformation might be sufficient for  normalizing the abundance data for the exploratory analyses
qplot(log10(rowSums(otu_table(ps2))),binwidth=0.2) +
  xlab("Logged counts-per-sample")

# # do the log transform??
# pslog <- transform_sample_counts(ps, function(x) log(1 + x))
# # ordination
# out_bray_log <- ordinate(pslog, method = "PCoA", distance = "bray")
# evals <- out_bray_log$Eigenvalues
# # plot
# plot_ordination(pslog, out_bray_log, color = "sex") +
#   labs(col = "sex") 


# Linear modeling
nes_richness <- estimate_richness(ps2, split = TRUE, measures = c("Shannon", "Simpson")) %>% 
  mutate(id = rownames(.)) %>% 
  # mutate(timepoint = case_when(
  #   str_detect(id, "T2") ~ "T2",
  #   str_detect(id, "T3") ~ "T3",
  #   !((str_detect(id, "T2") | (str_detect(id, "T3")))) ~ "T1"
  # )) %>% 
  left_join(nes_sampling, by = "id") %>% 
  mutate(id = str_replace(id, "T[23]", "")) %>% 
  filter(!(id == "17BEMa11Fec")) %>% 
  mutate(id = str_replace(id, "17BEMa", "E"))

nes_richness$Shannon_RF 

# load microsats ---------------------------------------------------------------
nes_msats <- read_xls("../data/nes_msats_cleaned.xls")
# calc het and g2
nes_geno <- convert_raw(nes_msats[2:ncol(nes_msats)])
g2_microsats(nes_geno)
# het df
nes_het <- data.frame(id = nes_msats$id, het = sMLH(nes_geno)) %>% 
  arrange(by = het)

# create df for modeling
nes_modeling <- left_join(nes_richness, nes_het, by = "id")

ggplot(nes_modeling, aes(x = sex, y = Shannon)) + geom_boxplot() + geom_jitter()
ggplot(nes_modeling, aes(x = het, y = Shannon)) + geom_point()

# some modeling
alpha_div_model <- lmer(formula = Shannon ~ sex + timepoint + het + (1|id), data = nes_modeling)
summary(alpha_div_model)

ps_dds <- phyloseq_to_deseq2(ps2, design2 = ~ timepoint + sex)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(ps_dds), 1, gm_mean)
ps_dds <- estimateSizeFactors(ps_dds, geoMeans = geoMeans)
diagdds <- DESeq(ps_dds, fitType="local")

ps_dds <- estimateDispersions(ps_dds)
abund <- getVarianceStabilizedData(ps_dds)


ord.nmds.bray <- ordinate(ps, method = "NMDS", distance="bray")
plot_ordination(ps, ord.nmds.bray, color="territory", title="Bray NMDS") +
  theme_minimal() +
  geom_point(size = 3)

# ranks
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))

abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")


# plot richness
plot_richness(ps, x="timepoint", measures=c("Shannon", "Simpson")) + theme_bw()
library(tidyr)
# transform to relative abundances
psra <- transform_sample_counts(ps, function(x){x / sum(x)})
# estimate alpha diversity and clean data.frame
nes_richness <- estimate_richness(pstest, split = TRUE, measures = c("Shannon", "Simpson")) %>% 
                  mutate(id = rownames(.)) %>% 
                  mutate(timepoint = case_when(
                    str_detect(id, "T2") ~ "T2",
                    str_detect(id, "T3") ~ "T3",
                    !((str_detect(id, "T2") | (str_detect(id, "T3")))) ~ "T1"
                  )) %>% 
                  mutate(id = str_replace(id, "T[23]", "")) %>% 
                  filter(!(id == "17BEMa11Fec")) %>% 
                  mutate(id = str_replace(id, "17BEMa", "E"))

plot(nes_richness$Shannon, sample_sums(ps)) 




# plot 
plot_abundance = function(physeq, ylabn = "",
                          Facet = "Order",
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

# transform to relative abundances
psra <- transform_sample_counts(ps, function(x){x / sum(x)})
plot_abundance(ps)
plot_abundance(psra)

# repeatability of diversity
library(rptR)
nes_rpt <- rptGaussian(Shannon ~ timepoint + (1|id), npermut = 1000, grname = "id", data = nes_richness)
summary(nes_rpt)

# load microsats ---------------------------------------------------------------
nes_msats <- read_xls("../data/nes_msats_cleaned.xls")
# calc het and g2
nes_geno <- convert_raw(nes_msats[2:ncol(nes_msats)])
g2_microsats(nes_geno)
# het df
nes_het <- data.frame(id = nes_msats$id, het = sMLH(nes_geno)) %>% 
  arrange(by = het)

# mean bacterial richness
nes_richness_mean <- nes_richness %>% group_by(id) %>% 
                  summarize(mean_shannon = mean(Shannon),
                            mean_simpson = mean(Simpson))
nes_richness_T1 <- nes_richness[nes_richness$timepoint == "T1", ]

nes_full <- left_join(nes_het, nes_richness_T1, by = "id")

ggplot(nes_full, aes(x=het, y=Shannon)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_smooth(method = "lm")

summary(lm(data = nes_full, mean_shannon ~ het))


# load sampling data -----------------------------------------------------------
nes_sampling <- read_xlsx("../data/sampling_data_processed.xlsx") %>% 
                  dplyr::rename(id = ID,
                         date = DATE,
                         territory = TERRITORY, 
                         sex = SEX) %>% 
                  mutate(id = str_replace(id, "17BEMA0", "17BEMa")) %>% 
                  mutate(id = str_replace(id, "0", ""))





# put together microbial and genetic diversity data
nes$id


# some ordination
ord_nmds_bray <- ordinate(ps, method="NMDS", distance="bray")
plot_ordination(ps, ord_nmds_bray, color = "timepoint", title="Bray NMDS")

 


# Construct the phylogenetic tree
seqs <- getSequences(seqtab_nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input")



