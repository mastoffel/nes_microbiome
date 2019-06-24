library(tidyverse)
library(readxl)
library(phyloseq)

# Load and prepare data --------------------------------------------------------

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

# health data inferred from blood cells
health_data_t1 <- read_xlsx("../data/processed/health_data.xlsx") %>% 
    mutate(ID = str_replace(ID, "MA", "Ma"),
           ID = paste0(ID, "T1"))
health_data_t2 <- read_xlsx("../data/processed/health_data.xlsx", sheet = 2) %>% 
    mutate(ID = str_replace(ID, "MA", "Ma"),
           ID = paste0(ID, "T3"))
health_data <- rbind(health_data_t1, health_data_t2) %>% 
            dplyr::rename(health_status = "Health status",
                          id = ID,
                          category = Category) %>% 
            mutate(health_status = ifelse(health_status  == "NA", NA, health_status ))

# read and process other northern elephant seal data
nes_sampling <- read_xlsx("../data/processed/sampling_data_processed.xlsx") %>% 
    dplyr::rename(id = ID,
                  date = DATE,
                  territory = TERRITORY, 
                  sex = SEX) %>% 
    mutate(id = str_replace(id, "17BEMA0", "17BEMa")) %>% 
    mutate(id = str_replace(id, "17BEMA", "17BEMa")) %>% 
    mutate(id = str_c(id, timepoint)) %>% 
    left_join(health_data, by = "id") %>% 
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

# save full phyloseq object
save(ps0, file = "../data/processed/ps0.RData")
