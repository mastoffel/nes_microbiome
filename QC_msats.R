# Quality control for microsatellite genotypes
# check microsat genotypes
library(inbreedR)
library(readxl)
library(stringr)
library(dplyr)
library(Demerelate)

# Part 1: Cleaning raw data ----------------------------------------------------
nes <- read_xlsx("../data/raw/elesealgenotypes_12_03.xlsx") 

# delete loci with %
nes <- nes[1:47]

# data cleaning
# naming
loci_names <- names(nes[seq(from = 2, to = ncol(nes), by = 2)]) %>% 
    str_replace("-", "") %>% 
    str_replace("\\.", "_") %>% 
    str_replace("\\*", "s") %>% 
    str_replace("__", "_") %>% 
    str_replace("%", "")

loci_names_a <- paste0(loci_names, "_a")
loci_names_b <- paste0(loci_names, "_b")

names(nes)[seq(from = 2, to = ncol(nes), by = 2)] <- loci_names_a
names(nes)[seq(from = 3, to = ncol(nes), by = 2)] <- loci_names_b

names(nes)
str(nes)

library(dplyr)
#change ** to NA

nes <- nes %>% 
    rename("id"= "X__1") #%>% 
  #  mutate(id = str_replace(id, "p1.fsa", "")) 
   # mutate_if(is.character, funs(replace(., .["**"], NA),. ))

replace_sym <- function(x){
    if (is.character(x)) {
        x <- na_if(x, "**")
    }
as.numeric(x)
}
nes[2:ncol(nes)] <- lapply(nes[2:ncol(nes)], replace_sym)
str(nes)
lapply(nes, table, useNA = "ifany")

# check if polymorphic
which(colSums(nes[seq(from = 2, to = ncol(nes), by = 2)] - nes[seq(from = 3, to = ncol(nes), by = 2)], na.rm = TRUE) == 0)

# filter out Mang09s locus
nes <- nes %>% select_if(!grepl("Mang09s", names(.)))
nes <- nes %>% select_if(!grepl("unclear", names(.)))

# all to integer
nes[, 2:ncol(nes)] <- lapply(nes[, 2:ncol(nes)], as.integer)
nes
# WriteXLS::WriteXLS(nes, "../data/raw/nes_msats_cleaned.xls")



# Part 2: HW checking
# check HW etc
library(pacman)
p_load(Demerelate, dplyr, pegas, ape, seqinr)

nes_msats <- as.data.frame(readxl::read_xls("../data/raw/nes_msats_cleaned.xls"))

gene_start <- 2
# create file with genotypes stored as "134/140" from two-column per locus format
change_geno_format <- function(seal_data_df) {
    loci_names <- names(seal_data_df[gene_start:ncol(seal_data_df)])
    # filter loci names and remove the allel add-on
    loci_names <- str_replace_all(loci_names[seq(from = 1, to = length(loci_names), by = 2)], "_a", "")
    short_geno <- data.frame(matrix(nrow = nrow(seal_data_df), ncol = length(loci_names)))
    names(short_geno) <- loci_names
    
    genotypes <- seal_data_df[, gene_start:ncol(seal_data_df)]
    length_data <- ncol(genotypes)
    
    col_num <- 1
    for (i in seq(from = 1, to = length_data, by = 2)) {
        short_geno[col_num] <- paste(genotypes[[i]], genotypes[[i+1]], sep = "/")
        col_num <- col_num + 1
    }
    
    short_geno[short_geno == "NA/NA"] <- NA
    geno_out <- cbind(seal_data_df[c(1)], short_geno)
    geno_out
}

nes_msats_genind <- change_geno_format(nes_msats)

check_na_genotypes <- function(df) {
    df[gene_start:ncol(df)] <- lapply(df[gene_start:ncol(df)], function(x) {
        if (sum(grepl("NA", x))>0) {
            x[which(grepl("NA", x))] <- NA
        }
        x
    })
    df
}

nes_msats_locus_format <- check_na_genotypes(nes_msats_genind)

#save and reload as genind object
if (dir.exists("output/genind_formatted")) {
    system(paste("rm -r output/genind_formatted"))
}
system("mkdir output/genind_formatted")

write.table(nes_msats_locus_format, file = paste("output/genind_formatted/nes_genind", ".txt", sep = ""
), sep = " ", quote = FALSE, row.names = FALSE)

nes_msats_pegas <- read.loci(file = paste("output/genind_formatted/nes_genind", ".txt", sep = ""
),  header = TRUE, loci.sep = " ", allele.sep = "/", col.loci = 2:ncol(nes_msats_locus_format))

# HW test
nes_hw <- hw.test(nes_msats_pegas, B = 1000)
nes_hw <- as.data.frame(nes_hw)
p.adjust(nes_hw$Pr.exact, method = "fdr") # none out of HW

