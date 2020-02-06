# prep data for ENA
library(tidyverse)
load("data/ps0.RData")

# prep spreadsheet for file submission
file_names <- list.files(path = "../raw_reads/primer_clipped_reads/")
fastq_df <- data.frame(forward_file_name = file_names[str_detect(file_names, "R1")],
                       reverse_file_name = file_names[str_detect(file_names, "R2")]) %>% 
            mutate(sample_name_for = str_match(forward_file_name, "17BE.+_R"),
                   sample_name_rev = str_match(reverse_file_name, "17BE.+_R")) %>% 
            mutate(sample_name = ifelse(sample_name_for == sample_name_rev, str_replace(sample_name_for, "\\_R", ""))) %>% 
            filter(!str_detect(sample_name, "17BEMa11Fec")) %>% 
            mutate(library_name = str_replace(forward_file_name, "\\_R..fastq.bz2", "")) %>% 
            mutate(instrument_model = "Illumina MiSeq",
                   library_source = "METAGENOMIC",
                   library_selection = "PCR",
                   library_strategy = "AMPLICON",
                   design_description = "Amplicon sequencing 300 bp paired-end reads of the 16S region with 341F-785R primers",
                   library_construction_protocol = "Illumina MiSeq Paired End",
                   file_type = "fastq",
                   library_layout = "PAIRED",
                   insert_size = 0
                   ) %>% 
    dplyr::rename(sample_alias = sample_name) %>% 
    select(sample_alias, instrument_model, library_name, library_source,
           library_selection, library_strategy, design_description,
           library_construction_protocol, insert_size, forward_file_name,
           reverse_file_name, file_type, library_layout)

write_tsv(fastq_df, "../ENA_submission/RUN_spreadsheet_ENAsubmission.tsv")



df <- ps0@sam_data %>% 
      as_tibble() %>% 
      dplyr::select(id) %>% 
      mutate(sample_alias = as.character(id)) %>% 
      mutate(instrument_model = "Illumina MiSeq",
             library_name = )





ncbi_dat <- ps0@sam_data %>% 
            as_tibble() %>% 
            dplyr::rename(sample_name = "id") %>% 
            mutate(sex = ifelse(sex == "F", "female", "male")) %>% 
            select(sample_name, sex) %>% 
            mutate(organism = "Mirounga angustirostris",
                   tissue = "rectum",
                   description = "Rectal swab using FLOQSwab")

WriteXLS::WriteXLS(ncbi_dat, "../ncbi_submission/BioSample.xls")                

