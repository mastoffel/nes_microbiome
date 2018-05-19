# check msat genotypes
library(inbreedR)
library(readxl)
library(stringr)
library(dplyr)
library(Demerelate)
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


library(readxl)
library(tibble)
library(dplyr)
library(Demerelate)
nes <- readxl::read_xls("../data/raw/nes_msats_cleaned.xls") %>% 
        add_column(factor(rep("SB", nrow(.))), .after = 1)
names(nes)[1:2] <- c("Sample-ID", "Population")

nes_rel <- Demerelate(as.data.frame(nes), value= "wang", file.output=FALSE, object = TRUE, pairs = 100)
hist(unlist(nes_rel$Empirical_Relatedness))

Loci.test(as.data.frame(nes), bt=1000, ref.pop="NA", object=TRUE, value= "wang", file.output=TRUE)
#Emp.calc(as.data.frame(nes), value="rxy",ref.pop=NA)

# prepare file for related package

library(related)
nes_rel <- coancestry(nes, wang =1)
nes_rel$relatedness

# simulations
sim <- familysim(nes_rel$freqs , 100)
output <- coancestry( sim , wang =1)
simrel <- cleanuprvals(output$relatedness , 100)
# plotting
relvalues <- simrel[, 6]
label1 <- rep("PO", 100)
label2 <- rep("Full", 100)
label3 <- rep("Half", 100)
label4 <- rep("Unrelated", 100)

labels <- c( label1 , label2 , label3 , label4 )
plot ( as.factor(labels) , relvalues , ylab =" Relatedness Value ", xlab =" Relatedness ")

Relationship <- labels
newdata <- as.data.frame(cbind(Relationship , relvalues))
newdata$relvalues <- as.numeric(as.character(newdata$relvalues))

ggplot(newdata, aes(x= relvalues, by = Relationship, colour = as.factor(Relationship))) + geom_density()

# which estimator is best?
compareestimators(nes_rel , 100)

test_rel <- coancestry ( GenotypeData , wang =2)

library(inbreedR)
nes_geno <- convert_raw(nes[2:ncol(nes)])
g2_microsats(nes_geno)
hist(sMLH(nes_geno))

plot(sMLH(nes_geno))

nes_het <- data.frame(id = nes$id, het = sMLH(nes_geno)) %>% 
                arrange(by = het)



