# prepare microsats for related

library(readxl)
library(readr)
nes_msats <- read_xls("../data/nes_msats_cleaned.xls")
write_delim(nes_msats, path = "../data/nes_msats_related.txt", delim = " ", na = "0", col_names = FALSE)

# test related
msats <- readgenotypedata("../data/nes_msats_related.txt")

## Simulate for 100 individuals to assess power of the analysis
sim <- familysim(msats$freqs, 100)
relsim <- coancestry(sim , wang = 1)

relsim <- cleanuprvals(relsim$relatedness , 100)
## Extract only the column containing the wang estimates
relvalues <- as.numeric(relsim[,"wang"])  
label1  <- rep("PO", 100)
label2  <- rep("Full", 100)
label3  <- rep("Half", 100)
label4  <- rep("Unrelated", 100)
labels <- c(label1 , label2 , label3 , label4)
relsimtab <- as.data.frame(cbind(relvalues,labels),stringsAsFactors=FALSE)
relsimtab$relvalues <- as.numeric(relsimtab$relvalues)

## Calculate relatedness (wang estimator) for the fur seal individuals
rel <- coancestry(msats$gdata, wang = 1)

relvals <- rel$relatedness[,c("pair.no","ind1.id","ind2.id","wang")]

ggplot(relsimtab, aes(x=labels, y=relvalues)) + 
  geom_boxplot(fill="lightgrey") + 
  theme(legend.position='none') +  
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +   
  ylab("Relatedness Estimate (Wang)")+
  xlab("Relatedness Category")+
  coord_cartesian(ylim = c(-0.3, 0.7))+
  theme(axis.text.x = element_text(size=10), axis.title.x = element_text(margin = margin(10, 0, 0, 0)))+
  theme(axis.title.y=element_text(margin = margin(0, 15, 0, 0)), axis.text.y = element_text(size=10))

## Make a boxplot for the pairs and unrelated categories and mark the outlier points with the number of the comparison.
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
relvals2[,"pairNames"] <-  paste(relvals2[,2], relvals2[,3], sep="")

b <-relvals %>%
  group_by(pairs) %>%
  mutate(outlier = ifelse(is_outlier(wang), pairNames , as.numeric(NA))) %>%
  ggplot(., aes(x = factor(pairs), y = wang)) +
  geom_boxplot(fill="lightgrey") +
  xlab("Relatedness Category")+
  ylab("")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +  
  coord_cartesian(ylim = c(-0.3, 0.7))+
  theme(axis.text.x = element_text(size=10), axis.title.x = element_text(margin = margin(10, 0, 0, 0)))+
  theme(axis.title.y=element_text(margin = margin(0, 15, 0, 0)), axis.text.y = element_text(size=10))+
  geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3,size=2)

grid.arrange(q,b, ncol=2)
