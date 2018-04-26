# create sheet with sample names
library(readxl)
id <- "17BEMa"
t1_samples <- paste0(id, 1:40)
# died or left: 12, 39
t2_samples <- paste0(id, 1:40, "(T2)")[-c(12, 39)]
# died or left: 12, 39, 17, 20, 33, 36
t3_samples <- paste0(id, 1:40, "(T3)")[-c(12, 39, 17, 20, 33, 36)]

feces_sample <- "17BEMa11Fec"

library(WriteXLS)
openxlsx::write.xlsx(data.frame(Box1 = c(t1_samples, t2_samples)), file = "box1samples.xlsx")
openxlsx::write.xlsx(data.frame(Box1 = c(t3_samples, feces_sample)), file = "box2samples.xlsx")

40+38+34 + 1
