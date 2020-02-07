R code for ‘Early sexual dimorphism in the developing gut microbiome of northern elephant seals’
------------------------------------------------------------------------------------------------

Stoffel, M.A., Acevedo-Whitehouse, K., Nami Morales-Durán, Grosser, S.,
Chakarov. N., Krüger, O., Nichols, H.J., Elorriaga-Verplancken, F.R.,
Hoffman, J.I.

The analysis is split into the following four parts:

*1\_reads\_to\_ASVs.R* This is the DADA2 pipeline which goes all the way
from raw reads to amplicon sequence variants (ASVs)

*2\_prep\_data.R* This script does some preprocessing on the ASVs and
integrates other sampling data, such as health.

*3\_filter\_ASVs.R* Quality control and filtering of the ASVs.

*4\_full\_analysis.R* This script contains the code to reproduce all
analyses and figures in the paper and the supplementary material.

### Other scripts

All other scripts are helper scripts which were used for the quality
control of the microsatellites, or contain code for plotting or
manipulation of microbiome data.

### Folders

The data folder contains all the data necessary to run script 4 with all
analyses.

-   health\_data.xlsx contains the health data derived from blood tests
-   nes\_msats\_cleaned.xls contains the microsatellite genotypes
-   ps0.RData is a phyloseq object containing all unfiltered microbiome
    ASVs
-   ps3.RData is the filtered phyloseq object used for most analyses
-   sampling\_data\_processed.xlsx contains all the other sampling data,
    such as date
