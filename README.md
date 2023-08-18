# QSeq

This package is intended to be used as part of R based microbiome analysis pipelines to transform compositional microbiome data into Quantitative Sequence data as described in (1, 2). For a detailed description of the context, see (3).

Some jargon defined:

OTU = Operational Taxonomic Unit. Used interchangeably with species concepts for microbiome research, typically created by clustering several sequence variants into one infered "operational taxonomic unit."

ASV = Amplicon Sequence Variant. A type of OTU where error inference is used to identify unique, biologically relevant variants of the amplicon that was sequenced.

QPCR = Quantitative Polymerase Chain Reaction. This is a lab method for counting the number of a particular gene in a sample.

This software uses phyloseq to organize the data, and requires a phyloseq object as input. Instructions on how to create phyloseq objects from OTU tables and sample data can be found [here](https://joey711.github.io/phyloseq/install.html#problem_with_r_core_version_number). Instructions on how to process marker gene sequence data into an OTU table in R using DADA2, with a hand-off to phyloseq can be found [here](https://benjjneb.github.io/dada2/index.html). 

This package requires two pieces of input data. The first is an OTU table. The second is a quantitative measure of total abundance. Both of these should be included in a phyloseq object. One option is to use QPCR to determine the total number of marker genes in a sample used for sequencing, but other methods could be used. Conceptually, this software uses sequencing to determine the composition of the sample, and an independent measure of abundance to determine total abundance for each sample. It operates by scaling the composition of each sample by the total abundance of that sample.

QSeq houses one primary function that operates as follows:

1) OTU count data is extracted from a phyloseq object and transformed into relative abundance.
2) The relative abundance values for each sample are multiplied by that sample's total abundance.
3) Sequence counts are rounded to the nearest whole integer to avoid partial sequence counts.
4) The new OTU table is merged back into phyloseq, and a phyloseq object is returned with the original metadata and taxonomy, but updated abundance values for each OTU in each sample.

QSeq is most useful for differential abundance testing, and correlation inference. It offers some benefit for alpha diversity testing; and depending on the question might enhance or detract from beta diversity analysis.

# Installation and Tutorial

To install the package from GitHub, first install devtools:
```{r}
install.packages("devtools")
```

Then to install QSeq from GitHub:
```{r}
library(devtools)
install_github("djeppschmidt/QSeq")
```

QSeq has one primary function that takes the following arguments:

ps = a phyloseq object with an OTU table, and sample data that includes a quantitative measure of total abundance for each sample (required).

abundance = name of column in the metadata that contains the abundance data (required).

To run:

```{r}
# ps is your phyloseq object
# replace "col.name" with the column name from your metadata

# not run:

library(phyloseq)
library(QSeq)

abundance<-"col.name"
Q.ps<-QSeq(ps, abundance)

# alternatively, can be run as:
Q.ps<-QSeq(ps, abundance="col.name")
```
You may now conduct the rest of your community analysis.

Please cite this work with: 

Epp Schmidt, DS, et al. (2023). QUANTITATIVE AMPLICON SEQUENCING IS NECESSARY TO IDENTIFY DIFFERENTIAL TAXA AND CORRELATED TAXA WHERE POPULATIONS SIZES DIFFER. Microbial Ecology. DOI: 10.1007/s00248-023-02273-z![image](https://github.com/Djeppschmidt/QSeq/assets/19291020/684950a1-346f-4ab7-b8e2-208bdd8481e8)

# NOTES

A note on rounding:
It is possible that rounding might generate zeroes in low abundance taxa. If this happens, I suggest revisiting your quantification method because it indicates a mismatch between the rank abundance curve and total abundance (more detected features than are measured in total abundance in sample). This effect indicates that your composition was either sampled from a larger population than what was sampled for quantification; or the abundance values for each sample were either transformed (for example, they are sometimes transformed to log scale), or not calculated appropriately (for example, ensure that QPCR has been calculated on an equal mass or volume bases to compare to the input for sequencing). In my work, this means ensuring that the QPCR abundance is calculated to represent the total abundance of bacteria that I expect from the mass of soil that was used for the extraction. In principle, calculating based on the total DNA input to the reaction should work as well. And, it should also be noted that rounding does add a small amount of variation to the dataset.

# Citations:

1) Jian et al., 2020. DOI: [10.1371/journal.pone.0227285](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0227285)
2) Epp Schmidt et al., 2022. DOI: [10.1016/j.apsoil.2022.104396](https://doi.org/10.1016/j.apsoil.2022.104396)
3) Epp Schmidt et al., (in prep) DOI: [10.1007/s00248-023-02273-z](https://link.springer.com/article/10.1007/s00248-023-02273-z)
