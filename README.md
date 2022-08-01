# QSeq

This package is intended to be used as part of R based microbiome analysis pipelines to transform compositional microbiome data into Quantitative Sequence data as described in (1, 2). For a detailed description of the context, see (3).

Some jargon defined:

OTU = Operational Taxonomic Unit. Used interchangeably with species concepts for microbiome research, typically created by clustering several sequence variants into one infered "operational taxonomic unit."
ASV = Amplicon Sequence Variant. A type of OTU where error inference is used to identify unique, biologically relevant variants of the amplicon that was sequenced.
QPCR = Quantitative Polymerase Chain Reaction. This is a lab method for counting the number of a particular gene in a sample.

This software uses phyloseq to organize the data, and requires a phyloseq object as input. Instructions on how to create phyloseq objects from OTU tables and sample data can be found [here](https://joey711.github.io/phyloseq/install.html#problem_with_r_core_version_number). Instructions on how to process marker gene sequence data into an OTU table in R using DADA2, with a hand-off to phyloseq can be found [here](https://benjjneb.github.io/dada2/index.html). 

This package requires two pieces of input data. The first is an OTU table. The second is a quantitative measure of total abundance. One option is to use QPCR to determine the total number of marker genes in a sample used for sequencing, but other methods could be used. Conceptually, this software uses sequencing to determine the composition of the sample, and an independent measure of abundance to determine total abundance for each sample. It operates by scaling the composition of each sample by the total abundance of that sample.

houses one primary function that operates as follows:

1) OTU count data is transformed into relative abundance.
2) The relative abundance values for each sample are multiplied by that sample's total abundance.
3) Sequence counts are rounded to the nearest whole integer to avoid partial sequence counts.
4) The new OTU table is merged back into phyloseq, and a phyloseq object is returned with the original metadata and taxonomy, but updated abundance values for each OTU in each sample.

In most cases, this QSeq operation should not change the composition of the data, and so it can be done before alpha diversity calculations. It is most useful for beta diversity and differential abundance testing, so I prefer to do alpha diversity calculations on unmodified data before QSeq transformation. 

It is possible that rounding might generate zeroes. If this happens, I suggest revisiting your quantification method. This effect indicates that your composition was either sampled from a larger population than what was sampled for quantification; or the abundance values for each sample were either transformed (for example, they are sometimes transformed to log scale), or not calculated appropriately (for example, ensure that QPCR has been calculated on an equal mass or volume bases to compare to the input for sequencing). In my work, this means ensuring that the QPCR abundance is calculated to represent the total abundance of bacteria that I expect from the mass of soil that was used for the extraction. In principle, calculating based on the total DNA input to the reaction should work as well.

A note on rarefaction:
Rarefaction is sometimes used to control the sampling effort for each sample. This is most appropriate for standardizing the data before calculating alpha diversity metrics. QSeq should not materially affect alpha diversity calculations, and rather is most useful for beta diversity and differential abundance analysis. From this perspective, I do not include any subsampling routine in this package. If you want to model to account for sampling error, I recommend capturing the total number of sequences per sample prior to implementing QSeq, since QSeq will modify the OTU table and obscure this information.

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

abundance<-"col.name"
Q.ps<-QSeq(ps, abundance)

# alternatively, can be run as:
Q.ps<-QSeq(ps, abundance="col.name")
```
You may now conduct the rest of your community analysis.

Please cite this package with: 

[insert citation when published]

# Citations:

1) Jian et al., 2020. DOI: []()
2) Epp Schmidt et al., 2022. DOI: []()
3) Epp Schmidt et al., (in prep) DOI: []()