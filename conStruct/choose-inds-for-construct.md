Find High Quality Gopher Tortoise Individuals for Construct Analysis
================
Alex Krohn[^1]
5/16/2022

## Introduction

This file chooses high quality Gopher Tortoise individuals to use for
construct analysis for the 9 subsets determined through hierarchical
structure modelling. The subsets are all sampling locations,
Stimpson+Perdido, Stimpson, Perdido, Conecuh+Geneva+Rayonier+SD,
Conecuh, SD, Geneva+Rayonier, Conecuh+SD.

The overall plan for creating datasets for conStruct analysis will be as
follows: 1) Find individuals from the metadata that meet the subsetting
requirements. 2) Filter sites by creating a new VCF with filtered sites
present in the subsetted individuals. 3) Take one random SNP per RAD
locus using this script from [Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py).
4) Convert the VCF to Structure format using PGD Spider 5) Import the
STR file, and filter individuals to make sure no individual has more
than 10% missing data. 6) For the remaining individuals, ensure that no
individual’s PC1 or PC2 scores correlate with their amount of missing
data. 7) Use that remaining dataset for conStruct analyses.

## All Sampling Locations —-

Note that conStruct requires many more loci than individuals, so we
cannot include all individuals. So, here I’ll subset to the 5
individuals with the most loci from each sampling location
(`site.fixed.simplified`).

``` r
library(tidyverse)

system("mkdir allinds")
setwd("allinds")


# Import metadata
metadata <- read_csv("../data/modified-metadata-for-seqd-individuals.csv") 

high.quality.inds.to.keep <- metadata %>%
  filter(!is.na(lat)) %>% # Make sure they also have lat/long data
  group_by(Site.fixed.simplified) %>%
  top_n(5, loci) %>%
  dplyr::select(vcfname)

write.table(high.quality.inds.to.keep[,2], "all-sites-inds-to-keep.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

From the vcf output by ipyrad, filter sites heavily to only include
biallelic sites present in 90% of individuals. Then use this script from
[Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py)
to randomly sample one SNP per RAD locus to minimize the effect of LD.

``` bash
vcftools --vcf ../data/merged-datasets.vcf --keep all-sites-inds-to-keep.txt --max-missing 0.9 --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out all-sites-maxmissing9-biallelic
```

``` bash
python single_snp.py all-sites-maxmissing9-biallelic.recode.vcf 
```

That leaves 7405 SNPs in the all sites dataset.

Then use [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) (Java
applet) to convert that final VCF to Structure format.

Next, import the STR file, check individual-level missingness, re-filter
as needed.

``` r
library(conStruct)

# Import the allele data for best inds
construct.data <- structure2conStruct(infile = "mod_all-sites-maxmissing9-biallelic.recode.vcf.str",
    onerowperind = FALSE, start.loci = 3, start.samples = 2, missing.datum = -9,
    outfile = "all-sites-dataset")

# OR load the existing data load('all-sites-dataset.RData') construct.data <-
# freqs rm(freqs)

# Import lat/long data.
metadata <- read_csv("../data/modified-metadata-for-seqd-individuals.csv")

# Recheck missingness of individuals
missingness <- apply(construct.data, 1, function(x) sum(is.na(x))/ncol(construct.data))
hist(missingness)

# How many inds have less than 20% missing data?
length(names(missingness[missingness < 0.2]))  #56

# Do they come from all populations?
metadata %>%
    filter(vcfname %in% names(missingness[missingness < 0.2])) %>%
    count(Site.fixed.simplified)  #yes, but it leaves only 1 individual from site 5

# Save the final list of individuals
final.inds.list <- metadata %>%
    filter(vcfname %in% names(missingness[missingness < 0.1])) %>%
    select(vcfname)
```

Finally, run a PCA on the final conStruct dataset to make sure neither
PC1 nor PC2 correlate with the level of missing data.

``` r
library(adegenet)

# Load the data using adegenet
allinds.genepop <- read.structure(file = "mod_all-sites-maxmissing9-biallelic.recode.vcf.str",
    n.ind = 60, n.loc = 7405, col.lab = 1, col.pop = 2, row.marknames = 1, col.others = NA,
    onerowperind = FALSE, ask = TRUE)

# Subset to just the construct individuals
final.construct.dataset <- allinds.genepop[indNames(allinds.genepop) %in% final.inds.list$vcfname]

# Scale NAs by mean value, then run PCA
all.sites.scaled <- tab(final.construct.dataset, freq = TRUE, NA.method = "mean")
first.pca <- dudi.pca(all.sites.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
    nf = 10)

# Calculate missingness again
missingness <- apply(final.construct.dataset$tab, 1, function(x) sum(is.na(x))/ncol(final.construct.dataset$tab))

# Calculate correlation between PC1 and PC2, by regression
summary(lm(first.pca$li[, 1] ~ missingness))  # P = 0.44
summary(lm(first.pca$li[, 2] ~ missingness))  # P = 0.2
```

## Stimpson and Perdido—–

``` r
library(tidyverse)

system("mkdir ../stimpson-perdido")
setwd("../stimpson-perdido")

# Filter metadata, keep the top 10 best individuals
high.quality.inds.to.keep <- metadata %>%
  filter(Site.fixed.simplified == "Stimpson" | Site.fixed.simplified == "Perdido") %>%
  filter(!is.na(lat)) %>% # Make sure they also have lat/long data
  group_by(Site.fixed.simplified) %>% 
  top_n(10, loci) %>%
  dplyr::select(vcfname)

write.table(high.quality.inds.to.keep[,2], "stimp-perdido-inds-to-keep.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

From the vcf output by ipyrad, filter sites heavily to only include
biallelic sites present in 90% of individuals. Then use this script from
[Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py)
to randomly sample one SNP per RAD locus to minimize the effect of LD.

``` bash
vcftools --vcf ../../../merged-datasets.vcf --keep stimp-perdido-inds-to-keep.txt --max-missing 0.9 --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out stimp-perdido-maxmissing9-biallelic
```

``` bash
python ~/Documents/NorthCarolina/TangledBank/tbc-hub/3rad-analyses/single_snp.py stimp-perdido-maxmissing9-biallelic.recode.vcf 
```

That leaves 6921 SNPs in the Stimpson/Perdido dataset.

Then use [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) (Java
applet) to convert that final VCF to Structure format.

Next, import the STR file, check individual-level missingness, re-filter
as needed.

``` r
# Import the allele data for best inds
genepop <- read.structure(file = "mod_stimp-perdido-maxmissing9-biallelic.recode.vcf.str",
    n.ind = 20, n.loc = 6921, col.lab = 1, col.pop = 2, row.marknames = 1, col.others = NA,
    onerowperind = FALSE, ask = TRUE)

# Recheck missingness of individuals
missingness <- apply(genepop$tab, 1, function(x) sum(is.na(x))/ncol(genepop$tab))
hist(missingness)

# How many inds have less than 20% missing data?
length(names(missingness[missingness < 0.2]))  # all 20!
```

Finally, run a PCA on the final conStruct dataset to make sure neither
PC1 nor PC2 correlate with the level of missing data.

``` r
# Scale NAs by mean value, then run PCA
genepop.scaled <- tab(genepop, freq = TRUE, NA.method = "mean")
first.pca <- dudi.pca(genepop.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
    nf = 10)

# Calculate correlation between PC1 and PC2, by regression
summary(lm(first.pca$li[, 1] ~ missingness))  # P = 0.23
summary(lm(first.pca$li[, 2] ~ missingness))  # P = 0.423
```

## Conecuh + Solon Dixon + Geneva + Rayonier —–

``` r
library(tidyverse)

system("mkdir ../conecuh-sd-geneva-rayonier")
setwd("../conecuh-sd-geneva-rayonier/")

# Filter metadata, keep the top 10 best individuals
high.quality.inds.to.keep <- metadata %>%
  filter(Site.fixed.simplified != "Stimpson" & Site.fixed.simplified != "Perdido") %>%
  filter(!is.na(lat)) %>% # Make sure they also have lat/long data
  group_by(Site.fixed.simplified) %>% 
  top_n(10, loci) %>%
  dplyr::select(vcfname)

write.table(high.quality.inds.to.keep[,2], "conecuh-geneva-rayonier-sd-inds-to-keep.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

From the vcf output by ipyrad, filter sites heavily to only include
biallelic sites present in 90% of individuals. Then use this script from
[Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py)
to randomly sample one SNP per RAD locus to minimize the effect of LD.

``` bash
vcftools --vcf ../../merged-datasets.vcf --keep conecuh-geneva-rayonier-sd-inds-to-keep.txt --max-missing 0.9 --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out conecuh-geneva-rayonier-sd-maxmissing9-biallelic
```

``` bash
python ~/Documents/NorthCarolina/TangledBank/tbc-hub/3rad-analyses/single_snp.py conecuh-geneva-rayonier-sd-maxmissing9-biallelic.recode.vcf 
```

That leaves 3691 SNPs in the Conecuh+Geneva+Rayonier+SolonDixon dataset.

Then use [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) (Java
applet) to convert that final VCF to Structure format.

Next, import the STR file, check individual-level missingness, re-filter
as needed.

``` r
# Import the allele data for best inds
genepop <- read.structure(file = "mod_conecuh-geneva-rayonier-sd-maxmissing9-biallelic.recode.vcf.str",
    n.ind = 100, n.loc = 3691, col.lab = 1, col.pop = 2, row.marknames = 1, col.others = NA,
    onerowperind = FALSE, ask = TRUE)

# Recheck missingness of individuals
missingness <- apply(genepop$tab, 1, function(x) sum(is.na(x))/ncol(genepop$tab))
hist(missingness)

# How many inds have less than 20% missing data?
length(names(missingness[missingness < 0.2]))  #90

# Do they come from all populations?
metadata %>%
    filter(vcfname %in% names(missingness[missingness < 0.2])) %>%
    count(Site.fixed.simplified)  #yes, but it leaves only 2 individuals from site 5

# Save the final list of individuals
final.inds.list <- metadata %>%
    filter(vcfname %in% names(missingness[missingness < 0.2])) %>%
    select(vcfname)
```

Finally, run a PCA on the final conStruct dataset to make sure neither
PC1 nor PC2 correlate with the level of missing data.

``` r
# Filter original datset to just the final inds
genepop.filtered <- genepop[indNames(genepop) %in% final.inds.list$vcfname]

# Recalculate missingness
missingness <- apply(genepop.filtered$tab, 1, function(x) sum(is.na(x))/ncol(genepop.filtered$tab))

# Scale NAs by mean value, then run PCA
genepop.scaled <- tab(genepop.filtered, freq = TRUE, NA.method = "mean")
first.pca <- dudi.pca(genepop.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
    nf = 10)

# Calculate correlation between PC1 and PC2, by regression
summary(lm(first.pca$li[, 1] ~ missingness))  # P = 0.99
summary(lm(first.pca$li[, 2] ~ missingness))  # P = 0.92
```

## Conecuh Only —–

``` r
library(tidyverse)

setwd("~/Documents/NorthCarolina/TangledBank/gopher-tortoise/data-analysis/construct/new-way-of-subsetting/conecuh-only/")

# Filter metadata, keep the top 20 best individuals
high.quality.inds.to.keep <- metadata %>%
  filter(Site.fixed.simplified != "Stimpson" & Site.fixed.simplified != "Perdido" & Site.fixed.simplified != "Geneva" & Site.fixed.simplified != "SolonDixon" & Site.fixed.simplified != "Rayonier") %>%
  filter(!is.na(lat)) %>% # Make sure they also have lat/long data
  group_by(Site.fixed.simplified) %>% 
  top_n(10, loci) %>%
  dplyr::select(vcfname)

write.table(high.quality.inds.to.keep[,2], "conecuh-only-inds-to-keep.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

From the vcf output by ipyrad, filter sites heavily to only include
biallelic sites present in 90% of individuals. Then use this script from
[Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py)
to randomly sample one SNP per RAD locus to minimize the effect of LD.

``` bash
vcftools --vcf ../../data/merged-datasets.vcf --keep conecuh-only-inds-to-keep.txt --max-missing 0.9 --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out conecuh-only-maxmissing9-biallelic
```

``` bash
python single_snp.py conecuh-only-maxmissing9-biallelic.recode.vcf 
```

That leaves 2694 SNPs in the Conecuh-only dataset.

Then use [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) (Java
applet) to convert that final VCF to Structure format.

Next, import the STR file, check individual-level missingness, re-filter
as needed.

``` r
# Import the allele data for best inds
genepop <- read.structure(file = "mod_conecuh-only-maxmissing9-biallelic.recode.vcf.str",
    n.ind = 70, n.loc = 2694, col.lab = 1, col.pop = 2, row.marknames = 1, col.others = NA,
    onerowperind = FALSE, ask = TRUE)

# Recheck missingness of individuals
missingness <- apply(genepop$tab, 1, function(x) sum(is.na(x))/ncol(genepop$tab))

# How many inds have less than 20% missing data?
length(names(missingness[missingness < 0.2]))  #62

# Do they come from all populations?
metadata %>%
    filter(vcfname %in% names(missingness[missingness < 0.2])) %>%
    count(Site.fixed.simplified)  #yes, but it leaves only 4 individuals from site 5

# Save the final list of individuals
final.inds.list <- metadata %>%
    filter(vcfname %in% names(missingness[missingness < 0.2])) %>%
    select(vcfname)
```

Finally, run a PCA on the final conStruct dataset to make sure neither
PC1 nor PC2 correlate with the level of missing data.

``` r
# Filter original datset to just the final inds
genepop.filtered <- genepop[indNames(genepop) %in% final.inds.list$vcfname]

# Recalculate missingness
missingness <- apply(genepop.filtered$tab, 1, function(x) sum(is.na(x))/ncol(genepop.filtered$tab))

# Scale NAs by mean value, then run PCA
genepop.scaled <- tab(genepop.filtered, freq = TRUE, NA.method = "mean")
first.pca <- dudi.pca(genepop.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
    nf = 10)

# Calculate correlation between PC1 and PC2, by regression
summary(lm(first.pca$li[, 1] ~ missingness))  # P = 0.60
summary(lm(first.pca$li[, 2] ~ missingness))  # P = 0.59
```

## Geneva + Rayonier

``` r
library(tidyverse)

system("mkdir ../geneva-rayonier")
setwd("../geneva-rayonier/")

# Filter metadata
high.quality.inds.to.keep <- metadata %>%
  filter(Site.fixed.simplified == "Geneva" | Site.fixed.simplified == "Rayonier") %>%
  filter(!is.na(lat)) %>% # Make sure they also have lat/long data
  group_by(Site.fixed.simplified) %>% 
  top_n(10, loci) %>%
  dplyr::select(vcfname)

write.table(high.quality.inds.to.keep[,2], "geneva-rayonier-inds-to-keep.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

From the vcf output by ipyrad, filter sites heavily to only include
biallelic sites present in 90% of individuals. Then use this script from
[Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py)
to randomly sample one SNP per RAD locus to minimize the effect of LD.

``` bash
vcftools --vcf ../../data/merged-datasets.vcf --keep geneva-rayonier-inds-to-keep.txt --max-missing 0.9 --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out geneva-rayonier-maxmissing9-biallelic
```

``` bash
python single_snp.py geneva-rayonier-maxmissing9-biallelic.recode.vcf 
```

That leaves 11096 SNPs in the Geneva-Rayonier dataset.

Then use [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) (Java
applet) to convert that final VCF to Structure format.

Next, import the STR file, check individual-level missingness, re-filter
as needed.

``` r
# Import the allele data for best inds
genepop <- read.structure(file = "mod_geneva-rayonier-maxmissing9-biallelic.recode.vcf.str",
    n.ind = 20, n.loc = 11096, col.lab = 1, col.pop = 2, row.marknames = 1, col.others = NA,
    onerowperind = FALSE, ask = TRUE)

# Recheck missingness of individuals
missingness <- apply(genepop$tab, 1, function(x) sum(is.na(x))/ncol(genepop$tab))

# How many inds have less than 20% missing data?
length(names(missingness[missingness < 0.2]))  # all 20!
```

Finally, run a PCA on the final conStruct dataset to make sure neither
PC1 nor PC2 correlate with the level of missing data.

``` r
# Scale NAs by mean value, then run PCA
genepop.scaled <- tab(genepop, freq = TRUE, NA.method = "mean")
first.pca <- dudi.pca(genepop.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
    nf = 10)

# Calculate correlation between PC1 and PC2, by regression
summary(lm(first.pca$li[, 1] ~ missingness))  # P = 0.053
summary(lm(first.pca$li[, 2] ~ missingness))  # P = 0.70
```

## Perdido only

``` r
library(tidyverse)

system("mkdir ../perdido only")
setwd("../perdido-only/")

# Filter metadata, keep the top 10 best individuals
high.quality.inds.to.keep <- metadata %>%
  filter(Site.fixed.simplified == "Perdido") %>%
  filter(!is.na(lat)) %>% # Make sure they also have lat/long data
  group_by(Site.fixed.simplified) %>%
  top_n(10, loci) %>%
  dplyr::select(vcfname)

write.table(high.quality.inds.to.keep[,2], "perdido-inds-to-keep.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

From the vcf output by ipyrad, filter sites heavily to only include
biallelic sites present in 90% of individuals. Then use this script from
[Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py)
to randomly sample one SNP per RAD locus to minimize the effect of LD.

``` bash
vcftools --vcf ../../data/merged-datasets.vcf --keep perdido-inds-to-keep.txt --max-missing 0.9 --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out perdido-only-maxmissing9-biallelic
```

``` bash
python ingle_snp.py perdido-only-maxmissing9-biallelic.recode.vcf 
```

That leaves 12816 SNPs in the Perdido-only dataset.

Then use [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) (Java
applet) to convert that final VCF to Structure format.

Next, import the STR file, check individual-level missingness, re-filter
as needed.

``` r
# Import the allele data for best inds
genepop <- read.structure(file = "mod_perdido-only-maxmissing9-biallelic.recode.vcf.str",
    n.ind = 10, n.loc = 12816, col.lab = 1, col.pop = 2, row.marknames = 1, col.others = NA,
    onerowperind = FALSE, ask = TRUE)

# Recheck missingness of individuals
missingness <- apply(genepop$tab, 1, function(x) sum(is.na(x))/ncol(genepop$tab))

# How many inds have less than 20% missing data?
length(names(missingness[missingness < 0.2]))  # all!
```

Finally, run a PCA on the final conStruct dataset to make sure neither
PC1 nor PC2 correlate with the level of missing data.

``` r
# Scale NAs by mean value, then run PCA
genepop.scaled <- tab(genepop, freq = TRUE, NA.method = "mean")
first.pca <- dudi.pca(genepop.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
    nf = 10)

# Calculate correlation between PC1 and PC2, by regression
summary(lm(first.pca$li[, 1] ~ missingness))  # P = 0.89
summary(lm(first.pca$li[, 2] ~ missingness))  # P = 0.46
```

## Stimpson only—-

``` r
library(tidyverse)

system("mdkir ../stimp-only")
setwd("../stimp-only/")

# Filter metadata, keep the top 10 best individuals
high.quality.inds.to.keep <- metadata %>%
  filter(Site.fixed.simplified == "Stimpson") %>%
  filter(!is.na(lat)) %>% # Make sure they also have lat/long data
  group_by(Site.fixed.simplified) %>%
  top_n(10, loci) %>%
  dplyr::select(vcfname)

write.table(high.quality.inds.to.keep[,2], "stimpson-inds-to-keep.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

From the vcf output by ipyrad, filter sites heavily to only include
biallelic sites present in 90% of individuals. Then use this script from
[Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py)
to randomly sample one SNP per RAD locus to minimize the effect of LD.

``` bash
vcftools --vcf ../../../merged-datasets.vcf --keep stimpson-inds-to-keep.txt --max-missing 0.9 --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out stimpson-only-maxmissing9-biallelic
```

``` bash
python single_snp.py stimpson-only-maxmissing9-biallelic.recode.vcf 
```

That leaves 7096 SNPs in the Perdido-only dataset.

Then use [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) (Java
applet) to convert that final VCF to Structure format.

Next, import the STR file, check individual-level missingness, re-filter
as needed.

``` r
# Import the allele data for best inds
genepop <- read.structure(file = "mod_stimpson-only-maxmissing9-biallelic.recode.vcf.str",
    n.ind = 10, n.loc = 7096, col.lab = 1, col.pop = 2, row.marknames = 1, col.others = NA,
    onerowperind = FALSE, ask = TRUE)

# Recheck missingness of individuals
missingness <- apply(genepop$tab, 1, function(x) sum(is.na(x))/ncol(genepop$tab))

# How many inds have less than 20% missing data?
length(names(missingness[missingness < 0.2]))  # all!
```

Finally, run a PCA on the final conStruct dataset to make sure neither
PC1 nor PC2 correlate with the level of missing data.

``` r
# Scale NAs by mean value, then run PCA
genepop.scaled <- tab(genepop, freq = TRUE, NA.method = "mean")
first.pca <- dudi.pca(genepop.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
    nf = 10)

# Calculate correlation between PC1 and PC2, by regression
summary(lm(first.pca$li[, 1] ~ missingness))  # P = 0.664
summary(lm(first.pca$li[, 2] ~ missingness))  # P = 0.06
```

## Solon Dixon only —–

``` r
library(tidyverse)

setwd("~/Documents/NorthCarolina/TangledBank/gopher-tortoise/data-analysis/construct/new-way-of-subsetting/sd-only/")

# Filter metadata, keep all the individuals
high.quality.inds.to.keep <- metadata %>%
  filter(Site.fixed.simplified == "SolonDixon") %>%
  filter(!is.na(lat)) %>% # Make sure they also have lat/long data
  # group_by(Site.fixed.simplified) %>%
  # top_n(10, loci) %>%
  dplyr::select(vcfname)

write.table(high.quality.inds.to.keep, "sd-inds-to-keep.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

From the vcf output by ipyrad, filter sites heavily to only include
biallelic sites present in 90% of individuals. Then use this script from
[Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py)
to randomly sample one SNP per RAD locus to minimize the effect of LD.

``` bash
vcftools --vcf ../../../merged-datasets.vcf --keep sd-inds-to-keep.txt --max-missing 0.9 --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out sd-only-maxmissing9-biallelic
```

``` bash
python ~/Documents/NorthCarolina/TangledBank/tbc-hub/3rad-analyses/single_snp.py sd-only-maxmissing9-biallelic.recode.vcf 
```

That leaves 6659 SNPs in the Perdido-only dataset.

Then use [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) (Java
applet) to convert that final VCF to Structure format.

Next, import the STR file, check individual-level missingness, re-filter
as needed.

``` r
# Import the allele data for best inds
genepop <- read.structure(file = "mod_sd-only-maxmissing9-biallelic.recode.vcf.str",
    n.ind = 20, n.loc = 6659, col.lab = 1, col.pop = 2, row.marknames = 1, col.others = NA,
    onerowperind = FALSE, ask = TRUE)

# Recheck missingness of individuals
missingness <- apply(genepop$tab, 1, function(x) sum(is.na(x))/ncol(genepop$tab))

# How many inds have less than 20% missing data?
length(names(missingness[missingness < 0.2]))  # all!
```

Finally, run a PCA on the final conStruct dataset to make sure neither
PC1 nor PC2 correlate with the level of missing data.

``` r
# Scale NAs by mean value, then run PCA
genepop.scaled <- tab(genepop, freq = TRUE, NA.method = "mean")
first.pca <- dudi.pca(genepop.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
    nf = 10)

# Calculate correlation between PC1 and PC2, by regression
summary(lm(first.pca$li[, 1] ~ missingness))  # P = 0.636
summary(lm(first.pca$li[, 2] ~ missingness))  # P = 0.045, very very close, and that's ALL inds! Gets worse as you subset
```

## Conecuh and Solon Dixon —–

``` r
library(tidyverse)

setwd("~/Documents/NorthCarolina/TangledBank/gopher-tortoise/data-analysis/construct/new-way-of-subsetting/conecuh-sd/")

# Filter metadata, keeping top 7 inds per site.fixed (had to keep top 7 to keep dataset positive definite)
high.quality.inds.to.keep <- metadata %>%
  filter(site.new.names == "Conecuh - Other Site" |
         site.new.names == "Conecuh - Site 1" |
         site.new.names == "Conecuh - Site 2" |
         site.new.names == "Conecuh - Site 3" |
         site.new.names == "Conecuh - Site 4" |
         site.new.names == "Conecuh - Site 5" |
         site.new.names == "Conecuh - Site 6" |
         site.new.names == "SolonDixon") %>%
  filter(!is.na(lat)) %>% # Make sure they also have lat/long data
  group_by(site.new.names) %>% 
  top_n(7, loci) %>%
  dplyr::select(vcfname)

write.table(high.quality.inds.to.keep[,2], "conecuh-sd-inds-to-keep.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

From the vcf output by ipyrad, filter sites heavily to only include
biallelic sites present in 90% of individuals. Then use this script from
[Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py)
to randomly sample one SNP per RAD locus to minimize the effect of LD.

``` bash
vcftools --vcf ../../data/merged-datasets.vcf --keep conecuh-sd-inds-to-keep.txt --max-missing 0.9 --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all --out conecuh-sd-maxmissing9-biallelic
```

``` bash
python single_snp.py conecuh-sd-maxmissing9-biallelic.recode.vcf 
```

That leaves 4280 SNPs from 56 inds in the Conecuh+SD dataset.

Then use [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) (Java
applet) to convert that final VCF to Structure format.

Next, import the STR file, check individual-level missingness, re-filter
as needed.

``` r
library(adegenet)
# Import the allele data for best inds
genepop <- read.structure(file = "mod_conecuh-sd-maxmissing9-biallelic.recode.vcf.str",
    n.ind = 56, n.loc = 4280, col.lab = 1, col.pop = 2, row.marknames = 1, col.others = NA,
    onerowperind = FALSE, ask = TRUE)

# Recheck missingness of individuals
missingness <- apply(genepop$tab, 1, function(x) sum(is.na(x))/ncol(genepop$tab))

# How many inds have less than 20% missing data?
length(names(missingness[missingness < 0.2]))  # 50 of 56

# Do they come from all populations?
metadata %>%
    filter(vcfname %in% names(missingness[missingness < 0.2])) %>%
    count(Site.fixed.simplified)

# Save the final list of individuals
final.inds.list <- metadata %>%
    filter(vcfname %in% names(missingness[missingness < 0.2])) %>%
    select(vcfname)

write.table(final.inds.list, "final-inds-list-for-construct-conecuh-sd.txt", col.names = FALSE,
    row.names = FALSE, quote = FALSE)
```

Finally, run a PCA on the final conStruct dataset to make sure neither
PC1 nor PC2 correlate with the level of missing data.

``` r
# Filter original datset to just the final inds
genepop.filtered <- genepop[indNames(genepop) %in% final.inds.list$vcfname]

# Recalculate missingness
missingness <- apply(genepop.filtered$tab, 1, function(x) sum(is.na(x))/ncol(genepop.filtered$tab))

# Scale NAs by mean value, then run PCA
genepop.scaled <- tab(genepop.filtered, freq = TRUE, NA.method = "mean")
first.pca <- dudi.pca(genepop.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
    nf = 10)

# Calculate correlation between PC1 and PC2, by regression
summary(lm(first.pca$li[, 1] ~ missingness))  # P = 0.675
summary(lm(first.pca$li[, 2] ~ missingness))  # P = 0.460
```

[^1]: <alex@tbconservation.org>
