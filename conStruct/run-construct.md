Run conStruct on a Subset of Gopher Tortoises
================

## Introduction

Use this file to subset a highly filtered Structure dataset (created
using `choosing-high-quality-inds.Rmd`) to individuals with less than
20% missing data. Then run conStruct.

For this example, we will use the subset that includes individuals from
all sites, but the principal is the same for each of the subsets.

``` r
library(conStruct)
library(tidyverse)

setwd("allinds")

# Import the allele data
construct.data <- structure2conStruct(infile = "mod_all-sites-maxmissing9-biallelic.recode.vcf.str",
    onerowperind = FALSE, start.loci = 3, start.samples = 2, missing.datum = -9,
    outfile = "all-sites-dataset")

# Import lat/long data.
metadata <- read_csv("../../data/metadata/modified-metadata-for-seqd-individuals.csv")

latlongdata <- metadata %>%
    dplyr::select(vcfname, lat, long)

## Find and remove individuals that have > 20% missing data Calculate
## missingness
missingness <- apply(construct.data, 1, function(x) sum(is.na(x))/ncol(construct.data))

# Save the names of individuals to keep
final.inds.list <- metadata %>%
    filter(vcfname %in% names(missingness[missingness < 0.2])) %>%
    select(vcfname)

final.inds.list

# Subset construct data
final.data <- construct.data[row.names(construct.data) %in% final.inds.list$vcfname,
    ]

# Data check. Should leave 62 inds
dim(construct.data)
dim(final.data)



# Subset metadata to just the inds in the final dataset
final.lat.longs <- latlongdata[latlongdata$vcfname %in% rownames(final.data), ]


# Make Coordinate dataframe
coords <- final.lat.longs %>%
    dplyr::select(lat, long) %>%
    as.matrix()

# Create geographic distance matrix
geoDist <- fields::rdist.earth(coords, miles = FALSE)

# Check to make sure all the datasets have the right numbers of individuals
setdiff(final.lat.longs$vcfname, rownames(final.data))
setdiff(rownames(final.data), final.lat.longs$vcfname)
dim(geoDist)
dim(final.data)
```

Run conStruct!

``` r
my.xvals <- x.validation(train.prop = 0.9, n.reps = 10, K = 1:7, freqs = final.data,
    geoDist = geoDist, coords = coords, prefix = "k1-7", n.iter = 10000, make.figs = TRUE,
    save.files = TRUE, parallel = TRUE, n.nodes = 20)
```
