Make IBD Plots and Run Mantel Test
================
Alex Krohn[^1]

## Introduction

For each pairwise comparison of individuals, this plots pairwise genetic
distance by pairwise geographic distance to evaluate the significance of
Isolation by Distance in the data. Then, we test for Isolation by
Distance using a Mantel test.

## Make the IBD Plot

``` r
library(RgoogleMaps)
library(raster)
library(gdata)
library(tidyverse)
library(SNPRelate)

# Load the data
snpgdsVCF2GDS(vcf.fn = ".../data/allinds-0.85maxmissing-onesnp.vcf", out.fn = "allinds-0.85maxmissing-onesnp.pcaVar.gds")
genofile <- snpgdsOpen("allinds-0.85maxmissing-onesnp.pcaVar.gds")
ibs.values <- snpgdsIBS(genofile, num.thread = 2, autosome.only = FALSE)

geno.matrix <- ibs.values$ibs

# Extract the individual's names from the VCF, and match the column name in the
# resulting df to match the column name from the metadata
individual.names <- radiator::extract_individuals_vcf("../data/allinds-0.85maxmissing-onesnp.vcf") %>%
    rename(vcfname = INDIVIDUALS)

# Load lat longs, and subset just to the individuals genotyped
metadata <- read_csv("../modified-metadata-for-seqd-individuals.csv") %>%
    dplyr::select(vcfname, lat, long)

final.lat.longs <- as.data.frame(left_join(individual.names, metadata, by = "vcfname"))

# Set coordinate system
coordinates(final.lat.longs) <- c("long", "lat")
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(final.lat.longs) <- crs.geo

# Create geographic distance matrix. Change the diagonal to be NA instead of 0
geodist <- pointDistance(final.lat.longs, longlat = TRUE)
geodist[which(geodist == 0)] <- NA

# Check to make sure the matrices are the same dimensions
dim(geodist)
dim(gen.distance.matrix)

# Plot
par(mar = c(5.1, 4.7, 4.2, 2.1))
plot(geodist/1000, as.matrix(1 - geno.matrix), ylab = expression(paste("Pairwise Genetic Distance (",
    pi, ")", sep = "")), xlab = "Pairwise Geographic Distance (km)", pch = 19, cex = 1.5,
    cex.lab = 1.3)
```

# Run the Mantel Test

``` r
library(vegan)
vegan::mantel(geodist, geno.matrix, permutations = 999)
```

[^1]: <alex@tbconservation.org>
