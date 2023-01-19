Run and plot sPCA
================

This tutorial runs an sPCA on all high quality Gopher Tortoises. By
starting with a different subset of individuals it can be run on any
subsets from the paper. In general this tutorial follows the standard
tutorial from
(adegenet)\[<https://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf>\].

``` r
# Filter original VCF fsing 9 max missing, biallelic, and one snp per locus in
# VCFTools (and using single_snp.py), then convert to STR using PGD Spider

system("vcftools --vcf ../data/merged-datasets.vcf --max-alleles 2 --min-alleles 2 --max-missing 0.9 --recode --recode-INFO-all --out biallelic-9maxmissing-allinds")

system("python ../fastSTRUCTURE/single_snp.py biallelic-9maxmissing-allinds.recode.vcf")

# Leaves 336 individuals with 548 SNPs
```

Next, convert this VCF to STR using [PGD
Spider](http://www.cmpg.unibe.ch/software/PGDSpider/).

Then use `adegenet` to run the sPCA

``` r
library(adegenet)
library(tidyverse)
library(paletteer)

# Import STR
allinds.genepop <- read.structure(file = "mod_biallelic-9maxmissing-allinds.recode.vcf.str",
    n.ind = 336, n.loc = 548, col.lab = 1, col.pop = 2, row.marknames = 1, col.others = NA,
    onerowperind = FALSE, ask = TRUE)

# Import lat/long data for each individual
ind.names <- indNames(allinds.genepop)

# Read in metadata
metadata <- read_csv("../data/modified-metadata-for-seqd-individuals.csv")


# Subset metadata to only inds in structurefile, then keep only individuals
# with lat and long
final.lat.longs <- metadata[match(ind.names, metadata$vcfname), ] %>%
    filter(!is.na(lat))

final.genepop <- allinds.genepop[indNames(allinds.genepop) %in% final.lat.longs$vcfname]

# Select only lat and long from metadata
final.lat.longs <- final.lat.longs %>%
    dplyr::select(lat, long)


# Plot Ho vs He to see if there are alleles that really vary from HWE
allinds.smry <- summary(final.genepop)

plot(allinds.smry$Hobs, allinds.smry$Hexp, xlab = "Ho", ylab = "He")
abline(0, 1, col = "red")  # Looks good!!
```

Next, make sure there isn’t a strong correlation between PC axes and
missing data. If there is, you’ll have to subset the data further to
remove individuals with high missingness.

``` r
# First replace missing values with mean allele freqs,
allinds.scaled <- tab(final.genepop, freq = TRUE, NA.method = "mean")

# Calculate the first 30 axes of the PCA
first.pca <- dudi.pca(allinds.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
    nf = 30)

# See if PC1 or PC2 correlate with missingness
missingness <- apply(final.genepop$tab, 1, function(x) sum(is.na(x))/ncol(final.genepop$tab))
hist(missingness)  # Most inds have very little, a few inds have a lot

summary(lm(first.pca$li[, 1] ~ missingness))  # Not significant
cor(y = first.pca$li[, 1], x = missingness)  # 0.027
```

Finally, having verified that the data is high quality, we can run and
plot the sPCA.

``` r
## Run sPCA
lat.longs <- dplyr::select(final.lat.longs, lat, long)

# Using genetic distance as connection network (type = 5) Unscaled with full
# dataset. d2 = 0.05 = ~5 km linear distance, a high estimate of dispersal for
# Gopher Tortoises
myspca <- spca(final.genepop, xy = lat.longs, scale = TRUE, type = 5, d1 = 0, d2 = 0.05)  # Keep first 20 axes of global and local axes

# Get variances for first two global axes
summary.spca(myspca)
# 5.43% from Axis 1 2.56% from Axis 2

# Attach metadata to the lagged scores in space ($ls)
spca.ls <- myspca$ls
spca.ls$vcfname <- rownames(myspca$ls)

spca.ls.final <- dplyr::select(pc.with.metadata, lat.x, long.x, vcfname, Site.fixed,
    Location) %>%
    left_join(spca.ls, ., by = "vcfname")

# Plot two global spca axes, colored by Location
allinds.spca <- ggplot(data = spca.ls.final) + geom_point(aes(x = `Axis 1`, y = `Axis 2`,
    color = Location), size = 4) + xlab("sPC1 (5.43% variance explained)") + ylab("sPC2 (2.56% variance explained)") +
    scale_color_manual(values = paletteer_d("ggsci::nrc_npg"), name = "Location") +
    theme_bw(base_size = 12)
```
