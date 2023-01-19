Calculate Genetic Diversity Statistics
================

If you haven’t already, be sure to open RStudio from Terminal. This uses
the same VCF as the fastSTRUCTURE analyses, so follow the instructions
in the `subset-inds-and-run-faststructure` tutorial to generate that VCF
if you haven’t already.

``` r
library(hierfstat)
library(vcfR)
library(tidyverse)


# Import pop data
pop.data <- read_delim("../data/modified-metadata-for-seqd-individuals.csv", delim = ",") %>%
    mutate(site.new.names = case_when(Site.fixed == "County Road 24" ~ "Conecuh - Other Site",
        Site.fixed == "FS 305" ~ "Conecuh - Other Site", Site.fixed == "FSR324" ~
            "Conecuh - Other Site", Site.fixed == "Hogfoot Road" ~ "Conecuh - Other Site",
        Site.fixed == "Nellie Pond Road" ~ "Conecuh - Other Site", Site.fixed ==
            "Open Pond Road" ~ "Conecuh - Other Site", Site.fixed == "1" ~ "Conecuh - Site 1",
        Site.fixed == "2" ~ "Conecuh - Site 2", Site.fixed == "3" ~ "Conecuh - Site 3",
        Site.fixed == "4" ~ "Conecuh - Site 4", Site.fixed == "5" ~ "Conecuh - Site 5",
        Site.fixed == "6" ~ "Conecuh - Site 6", TRUE ~ Site.fixed)) %>%
    mutate(pca.pops = case_when(Location == "Conecuh" | Location == "SolonDixon" ~
        "Conecuh + Solon Dixon", Location == "Perdido" | Location == "Stimpson" ~
        "Perdido + Stimpson", Location == "Geneva" | Location == "Rayonier" ~ "Geneva + Rayonier",
        TRUE ~ Location))

pop.data$Site.fixed[is.na(pop.data$Site.fixed)] <- "unknown"
pop.data$Location[is.na(pop.data$Location)] <- "unknown"

# Column with population name
pop.col <- "pca.pops"

# Set VCF
my.vcf <- "../fastSTRUCTURE/allinds-8maxmissing-onesnp.vcf"

# Read VCF
myvcf <- vcfR::read.vcfR(my.vcf)

# Convert to gen ID
genid <- vcfR2genind(myvcf)

# Match pop names to individuals
inds <- data.frame(vcfname = row.names(genid$tab))
pop.names <- left_join(inds, pop.data, by = "vcfname") %>%
    dplyr::select(all_of(pop.col))

# Numbers of samples per population
n <- table(pop.names)

# Convert to hierfstat df
df <- genind2hierfstat(genid, pop = pop.names)

# CALCULATE!
hs <- Hs(df)
ho <- Ho(df)
fis <- 1 - (ho/hs)

# file name
file.name <- "all-torts-gen-diversity"

# Write to disk
write_csv(x = data.frame(pop = names(ho), n = n, ho = round(ho, 4))[, c(1, 3, 4)],
    file = paste(file.name, "-ho-by-pca-pop.csv", sep = ""), col_names = TRUE)
write_csv(x = data.frame(pop = names(hs), n = n, hs = round(hs, 4))[, c(1, 3, 4)],
    file = paste(file.name, "-hs-by-pca-pop.csv", sep = ""), col_names = TRUE)
write_csv(x = data.frame(pop = names(fis), n = n, fis = round(fis, 4))[, c(1, 3,
    4)], file = paste(file.name, "-fis-by-pca-pop.csv", sep = ""), col_names = TRUE)

# final df
write_csv(x = data.frame(pop = names(fis), n = n, ho = round(ho, 4), hs = round(hs,
    4), fis = round(fis, 4))[, c(1, 3, 4:6)], file = paste(file.name, "-pop-gen-stats-final-by-pca-pop.csv",
    sep = ""), col_names = TRUE)

# Kinship derived FST
myvcf <- read.VCF(my.vcf, convert.chr = FALSE)

print(paste("Calculating Fst for ", file.name, ". Started at ", Sys.time(), sep = ""))
fs.dat <- fs.dosage(myvcf, pop = pop.names[, 1])

write.table(round(fs.dat$Fst2x2, 4), file = paste(file.name, "-fst-by-pca-pop.csv",
    sep = ""), col.names = TRUE, quote = FALSE, sep = ",")
```
