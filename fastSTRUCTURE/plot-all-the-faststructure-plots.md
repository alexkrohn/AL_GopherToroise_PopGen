Plot Many fastSTRUCTURE Admixture Plots
================

Plot all fastSTRUCTURE admixture plots from one fastSTRUCTURE run. This
script will recreate Supplementary Files 1, but could be used to create
File 3.

This script assumes you already ran fastSTRUCTURE using the
`subset-inds-and-runfastructure.md` tutorial, and now have all inputs
for the all individual run in the `allinds-output` folder.

``` r
library(tidyverse)

# Set working directory where *meanQ files are
setwd("allinds-output/")

# Set min and max k to plot
min.k <- 2
max.k <- 7

# find the files
mean.q.files <- list.files(pattern = paste("*[", min.k, "-", max.k, "].meanQ", sep = ""))

# Import each file as an item in as a list
meanq.list <- map(mean.q.files, ~read_delim(., delim = "  ", col_names = FALSE))


# Import individual names
individual.names <- radiator::extract_individuals_vcf("../mod_allinds-inds-to-keep.txt.recode.vcf") %>%
    rename(vcfname = INDIVIDUALS)

# Import the metadata
pop.df <- read.delim("../data/modified-metadata-for-seqd-individuals-new-names.csv",
    sep = ",", header = TRUE)

# Join metadata to VCF individuals
formatted.results <- left_join(individual.names, pop.df, by = "vcfname")

# Pick the column that should be the population
formatted.results$pop <- formatted.results$Site.fixed.simplified

pop.df.final <- formatted.results %>%
    dplyr::select(vcfname, pop)

# open the pdf file
pdf("allinds-fastructure-plots.pdf")

# Set the margins to be small, but in 1 column of 6 rows
par(mfrow = c(length(meanq.list), 1), mar = c(3.55, 4, 1, 1))

for (k in c(1:length(meanq.list))) {
    # Load the file, and attach metadata (in the same order as
    # individual.names)
    structure.results <- bind_cols(meanq.list[[k]], formatted.results)

    # To add lines separating pops, we need to know how many individuals of
    # each pop we have
    pop.df.final <- structure.results %>%
        dplyr::select(ind, pop)

    # Count the number of inds per pop
    pop.lengths <- count(pop.df.final, pop)

    # Arrange the pops as needed
    pop.lengths.rearranged <- pop.lengths %>%
        arrange(match(pop, c("SolonDixon", "1", "2", "3", "4", "5", "6", "Conecuh Unknown",
            "Geneva", "Rayonier", "Perdido", "Stimpson")))

    # Re order the structure results based on the pop lengths order
    structure.results$pop <- gdata::reorder.factor(pop.df.final$pop, new.order = pop.lengths.rearranged$pop)
    structure.results <- arrange(structure.results, pop)

    # Now re-do pop.lengths :''-D
    pop.lengths.rearranged <- count(structure.results, pop) %>%
        slice(match(unique(structure.results$pop), pop))  # Maintain the same order of pops as in formatted.results (dplyr sorts them alphabetically for some reason)

    # Find the last row number of each pop
    pop.lengths.rearranged$x.value.for.line <- cumsum(pop.lengths.rearranged$n)

    # And a value for 1/2 the distance between each group
    pop.lengths.rearranged$name.position <- pop.lengths.rearranged$x.value.for.line -
        (pop.lengths.rearranged$n)/2

    # Set the colors
    color.vector <- paletteer::paletteer_d("ggsci::default_igv")

    # Plot it all together!
    barplot(as.matrix(t(structure.results[, seq(from = 1, to = k + 1)])), col = color.vector[seq(from = 1,
        to = k + 1), 1], space = 0, border = NA, main = paste("K=", k + 1), ylab = "Admixture")
    text(x = pop.lengths.rearranged$name.position, -0.5, pop.lengths.rearranged$pop,
        xpd = T, srt = 90, cex = 0.8)
    abline(v = pop.lengths.rearranged$x.value.for.line, lty = 5, col = "white")
}
dev.off()
```
