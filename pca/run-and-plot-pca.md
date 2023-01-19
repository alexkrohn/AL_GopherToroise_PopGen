Run and Plot PCA
================

First, run the `make-pca-from-hdf5.py` script on the main
`merged-dataset.vcf`. This is run on all individuals, but can be applied
to any of the subsets. You would only have to change the metadata input
to the python file to a file that contained only the individuals present
in the subset of your choice (see `subset-inds-and-run-fastSTRUCTURE`
tutorial).

``` bash
python make-pca-from-hdf5.py -o all-inds-pca -d ../data/merged-dataset.vcf -m ../data/modified-metadata-for-seqd-individuals.csv -i vcfname -p dummypop
```

Next, plot the results in R.

``` r
library(tidyverse)

# Input PCA data
pca.data <- read.table("all-inds-pca.csv", sep = ",", skip = 1) %>%
    rename(vcfname = V1)


# Import pop data
pop.data <- read_csv("../data/modified-metadata-for-seqd-individuals.csv")

# Join them by ind
fulldataset <- left_join(pca.data, pop.data, by = "vcfname")


# Plot! (Get the % variance from the pdf output from the python script)
allinds.pca <- ggplot(data = fulldataset) + geom_point(aes(x = V2, y = V3, color = Location),
    size = 4) + xlab("PC1 (2.1% variance explained)") + ylab("PC2 (1.2% variance explained)") +
    scale_color_manual(values = paletteer_d("ggsci::nrc_npg"), name = "Location") +
    theme_bw(base_size = 12)
```
