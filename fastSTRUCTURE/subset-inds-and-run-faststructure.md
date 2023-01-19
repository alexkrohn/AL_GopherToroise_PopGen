Create Input files for fastSTRUCTURE
================

## Subset metadata

First, subset the metadata to each population grouping (found by running
all individuals through fastSTRUCTURE).

``` r
library(tidyverse)

# Import overall metadata ----
metadata <- read_csv("../data/modified-metadata-for-seqd-individuals-new-names.csv")

# All individuals
all.inds <- metadata %>%
    select(vcfname)

write.table(all.inds, "allinds-inds-to-keep.txt", col.names = FALSE, row.names = FALSE,
    quote = FALSE)

# Conecuh + Geneva + Rayonier ----
conecuh_geneva_rayonier <- metadata %>%
    filter(pca.pops == "Conecuh + Solon Dixon" | pca.pops == "Geneva + Rayonier") %>%
    select(vcfname)

write.table(conecuh_geneva_rayonier, "conecuh-genva-rayonier-inds-to-keep.txt", col.names = FALSE,
    row.names = FALSE, quote = FALSE)

# Perdido + Stimpson ----
perdido_stimpson <- metadata %>%
    filter(pca.pops == "Perdido + Stimpson") %>%
    select(vcfname)

write.table(perdido_stimpson, "perdido-stimpson-inds-to-keep.txt", col.names = FALSE,
    row.names = FALSE, quote = FALSE)

# Just conecuh ----
conecuh <- metadata %>%
    filter(pca.pops == "Conecuh + Solon Dixon") %>%
    select(vcfname)

write.table(conecuh, "conecuh-solondixon-inds-to-keep.txt", col.names = FALSE, row.names = FALSE,
    quote = FALSE)

# Geneva/Rayonier ----
geneva_rayonier <- metadata %>%
    filter(pca.pops == "Geneva + Rayonier") %>%
    dplyr::select(vcfname)

write.table(geneva_rayonier, "geneva-rayonier-inds-to-keep.txt", col.names = FALSE,
    row.names = FALSE, quote = FALSE)
```

## Subset the VCF

Next, for each txt file that you just created, use
[VCFtools](https://vcftools.sourceforge.net/) to select SNPs present in
at least 80% of individuals, that are biallelic.

``` bash
for i in $(ls -l *.txt) ;
do vcftools --vcf ../data/merged-datasets.vcf --keep "$i" --max-missing 0.8 --max-alleles 2 --min-alleles 2 --recode-INFO-all --recode --out "$i" ; 
done
```

Finally, select one SNP per RAD locus using this script from [Ben
Anderson](https://github.com/bmichanderson/RAD_scripts/blob/6a61611c3db3ec08ca572e7975b3fd9a0b941ebe/single_snp.py).

``` bash
for i in $(ls -l *recode.vcf) ; 
do python single_snp.py "$i" ; 
done
```

That will leave you with five VCFs to run through fastSTRUCTURE

## Running fastSTRUCTURE

You’ll first need to conver the VCF to a special fastSTRUCTURE formatted
STR file. I used [PGD
Spider](http://www.cmpg.unibe.ch/software/PGDSpider/) for this. For the
purposes of this script, all of the structure inputs will have the same
file name, but with `.str` appended on. For example, the VCF with all
individuals `allinds-inds-to-keep.txt.recode.vcf` will be converted to
`allinds-inds-to-keep.txt.recode.vcf.str`.

After converting, run
[fastSTRUCTURE](https://rajanil.github.io/fastStructure/) on each STR
file, and use the `chooseK.py` script to find the K that maximizes the
marginal likelihood.

For all individuals:

``` bash
mkdir allinds-output
for i in {1..7} ; 
do structure.py -K $i --input=allinds-inds-to-keep.txt.recode.vcf --output=allinds-output/k$i --format=str --prior=logistic; done
chooseK.py --input=allinds-output/k*
```

For all individuals except Perdido and Stimpson (only including Conecuh,
Solon Dixon, Geneva and Rayonier):

``` bash
mkdir conecuh-genva-rayonier-output
for i in {1..7} ; 
do structure.py -K $i --input=conecuh-genva-rayonier-inds-to-keep.txt.recode.vcf --output=conecuh-genva-rayonier-output/k$i --format=str --prior=logistic; done
chooseK.py --input=conecuh-genva-rayonier-output/k*
```

For just individuals from Perdido and Stimpson:

``` bash
mkdir perdido-stimpson-output
for i in {1..4} ; 
do structure.py -K $i --input=perdido-stimpson-inds-to-keep.txt.recode.vcf --output=perdido-stimpson-output/k$i --format=str --prior=logistic; done
chooseK.py --input=perdido-stimpson-output/k*
```

For just individuals from Conecuh and Solon Dixon:

``` bash
mkdir conecuh-solondixon-output
for i in {1..5} ; 
do structure.py -K $i --input=conecuh-solondixon-inds-to-keep.txt.recode.vcf --output=conecuh-solondixon-output/k$i --format=str --prior=logistic; done
chooseK.py --input=conecuh-solondixon-output/k*
```

For just individuals from Geneva and Rayonier:
geneva-rayonier-inds-to-keep.txt

``` bash
mkdir geneva-rayonier-output
for i in {1..4} ; 
do structure.py -K $i --input=geneva-rayonier-inds-to-keep.txt.recode.vcf --output=geneva-rayonier-output/k$i --format=str --prior=logistic; done
chooseK.py --input=geneva-rayonier-output/k*
```
