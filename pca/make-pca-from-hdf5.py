#!/usr/bin/env python

## Convert HDF5 to PCA using ipyrad
# Inputs are HDF5, population map (with columns of 'pop' and 'ind') as a CSV, and output full path with no file extension

import ipyrad.analysis as ipa
import pandas as pd
import numpy as np
import sys

# Point to dataset
data = input("What HDF5 dataset should we use? ")

# Load data
data = data

popfile = input("What population map should we use? ")

# Import population map, then convert to dictionary
df = pd.read_csv(popfile, sep = ",")

# Sort population map by inds and pops
pops = input("What column in the population map has the population's names? ")
inds = input("What column in the population map has the individual's names? ")

imap = imap = df.groupby(pops)[inds].apply(list).to_dict()

output = input("What is the path to the output (excluding file extensions)? ")

# Set missing data thresholds per population (percent of a pop that must have data)
minmap = {i: 0.1 for i in imap}

# Mincov filters SNPs that are shared across less than 50% of the individuals
pca = ipa.pca(data = data,
             imap = imap,
             minmap = minmap,
             mincov = 0.5,
             impute_method = "sample")

# Run the PCA
print("Running PCA...")

pca.run(nreplicates = 25)

print("Done")

# Save a plot of axes 1 and 2 to get the percents explained
print("Drawing PCA...")

pca.draw(outfile = "%s.pdf" % output)

print("Done")

# Save PCA dataset to csv
print("Saving PCA dataset...")

pd.DataFrame(pca.pcaxes[0], index = pca.names).to_csv("%s.csv" % output)

print("Done!")
