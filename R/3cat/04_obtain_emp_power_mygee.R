# this file was run on ACCRE cluster
## output: results/3cat/rankICCs_from_superpop.csv

library(dplyr)
library(rms)
library(readr)
library(ordinal)
library(bridgedist)

invisible( lapply( list.files("../functions",
                              pattern = "\\.R$",
                              full.names = TRUE), source ) )
