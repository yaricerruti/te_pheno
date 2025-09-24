#!/usr/bin/Rscript

args = commandArgs(TRUE)
infile_phylop = args[1]
infile_tags = args[2]
outfile = args[3]
w = args[4]

packages = c("readr")

for(package in packages)
  library(package, character.only = T)

################
## READING DATA
###############

## debug -> scores <- read_tsv(file = "/mnt/filippo/prj/te_enhancer_conservation/predictors_dataset/v1/phyloP/dELS_real_ID_phyloP.tab", col_names = F)
scores <- readr::read_tsv(file = infile_phylop, col_names = F)
## debug -> tags <- read_tsv(file = "/mnt/filippo/prj/te_enhancer_conservation/predictors_dataset/v1/overlap_output/dELS_real_te_func_ovl_tags.bed", col_names = F)
tags <- readr::read_tsv(file = infile_tags, col_names = F)

#############
## MAIN
############

## transforming into dataframes and renaming columns for subsequent merge
scores = as.data.frame(scores)
tags = as.data.frame(tags)

colnames(scores) = c("ID", "width", "width.covered", "phyloP.sum", "phyloP.mean0", "phyloP.mean")
colnames(tags) = c("ID", "TAG")

## merging the phyloP and tags dataframes
merged <- merge(scores, tags, by="ID")
merged$ID = NULL
merged = merged[, c(6, 1, 2, 3, 4, 5)]
merged$GROUP = gsub("[0-9]", "", unlist(strsplit(w, "_"))[2])

readr::write_tsv(merged, file = outfile, col_names = F, quote = "none")

##