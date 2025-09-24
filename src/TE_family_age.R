#!/usr/bin/Rscript

args = commandArgs(T)

packages = c("readr")

for(package in packages){
  library(package, character.only = T)}


############
## VARIABLES
############

infile = args[1]
outfile = args[2]
# infile = "/mnt/filippo/prj/te_enhancer_conservation/age_dataset/v2/rmsk_age_processed.bed"
# outfile = "/mnt/filippo/prj/te_enhancer_conservation/age_dataset/v2/te_family_age_test.tsv"


###############
## READING DATA
###############

bed = readr::read_tsv(file = infile, col_names = F)
bed = as.data.frame(bed)
colnames(bed) = c("chr", "start", "end", "class", "age")

#######
## MAIN
#######

bed_list = split(bed, f=bed$class)

## computing family-wise median age
median_age = as.data.frame(sapply(bed_list, function(x) return(median(x$age))))
median_age$CLASS = row.names(median_age)
row.names(median_age) = NULL
median_age = median_age[, c(2, 1)]
colnames(median_age) = c("CLASS", "median_age")

## computing family-wise 3rd quartile of age distribution
quartile_age = as.data.frame(sapply(bed_list, function(x) return(quantile(x$age)[4])))
quartile_age$CLASS = row.names(quartile_age)
row.names(quartile_age) = NULL
quartile_age = quartile_age[, c(2, 1)]
colnames(quartile_age) = c("CLASS", "q3_age")
quartile_age$CLASS = sapply(strsplit(quartile_age$CLASS, "\\."), getElement, 1)

## computing family-wise 90th percentile of age distribution
percentile_age = as.data.frame(sapply(bed_list, function(x) return(quantile(x$age, probs = seq(0, 1, 0.1))[10])))
percentile_age$CLASS = row.names(percentile_age)
row.names(percentile_age) = NULL
percentile_age = percentile_age[, c(2, 1)]
colnames(percentile_age) = c("CLASS", "p90_age")
percentile_age$CLASS = sapply(strsplit(percentile_age$CLASS, "\\."), getElement, 1)

## computing mode (most frequent value) of age distribution
mode_age = as.data.frame(sapply(bed_list, function(x) return(as.numeric(names(sort(table(x$age), decreasing = T)[1])))))
mode_age$CLASS = row.names(mode_age)
row.names(mode_age) = NULL
mode_age = mode_age[, c(2, 1)]
colnames(mode_age) = c("CLASS", "mode_age")

## merging the two results
TE_family_age = merge(median_age, quartile_age, by = "CLASS")
TE_family_age = merge(TE_family_age, percentile_age, by = "CLASS")
TE_family_age = merge(TE_family_age, mode_age, by = "CLASS")
TE_family_age = TE_family_age[order(TE_family_age$CLASS), ]

## rounding
TE_family_age$median_age = round(TE_family_age$median_age, digits = 0)
TE_family_age$q3_age = round(TE_family_age$q3_age, digits = 0)
TE_family_age$p90_age = round(TE_family_age$p90_age, digits = 0)

## output
write.table(TE_family_age, file = outfile, col.names = T, row.names = F, quote = F, sep = "\t")

##