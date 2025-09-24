#!/usr/bin/Rscript

args=commandArgs(T)

packages = c("readr")

for(package in packages){
  library(package, character.only = T)}

options(scipen = 999)

############
## FUNCTIONS
############

"%notin%" = Negate("%in%")

############
## VARIABLES
############

infile_rmsk_prep = args[1]
infile_bed = args[2]
outfile = args[3]


###############
## READING DATA
###############

## reading bed files
bed = readr::read_tsv(file = infile_bed, col_names = F)
## debug -> bed = readr::read_tsv(file = "/mnt/filippo/prj/te_enhancer_conservation/revision_dataset/v1/chrom_acc/rmsk_ccre_overlap.bed" , col_names = F)
bed = as.data.frame(bed)
colnames(bed) = c("chr", "start", "end", "CLASS", "ovl_length")

## reading preprocessed repeatmasker
rmsk_prep = readr::read_tsv(file = infile_rmsk_prep, col_names = F)
## debug -> rmsk_prep = readr::read_tsv(file = "/mnt/filippo/prj/te_enhancer_conservation/revision_dataset/v1/preprocessing/rmsk_filtered.bed", col_names = F)
rmsk_prep = as.data.frame(rmsk_prep)
colnames(rmsk_prep) <- c("chr", "start", "end", "repName")

#######
#######
## MAIN
#######
#######

## computing genomic coverage and count for each repName
rmsk_prep$width = rmsk_prep$end-rmsk_prep$start
repName_cov = sapply(split(rmsk_prep, f=rmsk_prep$repName), function(x) return(sum(x$width)))
repName_count = sapply(split(rmsk_prep, f=rmsk_prep$repName), function(x) return(nrow(x)))
cov_count = merge(as.data.frame(repName_cov), as.data.frame(repName_count), by=0)
colnames(cov_count)[1] = "CLASS"

## preparing bed data for ratio computation (coordinates are not necessary)
bed = bed[, c("CLASS", "ovl_length")]

## computing OVL/total ratio:
## bed is split according to functional region, in each functional region the
## number of bases overlapping a regulatory element is divided by the total number
## of bases (coverage) of a TE family, separately for each family

ratio = data.frame()
for(c in sort(unique(bed$CLASS))) {
  coverage = cov_count[cov_count$CLASS == c, "repName_cov"]
  total_count = cov_count[cov_count$CLASS == c, "repName_count"]
  ovl_bases = sum(bed[bed$CLASS == c, ]$ovl_length)
  ovl_count = nrow(bed[bed$CLASS == c & bed$ovl_length > 0, ])
  ratio = rbind(ratio, data.frame("CLASS" = c,
                                  "coverage" = coverage,
                                  "ovl_bases" = ovl_bases,
                                  "base_ratio" = ovl_bases/coverage,
                                  "total_count" = total_count,
                                  "ovl_count" = ovl_count,
                                  "count_ratio" = ovl_count/total_count))
}


## output
ratio = ratio[order(ratio$CLASS), ]
write.table(ratio, file = outfile, col.names = T, row.names = F, quote = F, sep = "\t")

##
