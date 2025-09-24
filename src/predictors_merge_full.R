#!/usr/bin/Rscript

args = commandArgs(TRUE)

outfile = args[1]
infile_rmsk = args[2]
infile_clade = args[3]
infile_ratio = args[4]
infile_score = args[5]
infile_age = args[6]
infile_multiple_pred = args[6:length(args)]

packages = c("readr")

for(p in packages)
  library(p, character.only = T)

## removing scientific notation
options(scipen = 100)

###############
## READING DATA
###############

## reading repeatmasker
rmsk = readr::read_tsv(file = infile_rmsk, col_names = F)
colnames(rmsk) <- c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                    "genoName", "genoStart", "genoEnd", "genoLeft", "strand",
                    "repName", "repClass", "repFamily", "repStart", "repEnd",
                    "repLeft", "id")

## reading clade
clade = read.table(file = infile_clade, header = T)

## reading chromatin accessibility
ratio = read.table(file = infile_ratio, header = T)

## reading score
score = read.table(file = infile_score, header = T)

## reading age
age = read.table(file = infile_age, header = T)

## reading phyloP
infile_phyloP = infile_multiple_pred[grepl("phyloP", infile_multiple_pred)]

phylo = data.frame()
for(f in infile_phyloP) {
  temp = read.table(file = f, header = T)
  temp$FUNC = unlist(strsplit(f, "_"))[2]
  phylo = rbind(phylo, temp)
}
rm(temp)

## reading DRS
infile_DRS = infile_multiple_pred[grepl("DRS", infile_multiple_pred)]

DRS = data.frame()
for(f in infile_DRS) {
  temp = read.table(file = f, header = T)
  temp$FUNC = unlist(strsplit(f, "_"))[2]
  DRS = rbind(DRS, temp)
}
rm(temp)

#######
## MAIN
#######

# merging chromatin accessibility and age
m1 = merge(ratio[, c("CLASS", "coverage", "ovl_bases", "base_ratio")], age[, c("CLASS", "median_age")], by = "CLASS", all = T)

# merging clades
m2 = merge(m1, clade, by.x = "median_age", by.y = "age", all.x = T)

## merging sequence score
m3 = merge(m2, score, by = "CLASS", all = T)

## merging DRS
m4 = merge(m3, DRS[, c("GROUP", "CLASS", "FUNC", "TYPE", "median_DRS")], by = "CLASS", all = T)

## merging phyloP
m5 = merge(m4, phylo[, c("GROUP", "CLASS", "FUNC", "TYPE", "phyloP.class.average")], by = c("GROUP", "CLASS", "FUNC", "TYPE"), all = T)

## merging repeatmasker repClass and repFamily info
rmsk$repName = gsub("_", "-", rmsk$repName)
rmsk = rmsk[rmsk$genoName %in% paste("chr", 1:22, sep = ""), ]
m5 = merge(m5, unique(rmsk[, c("repName", "repClass")]), by.x = "CLASS", by.y = "repName", all.x = T)
m5 = merge(m5, unique(rmsk[, c("repName", "repFamily")]), by.x = "CLASS", by.y = "repName", all.x = T)

## predictors matrix
predictors = m5[, c("GROUP", "FUNC", "TYPE", "CLASS", "repFamily", "repClass", "median_age", "clade", "coverage", "ovl_bases", "base_ratio", "sequence_score", "seq", "median_DRS", "phyloP.class.average")]
colnames(predictors)[colnames(predictors) == "seq"] = "sequence_type"
colnames(predictors)[colnames(predictors) == "CLASS"] = "repName"
colnames(predictors)[colnames(predictors) == "phyloP.class.average"] = "phyloP"
colnames(predictors)[colnames(predictors) == "median_DRS"] = "DRS"

## selecting duplicated because of ambiguous nomenclature
dup = predictors[duplicated(predictors$repName) | duplicated(predictors$repName, fromLast = T), ]  

tokeep = data.frame()
for(n in unique(dup$repName)) {
## keeping repfamily annotation that has the most copies
  keep = names(sort(sapply(split(rmsk[rmsk$repName == n, ], f=rmsk[rmsk$repName == n, ]$repFamily), function(x) return(nrow(x))), decreasing = T))[1]
  tokeep = rbind(tokeep, dup[dup$repName == n & dup$repFamily == keep, ])
}

## removing duplicates and adding correct classes
predictors = predictors[!(duplicated(predictors$repName) | duplicated(predictors$repName, fromLast = T)), ]  
predictors = rbind(predictors, tokeep)


## output
write.table(predictors, file = outfile, col.names = T, row.names = F, quote = F, sep = "\t")

##

