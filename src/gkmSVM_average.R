#!/usr/bin/Rscript

args = commandArgs(T)

packages = c("readr")

for(package in packages) {
  library(package, character.only = T)
}

############
## VARIABLES
############

input_consensus=args[1]
input_instances=args[2]
output=args[3]

###############
## READING DATA
###############

# input_consensus = "/mnt/filippo/prj/te_enhancer_conservation/gkmSVM_dataset/v6/ccre_consensus_SVMmodel_weights.out"
# input_instances = "/mnt/filippo/prj/te_enhancer_conservation/gkmSVM_dataset/v6/ccre_instances_SVMmodel_weights.out"

consensus = readr::read_tsv(file = input_consensus, col_names = F)
consensus = as.data.frame(consensus)
consensus$seq = "consensus"
colnames(consensus) = c("CLASS", "sequence_score", "seq")

instances = readr::read_tsv(file = input_instances, col_names = F)
instances = as.data.frame(instances)
instances$seq = "instances"
colnames(instances) = c("CLASS", "sequence_score", "seq")


#######
#######
## MAIN
#######
#######

instances$CLASS = gsub("_.*", "", instances$CLASS)
avg = as.data.frame(tapply(instances$sequence_score, instances$CLASS, FUN=mean))
avg$CLASS = row.names(avg)
row.names(avg) = NULL
colnames(avg)[1] = c("sequence_score")
avg = avg[, c(2, 1)]
avg$seq = "instances"
out = rbind(consensus, avg)
write.table(out, file = output, col.names = T, row.names = F, quote = F, sep = "\t")

##



