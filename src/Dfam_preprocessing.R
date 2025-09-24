#!/usr/bin/Rscript

args=commandArgs(T)

library(seqinr)

############
## VARIABLES
############

infile_seqs = args[1]
infile_ratio = args[2]
outfile_matching = args[3]
outfile_missing = args[4]

###############
## READING DATA
###############

## debug -> dfam = readLines("/mnt/filippo/prj/te_enhancer_conservation/local/data/Dfam-RepeatMasker.lib")
dfam = readLines(infile_seqs)

## debug -> ratio = read.table(file = "/mnt/filippo/prj/te_enhancer_conservation/ratio_dataset/v2/te_base_ovl_ratio.tsv", header = T)
ratio = read.table(infile_ratio, header = T)


#######
#######
## MAIN
#######
#######


########################################################
## generating Dfam consensus list with repNames as names

consensus_list = list()

for(element in dfam){
  if(startsWith(element, ">")) {
    key = gsub("_", "-", unlist(strsplit(gsub("^>", "", element), "#"))[1])
    consensus_list[[key]] = c()
  } else { ## possible since there are no duplicated keys
    consensus_list[[key]] = append(consensus_list[[key]], element)
  }
}

consensus_list = lapply(consensus_list, paste, collapse = "")

perfect_match = sort(intersect(unique(ratio$CLASS), unique(names(consensus_list))))
missing = sort(setdiff(unique(ratio$CLASS), unique(names(consensus_list))))

consensus_list = consensus_list[perfect_match]

## output matching
seqinr::write.fasta(sequences = consensus_list,
                    names = names(consensus_list),
                    file.out = outfile_matching)

## output missing
write.table(as.data.frame(missing), file = outfile_missing, col.names = F, row.names = F, quote = F, sep = "\t")

# ########
# ## PLOTS
# ########
# 
# ################################
# ## sequence length distributions
# ## comparison between consensus and chromHMM enhancers
# 
# library(ggplot2)
# library(ggpubr)
# 
# consensus_length = as.data.frame(sapply(consensus_list, nchar))
# consensus_length$CLASS = row.names(consensus_length)
# colnames(consensus_length) = c("length", "CLASS")
# consensus_length = consensus_length[, c(2, 1)]
# row.names(consensus_length) = NULL
# 
# rmsk = readr::read_tsv(file = "/mnt/filippo/prj/te_enhancer_conservation/local/data/rmsk_hg38.txt", col_names = F)
# rmsk = as.data.frame(rmsk)
# rmsk = unique(rmsk[, c(11, 12)])
# rmsk$X11 = gsub("_", "-", rmsk$X11)
# consensus_length = merge(consensus_length, rmsk, by.x = "CLASS", by.y = "X11", all.x = T)
# colnames(consensus_length)[3] = "repClass"
# table(consensus_length$repClass)
# 
# hmm = read.table(file = "/mnt/filippo/prj/te_enhancer_conservation/gkmSVM_dataset/v1/chromHMM_enhancers_subset.bed")
# hmm$length = hmm$V3-hmm$V2
# 
# histogram_consensus_facet = ggplot(consensus_length, aes(x=length)) +
#   theme_pubr(base_size = 10) +
#   geom_histogram(bins = 100) +
#   facet_wrap(~repClass) +
#   labs(title = "consensus sequences length distribution")
# 
# histogram_consensus_facet
# 
# histogram_consensus_vs_hmm = ggarrange(ggplot(consensus_length, aes(x=length)) +
#                                          theme_pubr(base_size = 10) +
#                                          geom_histogram(bins = 200) +
#                                          scale_x_continuous(limits = c(0, 10000), breaks = seq(0, 10000, 1000)) +
#                                          labs(title = "consensus sequences"),
#                                        ggplot(hmm, aes(x=length)) +
#                                          theme_pubr(base_size = 10) +
#                                          geom_histogram(bins = 200) +
#                                          scale_x_continuous(limits = c(0, 10000)) +
#                                          scale_x_continuous(limits = c(0, 10000), breaks = seq(0, 10000, 1000)) +
#                                          labs(title = "chromHMM enhancers"),
#                                        nrow = 2)
# 
# summary(consensus_length$length)
# summary(hmm$length)

##

