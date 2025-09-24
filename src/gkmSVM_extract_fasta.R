#!/usr/bin/Rscript

args = commandArgs(T)

packages = c("readr",
             "seqinr",
             "GenomicRanges",
             "BSgenome.Hsapiens.UCSC.hg38")

for(package in packages) {
  library(package, character.only = T)
}

############
## VARIABLES
############

## debug -> infile_missing = "/mnt/filippo/prj/te_enhancer_conservation/gkmSVM_dataset/v6/missing.txt"
## debug -> infile_rmsk = "/mnt/filippo/prj/te_enhancer_conservation/local/data/rmsk_hg38.txt"

infile_missing = args[1]
infile_rmsk = args[2]
outfile = args[3]

###############
## READING DATA
###############

missing = scan(file = infile_missing, character())

rmsk = readr::read_tsv(file = infile_rmsk, col_names = F)
rmsk = as.data.frame(rmsk)
rmsk = rmsk[, c(6, 7, 8, 11)]
colnames(rmsk) = c("chr", "start", "end", "CLASS")


#######
#######
## MAIN
#######
#######


rmsk$CLASS = gsub("_", "-", rmsk$CLASS)

instances = rmsk[rmsk$CLASS %in% missing, ]

instances_list = lapply(split(instances, f=instances$CLASS), function(x) {
  
  ## no sampling, testing model on all copies
  # ## sampling 1000 elements
  # x = x[sample(row.names(x), size = as.numeric(params_size)), ]
  
  ## reordering
  x = x[order(x$chr, x$start, x$end), ]
  
  ## transforming to GRanges
  x = makeGRangesFromDataFrame(x,
                               seqnames.field = "chr",
                               start.field = "start",
                               end.field = "end",
                               keep.extra.columns = T,
                               ignore.strand = T,
                               starts.in.df.are.0based = T)
  
})

## retrieving fasta sequences
seqs = list()
for(n in names(instances_list)) {
  
  ## extracting DNA sequences from hg38
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg38, instances_list[[n]])
  
  ## generating list where each element is a TE sequence
  for(i in 1:length(seq)) {
    
    ## transforming DNAstring obj to normal string for fasta output
    seqs[[paste0(n, "_", i)]] = as.character(seq[[i]])
  }
}

seqinr::write.fasta(sequences = seqs,
                    names = names(seqs),
                    file.out = outfile)

##
