args = commandArgs(TRUE)
infile_bed = args[1]
outfile <- args[2]

library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

debug <- FALSE

bed <- read.delim(infile_bed, header = FALSE)
colnames(bed) <- c("chr", "start", "end", "CLASS")
bed$chr <- paste0("chr", bed$chr)

bed <- split(bed, f = bed$CLASS)
bed <- lapply(bed, makeGRangesFromDataFrame)

if (debug) bed <- tail(bed)

compute_gc <- function(gr) {
    seq <- getSeq(Hsapiens, gr)
    base <- alphabetFrequency(seq, baseOnly = TRUE)
    base <- as.data.frame(base)
    base$total <- apply(base[, c("A", "C", "G", "T", "other")], 1, sum)
    base$gc <- apply(base[, c("C", "G")], 1, sum)
    sum(base$gc)/sum(base$total)
}

gc <- as.data.frame(sapply(bed, compute_gc))
colnames(gc) <- "gc_content"

list <- c("debug", "gc")
save(list = list, file = outfile)

