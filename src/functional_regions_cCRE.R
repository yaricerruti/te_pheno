#!/usr/bin/Rscript

args = commandArgs(TRUE)
infile = args[1] 
outfile = args[2]
inclassif = args[3]
inclass = args[4]

packages = c("GenomicRanges", "readr")

for(package in packages)
  library(package, character.only = T)

## removing scientific notation
options(scipen = 100)

########################
## DEFINING FUNCTIONS
########################

range.reduce <- function(df, select_classification, select_class){
  ############# FUNC DESCRIPTION
  ## df : dataframe in the form of a bed file with 4 columns: chrom, start, end, metadata column
  ## colnames are standardized to GRanges default values ["seqnames", "start", "end", "metadata-specific name"]
  ## ranges are subsetted, in this case for strong enhancer selection,
  ## then ranges on sex chromosomes and mitochondrial chromosome are removed
  ## then ordered according to chrom, start end,
  ## then transformed into a GRobject for reduction
  ## a dataframe in the form of the input bed file is returned with reduced ranges
  ############# FUNC DESCRIPTION
  
  "%notin%" <- Negate("%in%")
  df <- df[df[, select_classification] %in% select_class, ]
  df <- df[df$chrom %notin% c("chrX", "chrY", "chrMT"), ]
  df <- df[order(df$chrom, df$chromStart, df$chromEnd), ]
  df <- makeGRangesFromDataFrame(df, seqnames.field = "chrom",
                                 start.field = "chromStart",
                                 end.field = "chromEnd",
                                 keep.extra.columns = F,
                                 ignore.strand = T,
                                 starts.in.df.are.0based = T)
  df <- as.data.frame(reduce(df)) ## !!! RANGES ARE NOW 1-BASED !!!
  df$start=df$start-1 ## !!! RANGES ARE NOW 0-BASED !!!
  return(df[, -c(4, 5)])
}

#######################
## READING DATA
#######################

## reading file
## debug -> CCRE <- readr::read_tsv(file = "/mnt/filippo/prj/te_enhancer_conservation/dataset/v14/encodeCcreCombined.bed", col_names = F)
CCRE <- readr::read_tsv(file = infile, col_names = F)

#######################
## MAIN
#######################
CCRE <- as.data.frame(CCRE)
colnames(CCRE) = c("chrom",	"chromStart",	"chromEnd",	"name",	"score", "strand",
                   "thickStart", "thickEnd", "reserved", "ccre", "encodeLabel",
                   "zScore", "ucscLabel", "accessionLabel",	"description")

## reducing ranges
classification = as.character(inclassif) ## must be a string
classes = c(as.character(inclass)) ## must be a vector
CCRE <- range.reduce(df=CCRE,
                     select_classification = classification,
                     select_class=classes) ## !!! RANGES MUST BE 0-BASED IN OUTPUT!!!

## renaming chromosomes to standard .bed
CCRE$seqnames = gsub("chr", "", CCRE$seqnames)

## output
readr::write_tsv(CCRE, file=outfile, col_names = F, quote = "none")

##
