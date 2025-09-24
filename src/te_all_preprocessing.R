#!/usr/bin/Rscript

args = commandArgs(TRUE)

rmsk_infile = args[1]
mhc_infile = args[2]
outfile = args[3]
w = args[4]

packages = c("GenomicRanges", "regioneR", "readr")

for(package in packages)
  library(package, character.only = T)

## removing scientific notation
options(scipen = 100)


########################
## DEFINING FUNCTIONS
########################


###
rep.reduce <- function(df_melt){
  ############# FUNC DESCRIPTION
  ## df_melt : df in long format for each rmsk classification type
  ## class-wise reduction is performed:
  ## the melt dataframe is split into a list where each element is a
  ## dataframe containing ranges of a specific class
  ## each element is transformed into a GRobject, reduced (class-wise reduction)
  ## and returned as a 0-based dataframe
  ############# FUNC DESCRIPTION
  
  df_list = split(df_melt, f=df_melt$class.name) ## each element of the list is a TE class in order to perform a class-wise reduction
  gr_list = lapply(df_list, makeGRangesFromDataFrame, ## each element is transformed into a GRobject
                   seqnames.field = "seqnames",
                   start.field = "start",
                   end.field = "end",
                   ignore.strand = T,
                   keep.extra.columns = T,
                   starts.in.df.are.0based = T)
  gr_list = lapply(gr_list, reduce) ## each element is reduced (class-wise reduction)
  gr_list = lapply(gr_list, as.data.frame) ## back to dataframes !!! RANGES ARE NOW 1-BASED !!!
  for(name in names(gr_list)) ## class names are retained as list names, reassigned to each dataframe
    gr_list[[name]]$class.name = name
  df_list = do.call(rbind, gr_list)
  df_list = df_list[, -c(4, 5)]
  row.names(df_list) = NULL
  df_list$start = df_list$start-1 ## !!! RANGES ARE NOW 0-BASED !!!
  return(df_list) ## final dataframe is returned
}


###
MHC.remove <- function(ranges, to_remove){
  ############# FUNC DESCRIPTION
  ## ranges : a dataframe of ranges
  ## to_remove : a dataframe of coordinates to remove
  ## func takes two sets of ranges as input dataframes
  ## and removes the coordinates submitted as to_remove
  ## used to remove MHC coordinates from TEs and functional regions
  ############# FUNC DESCRIPTION
  
  if("strand" %in% colnames(ranges))
    ranges$strand = NULL
  ranges <- makeGRangesFromDataFrame(ranges,
                                     seqnames.field = "seqnames",
                                     start.field = "start",
                                     end.field = "end",
                                     keep.extra.columns = T,
                                     ignore.strand = T,
                                     starts.in.df.are.0based = T)
  colnames(to_remove) = c("seqnames", "start", "end")
  to_remove <- makeGRangesFromDataFrame(to_remove,
                                        seqnames.field = "seqnames",
                                        start.field = "start",
                                        end.field = "end",
                                        keep.extra.columns = T,
                                        ignore.strand = T,
                                        starts.in.df.are.0based = T) ## MHC ranges are 0-based
  mhc_ovl <- findOverlaps(ranges, to_remove)
  ranges_nomhc <- as.data.frame(ranges[-queryHits(mhc_ovl)])
  ranges_nomhc$strand = NULL
  ranges_nomhc$start = ranges_nomhc$start-1 ## !!! RANGES ARE 0-BASED !!!
  return(ranges_nomhc)
}



###
rep.randomize <- function(df_melt){
  ############# FUNC DESCRIPTION
  ## df_melt : df in long format for each rmsk classification type
  ## func circularly randomizes tes in melt form in a class-wise manner
  ## real tes are split into a list class-wise, trasformed into GRanges
  ## and circularly randomized class-by-class with regioneR::circularRandomizeRegions
  ## and standard parameters (hg19 genome, NA mask)
  ## a dataframe with randomized ranges is returned
  ## +1 is added to all randomized starts in order to avoid obtaining negative coordinates
  ## when retrasforming to 0-based before output
  ############# FUNC DESCRIPTION
  
  df_melt$seqnames = paste0("chr", df_melt$seqnames)
  df_list = split(df_melt, f=df_melt$class.name)
  gr_list = lapply(df_list, makeGRangesFromDataFrame,
                   seqnames.field = "seqnames",
                   start.field = "start",
                   end.field = "end",
                   keep.extra.columns = T,
                   ignore.strand = T,
                   starts.in.df.are.0based = F)
  gr_list = lapply(gr_list, regioneR::circularRandomizeRegions, genome = "hg38", mask = NA)
  gr_list = lapply(gr_list, as.data.frame)
  df_rand = do.call(rbind, gr_list)
  df_rand$seqnames = gsub("chr", "", df_rand$seqnames)
  df_rand = df_rand[df_rand$width > 2, ] ## removing all ranges with width 2 in order to avoid generating a 0-width range...
  df_rand$start=df_rand$start+1 ### !!! +1 to all randomized starts | ...here (width = at least 2, so start + 1 -> width = 1)
  rownames(df_rand) = NULL
  return(df_rand[, -c(4, 5)])
}


###
###############
## READING DATA
###############

## reading repeatmasker
## debug -> rmsk <- readr::read_tsv(file="/mnt/filippo/prj/te_enhancer_conservation/predictors_dataset/v2/rmsk_preprocessed.bed", col_names=F)
rmsk <- readr::read_tsv(file=rmsk_infile, col_names=F) ## !!! REPEATMASKER RANGES ARE 0-BASED !!!

## reading MHC coordinates
## debug -> mhc <- read.delim(file = "/mnt/filippo/prj/te_enhancer_conservation/local/data/MHC.bed", header = F)
mhc <- readr::read_tsv(file=mhc_infile, col_names=F)

#######
## MAIN
#######

## transforming tibbles into dataframes
rmsk <- as.data.frame(rmsk)
mhc <- as.data.frame(mhc)

## colnames
colnames(rmsk) = c("seqnames", "start", "end", "class.name")
colnames(mhc) = c("seqnames", "start", "end")

if(w != "real")
  rmsk <- rep.randomize(rmsk)

## merging overlapping TEs
rmsk_reduced <- rep.reduce(rmsk)

## removing MHC coordinates from TEs
rmsk_reduced <- MHC.remove(rmsk_reduced, mhc)

## adding @ALL tag
rmsk_reduced$class.name = paste0(rmsk_reduced$class.name, "@ALL")

## preparing output
out = rmsk_reduced
out = out[order(out$seqnames, out$start, out$end), ]
row.names(out) = NULL
out$width = NULL

##########################
## !!! OUTPUT IS 0 BASED !!!
##########################
readr::write_tsv(out, file=outfile, append = F, col_names = F, quote = "none")

##