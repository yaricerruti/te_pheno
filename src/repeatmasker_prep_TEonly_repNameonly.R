#!/usr/bin/Rscript

args = commandArgs(TRUE)

# infile = snakemake@input[[1]]
# outfile = snakemake@output[[1]]
infile = args[1]
outfile = args[2]

packages = c("reshape2", "readr")

for(package in packages)
  library(package, character.only = T)

## removing scientific notation
options(scipen = 100)

########################
## DEFINING FUNCTIONS
########################

###
rep.filter <- function(df_melt){
  ############# FUNC DESCRIPTION
  ## func only selects general classes width total
  ## genome-wide length > 100kbps (width is calculated within the scope of func)
  ## df_melt : df in long format for each rmsk classification type
  ## returns a melt df
  ############# FUNC DESCRIPTION
  
  df_melt$width = df_melt$genoEnd-df_melt$genoStart ## calculating TE width !!! RANGES MUST BE 0 BASED !!!
  cw_df = df_melt[, c(4, 5)] ## only selecting dataframe cols with class names and width
  width_by_classname = tapply(cw_df$width, cw_df$repName, FUN=sum) ## class-wise width sum
  width_by_classname = width_by_classname[width_by_classname>100000] ## only selecting classes over 100kbps genome wide (named vector)
  ## return(names(width_by_classname))
  
  df_melt = df_melt[df_melt[, "repName"] %in% names(width_by_classname), ]
  return(df_melt[, -c(6)]) ## returning filtered repeatmasker removing width column
}


###############
## READING DATA
###############

## reading repeatmasker
## debug -> rmsk <- readr::read_tsv(file = "/mnt/filippo/prj/te_enhancer_conservation/local/data/rmsk_hg38.txt", col_names=F)
rmsk <- readr::read_tsv(file = infile, col_names=F) ## !!! REPEATMASKER RANGES ARE 0-BASED !!!
rmsk = as.data.frame(rmsk)

#######
## MAIN
#######

## repeatmasker data manipulation - removing useless columns, reordering columns
colnames(rmsk) <- c("bin", "swScore", "milliDiv", "milliDel", "milliIns",
                    "genoName", "genoStart", "genoEnd", "genoLeft", "strand",
                    "repName", "repClass", "repFamily", "repStart", "repEnd",
                    "repLeft", "id")
rmsk <- rmsk[, c("genoName", "genoStart", "genoEnd", "repName", "repFamily", "repClass")]

## selecting autosomes
rmsk <- rmsk[rmsk[, "genoName"] %in% paste0("chr", 1:22), ] ## must be a classic dataframe, tibble is extremely slow for some reason

## selecting only repClass classes that refer to annotated TRANSPOSONS, removing otherwise non-transposon repetitive DNA
rmsk <- rmsk[rmsk$repClass %in% c("LINE",
                                  "SINE",
                                  "LTR",
                                  "Retroposon",
                                  "DNA"), ]

## keeping only repName annotation
rmsk = rmsk[, c("genoName", "genoStart", "genoEnd", "repName")]

######################################
## !!! REMOVING 100KBPS COVERAGE FILTER
## since now TajD is computed for a TE class in a genome-wide manner,
## also TE classes covering < 100kpbs may have enough SNPs for TajD
## genome-wide computation (no problem for phyloP or DRS computation)
## only selecting general classes width total genome-wide length > 100kbps
## rmsk_filt <- rep.filter(rmsk)

## renaming TE species (transforming "_" in "-")
rmsk$repName = gsub("_", "-", rmsk$repName)

## renaming chromosomes
rmsk$genoName = gsub("chr", "", rmsk$genoName)
rmsk$genoName = as.integer(rmsk$genoName)

## preparing ranges for GRanges operations:
## renaming columns to GRanges standard
colnames(rmsk)[1:3] <- c("seqnames", "start", "end")

## output (0-based)
readr::write_tsv(rmsk, file=outfile, col_names = F, quote = "none")

##