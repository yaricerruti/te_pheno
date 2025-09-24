#!/usr/bin/Rscript

args = commandArgs(TRUE)

infile = args[1]
outfile = args[2]
threshold = args[3]
logpath = args[4]

packages = c("reshape2", "readr", "logr")

for(package in packages)
  library(package, character.only = T)

## removing scientific notation
options(scipen = 100)

########################
## DEFINING FUNCTIONS
########################

rep.filter <- function(df, threshold){

  ############# FUNC DESCRIPTION
  ## func only selects copies whose length is above the threshold
  ## if threshold == "none" -> all copies retained
  ## if threshold == "median" -> copies above family median length retained
  ############# FUNC DESCRIPTION
  
  if(threshold == "median") {
    
    df$width = df$end-df$start ## calculating TE width !!! RANGES MUST BE 0 BASED !!!
    
    ## computing family medians
    medians = tapply(df$width, df$repName, FUN=median)
    
    ## selecting TE copies above family median length
    res = do.call(rbind,
                  lapply(split(df, f=df$repName), function(x, medians) {
                    return(x[x$width >= medians[names(medians) == unique(x$repName)], ])
                    }, medians))
    row.names(res) = NULL
    
    return(res[, -c(5)])
  } else {
    return(df)
  }
}


###############
## READING DATA
###############

## reading repeatmasker
## debug -> rmsk <- readr::read_tsv(file = "/media/yari-msi/DATA1/rmsk_preprocessed.bed", col_names=F)
rmsk <- readr::read_tsv(file = infile, col_names=F) ## !!! REPEATMASKER RANGES ARE 0-BASED !!!
rmsk = as.data.frame(rmsk)
colnames(rmsk) = c("chrom", "start", "end", "repName")

#######
## MAIN
#######

rmsk_filtered = rep.filter(rmsk, threshold)

frac_lost_all = round(nrow(rmsk_filtered)/nrow(rmsk), 2)

frac_lost = unlist(Map(function(x, y) nrow(x)/nrow(y), split(rmsk_filtered, f=rmsk_filtered$repName), split(rmsk, f=rmsk$repName)))

## output (0-based)
readr::write_tsv(rmsk_filtered, file=outfile, col_names = F, quote = "none")

## log
tmp = file.path(logpath)
logfile = log_open(tmp)
log_print("fraction of ALL TEs retained:\n")
log_print(frac_lost_all)
log_print("fraction of TEs retained by repName:\n")
log_print(frac_lost)
log_close()

##

