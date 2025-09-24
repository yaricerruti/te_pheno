#!/usr/bin/Rscript

args = commandArgs(T)
infile = args[2:length(args)]
outfile = args[1]

## removing scientific notation
options(scipen=999)

packages = c("readr")

for(package in packages)
  library(package, character.only = T)

################
## READING DATA
###############

## debug -> infile = c("DRS/dELS_real_te_func_ovl_weighted_DRS_pseudochroms.bed", "DRS/dELS_rand1_te_func_ovl_weighted_DRS_pseudochroms.bed")

DRS = data.frame()
for(f in infile){
  temp = readr::read_tsv(file = f, col_names = F)
  temp = as.data.frame(temp)
  temp$GROUP = gsub("[0-9]", "", unlist(strsplit(f, "_"))[2])
  DRS = rbind(DRS, temp)
}
colnames(DRS) = c("chr", "start", "end", "TAG", "DRS", "GROUP")

##############
## MAIN
##############

## computing width
DRS$width = DRS$end-DRS$start

## removing coordinates
DRS = DRS[, c("TAG", "DRS", "GROUP", "width")]

DRS = split(DRS, f=DRS$GROUP)
DRS_list = lapply(DRS, function(x) {
  l = split(x, f=x$TAG)
  return(l)})

#######################
## computing DRS median

## since a DR score has been assigned to each TE, family-wise median DRS is computed
DRS_mean <- data.frame()
for(n in names(DRS_list)){
  temp_df <- lapply(DRS_list[[n]], function(x){
    df = data.frame("TAG" = unique(x$TAG),
                    "element_count" = nrow(x),
                    "coverage" = sum(x$width),
                    "median_DRS" = median(x$DRS))
    lillie = tryCatch(expr = {l = nortest::lillie.test(x$DRS)$p.value},
                      error = function(e) {print(paste0("error when testing ", unique(x$TAG), " : --- ", conditionMessage(e))) ## sample size must be > 4
                        return(NULL)})
    if(is.null(lillie))
      df$Lillie_p = NA
    else
      df$Lillie_p = round(lillie, 4)
    return(df)
  })
  temp_df = do.call(rbind, temp_df)
  row.names(temp_df) = NULL
  temp_df$Lillie_bonf = p.adjust(temp_df$Lillie_p, method = "bonferroni")
  temp_df$GROUP = n
  temp_df = temp_df[, c("GROUP", "TAG", "element_count", "coverage", "median_DRS", "Lillie_p", "Lillie_bonf")]
  DRS_mean = rbind(DRS_mean, temp_df)
}
rm(temp_df)

newcols <- c("CLASS", "TYPE")
i=1
for(n in newcols){
  DRS_mean[, n] = sapply(strsplit(DRS_mean[, "TAG"], "@"), getElement, i)
  i=i+1
}
DRS_mean$TAG = NULL

DRS_mean = DRS_mean[, c("GROUP", "CLASS", "TYPE", "element_count", "coverage", "median_DRS", "Lillie_p", "Lillie_bonf")]

readr::write_tsv(DRS_mean, file = outfile, col_names = T, quote = "none")

##
