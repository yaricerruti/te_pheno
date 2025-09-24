#!/usr/bin/Rscript

args = commandArgs(T)
infile = args[1]
outfile = args[2]

## removing scientific notation
options(scipen=999)

packages = c("readr")

for(package in packages)
  library(package, character.only = T)

################
## READING DATA
###############

## debug -> phylo <- readr::read_tsv(file = "/mnt/filippo/prj/te_enhancer_conservation/predictors_dataset/v1/phyloP_dELS_concatenated.tsv", col_names = F)
phylo <- readr::read_tsv(file = infile, col_names = F)
phylo = as.data.frame(phylo)

##############
## MAIN
##############

colnames(phylo) = c("TAG", "width", "width.covered", "phyloP.sum", "phyloP.mean0", "phyloP.mean", "GROUP")
phylo = phylo[, c("GROUP", "TAG", "width", "width.covered", "phyloP.sum", "phyloP.mean")]

## computing coverage and adding it as column, without removing TEs
# phylo$coverage = phylo$width.covered/phylo$width ## NOT NECESSARY NOW, CLASS-WISE COVERAGE WILL BE COMPUTED LATER

## not filtering TEs
# phylo = data.table(phylo)
# phylo = phylo[phylo$base.coverage.ratio > 0.9]
# phylo = as.data.frame(phylo)
# phylo$base.coverage.ratio = NULL
# phylo$width = NULL

## mean and median of all elements belonging to a class, OVL or NONOVL, in each pseudochromosome
## splitting data in a list of lists:
##               - CLASS1@OVL@1_1
##       -real -|
##      |        - CLASS2@OVL@1_1 ...
## data-|
##      |        - CLASS1@OVL@1_1
##       -rand -|
##               - CLASS2@OVL@1_1 ...
##

phylo = split(phylo, f=phylo$GROUP)
phylo_list = lapply(phylo, function(x) {
  l = split(x, f=x$TAG)
  return(l)})

## bigWigAverageOverBed outputs a dataframe in the form of phylo (see above).
## it contains the nucleotide-based phyloP values averaged over each instance of
## a bed file containing annotated OVL vs NONOVL transposable elements.
## the phylo column "phyloP.sum" contains the sum of phyloP values of each base
## contained within each TE (each line of the file).
##the phylo column "width.covered" contains the number of bases covered by the
## phyloP count, in each TE, and it (usually) corresponds to the length of the
## TE itself. We are averaging a "class-based" phyloP sum across all TEs
## belonging to a class. instead of simply computing the mean of the TE-based
## phyloP, which would effectively be a "mean of means" (nt-based -> TE based -> class-based)
## we sum the TE-based phyloP.sum values of all TEs belonging to a specific class
## (which already are the sum of the phyloP value of all bases belonging to a TE
## instance), then we sum the total length of TEs belonging to a class, then divide
## the class-based phyloP sum by the class-based total length, effectively obtaining
## the average phyloP for each class instance.
## our instances are "TE-class@overlap-status@pseudochromosome" objects.

phylo_mean <- data.frame()
for(n in names(phylo_list)){
  temp_df <- lapply(phylo_list[[n]], function(x){
    df = data.frame("TAG" = unique(x$TAG),
                    "n.elements" = nrow(x),
                    "n.bases" = sum(x$width),
                    "n.covered.bases" = sum(x$width.covered),
                    "phyloP.class.sum" = sum(x$phyloP.sum)) ## class-based phyloP sum over all TE instances of a TE-class@overlap-status@pseudochromosome obj
    lillie = tryCatch(expr = {l = nortest::lillie.test(x$phyloP.mean)$p.value}, ## variant of Kolmogorov-Smirnov test to check whether the phyloP mean values are normally distributed
                      error = function(e) {print(paste0("error when testing ", unique(x$TAG), " : --- ", conditionMessage(e))) ## sample size must be > 4
                        return(NULL)})
    if(is.null(lillie))
      df$Lilliefors.pvalue = NA
    else
      df$Lilliefors.pvalue = round(lillie, 8)
    return(df)
  })
  temp_df = do.call(rbind, temp_df)
  row.names(temp_df) = NULL
  temp_df$Lilliefors.bonf = p.adjust(temp_df$Lilliefors.pvalue, method = "bonferroni")
  temp_df$GROUP = n
  temp_df = temp_df[, c("GROUP", "TAG", "n.elements", "n.bases", "n.covered.bases", "phyloP.class.sum", "Lilliefors.pvalue", "Lilliefors.bonf")]
  phylo_mean = rbind(phylo_mean, temp_df)
}
rm(temp_df)

## computing the average class-based phyloP
phylo_mean$phyloP.class.average = phylo_mean$phyloP.class.sum/phylo_mean$n.covered.bases

# computing the class-wise coverage (count of all bases covered in a class / all bases of a class)
phylo_mean$class.coverage = phylo_mean$n.covered.bases/phylo_mean$n.bases

newcols <- c("CLASS", "TYPE")
i=1
for(n in newcols){
  phylo_mean[, n] = sapply(strsplit(phylo_mean[, "TAG"], "@"), getElement, i)
  i=i+1
}
phylo_mean$TAG = NULL

phylo_mean = phylo_mean[, c("GROUP", "CLASS", "TYPE", "n.elements", "n.bases", "n.covered.bases", "class.coverage", "phyloP.class.sum", "phyloP.class.average", "Lilliefors.pvalue", "Lilliefors.bonf")]

readr::write_tsv(phylo_mean, file = outfile, col_names = T, quote = "none")

##
