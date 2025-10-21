## script to collect all data (predictors and outcomes) needed for linear regression models
args = commandArgs(TRUE)
infile_predictors = args[1]
infile_gc <- args[2]
infile_h2 <- args[3]
infile_egenes <- args[4]
infile_te_list <- args[5]
outfile <- args[6]

library(readxl)
library(plyr)

min_coverage <- 10000

data <- read.delim(infile_predictors)
load(infile_gc)
data <- merge(data, gc, by.x = "repName", by.y = 0, all = TRUE)

## outcomes
herit <- read_xlsx(infile_h2)
herit$TE_Name <- gsub("_", "-", herit$TE_Name)
herit <- herit[, c("TE_Name", "Enr")]
# remove ambiguous (repName or repFamily)
amb <- c("L2", "MIR")
herit <- herit[!herit$TE_Name %in% amb, ]

data <- merge(data, herit, by.x = "repName", by.y = "TE_Name", all.x = TRUE)

load(infile_egenes)
data <- merge(data, egene[, c("CLASS", "egenes_per_snp_common")], by.x = "repName", by.y = "CLASS", all.x = TRUE)

te_list <- scan(infile_te_list, what = character(0))
data <- data[data$repName %in% te_list, ]
data$DRS <- 1 - data$DRS

rownames(data) <- data$repName
data$repName <- NULL

list = c("data")
save(list = list, file = outfile)
