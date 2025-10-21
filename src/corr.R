## script to collect all data (predictors and outcomes) needed for linear regression models
args = commandArgs(TRUE)
infile_data = args[1]
outfile <- args[2]

load(infile_data)
data <- data[, c("median_age", "base_ratio", "sequence_score", "DRS", "phyloP", "gc_content", "Enr", "egenes_per_snp_common")]
corr <- cor(data, use = "pairwise", method = "spearman")

list = c("corr")
save(list = list, file = outfile)
