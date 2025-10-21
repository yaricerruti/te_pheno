## linear regression
args = commandArgs(TRUE)
infile_data = args[1]
outfile <- args[2]

library(car)
library(plyr)

my_merge <- function(df1, df2, ...) {
    df <- merge(df1, df2, by = 0, ...)
    rownames(df) <- df[,1]
    df[,1] <- NULL
    df
}

do_invn <- function(x) {
    qnorm((rank(x,na.last="keep", ties.method = "random")-0.5)/sum(!is.na(x)))
}
do_all_invn <- function(data) {
    data[, -1] <- apply(data[, -1], 2, do_invn)
    data
}
do_uni <- function(pr, pr_name, out) {
    pr <- pr[, pr_name, drop = FALSE]
    df <- na.omit(my_merge(pr, out))
    df[, 1] <- scale(df[, 1])
    df[, 2] <- do_invn(df[, 2])
    fm <- as.formula(paste(colnames(out), "~", pr_name))
    md <- lm(formula = fm, data = df)
}
do_all <- function(pr, out) {
    ## various regressions of one outcome on several predictors
    ## pr: a dataframe of predictors
    ## out: a dataframe with a single column of outcomes
    ## matching is done by rownames
    
    ## univariable regressions
    prl <- colnames(pr)
    names(prl) <- prl
    uni <- lapply(prl, do_uni, pr = pr, out = out)
    uni_sum <- lapply(uni, summary)
    uni_coeff <- ldply(lapply(uni_sum, function(x) x$coefficients[-1, ]), .id = "predictor")
    rownames(uni_coeff) <- uni_coeff$predictor
    uni_coeff$predictor <- NULL
    uni_coeff$r2 <- sapply(uni_sum, `[[`, "r.squared")
    
    ## multivariable regression without interactions
    df <- na.omit(my_merge(pr, out))
    df[, colnames(pr)] <- scale(df[, colnames(pr)])
    df[, colnames(out)] <- do_invn(df[, colnames(out)])
    fm <- as.formula(paste(colnames(out), "~ ."))
    multi <- lm(fm, data = df)
    multi_sum <- summary(multi)
    multi_coeff <- multi_sum$coefficients[-1, ]
    multi_r2 <- multi_sum$r.squared
    multi_adjr2 <- multi_sum$adj.r.squared
    multi_vif <- vif(multi)
    
    ## multivariable regression with all pairwise interactions
    ints <- apply(combn(prl, 2), 2 , paste, collapse = "*")
    fmr <- paste(ints, collapse = " + ")
    fml <- colnames(out)
    fm <- as.formula(paste(fml, fmr, sep = "~"))
    ints <- lm(fm, data = df)
    ints_sum <- summary(ints)
    ints_coeff <- ints_sum$coefficients[-1, ]
    ints_r2 <- ints_sum$r.squared
    ints_adjr2 <- ints_sum$adj.r.squared
    ints_vif <- vif(ints, type = "predictor")
    
    list(uni = uni, uni_sum = uni_sum, uni_coeff = uni_coeff, multi = multi, multi_coeff = multi_coeff, multi_r2 = multi_r2, 
         multi_adjr2 = multi_adjr2, multi_vif = multi_vif, ints = ints, ints_coeff = ints_coeff, ints_r2 = ints_r2, 
         ints_adjr2 = ints_adjr2, ints_vif = ints_vif)
}

load(infile_data)
data <- data[, c("median_age", "base_ratio", "sequence_score", "DRS", "phyloP", "gc_content", 
                 "Enr", "egenes_per_snp_common")]

outl <- c("Enr", "egenes_per_snp_common", "base_ratio")
names(outl) <- outl

pred <- data[, c("median_age", "sequence_score", "DRS", "phyloP", "gc_content")]
outcome <- data[, outl]
do_model <- function(out) {
    outdf <- outcome[, out, drop = FALSE]
    do_all(pr = pred, out = outdf)
}

model <- lapply(outl, do_model)

list = c("model")
save(list = list, file = outfile)

