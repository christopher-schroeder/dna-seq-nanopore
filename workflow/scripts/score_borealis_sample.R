suppressPackageStartupMessages({
    library(data.table)
    library(gamlss.dist)
})

wd         <- getwd()
cov_file   <- file.path(wd, snakemake@input[["cov"]])
model_file <- file.path(wd, snakemake@input[["model"]])
out_file   <- file.path(wd, snakemake@output[["outliers"]])

chrs        <- c(as.character(seq_len(22)), "X", "Y")
minObsDepth <- 10L

model <- fread(model_file)
setkeyv(model, c("chr", "pos"))

# Bismark cov.gz: chr, start, end, pct_meth, count_meth, count_unmeth
dt <- fread(cov_file, select = c(1L, 2L, 5L, 6L),
            col.names = c("chr", "pos", "meth", "unmeth"))
dt <- dt[chr %in% chrs]
dt[, x := as.integer(meth)]
dt[, n := as.integer(meth + unmeth)]
dt[, c("meth", "unmeth") := NULL]
dt <- dt[n >= minObsDepth]
setkeyv(dt, c("chr", "pos"))

# Inner join with model
dt <- model[dt, nomatch = 0]
setcolorder(dt, c("chr", "pos", "x", "n", "mu", "theta"))

# Compute p-values (replicates borealis:::computePvalsAndEffSize)
dt[, ltp := pmax(0, pBB(x, mu = mu, sigma = theta, bd = n))]
dt[, rtp := pmax(0, pBB(pmax(x - 1L, 0L), mu = mu, sigma = theta,
                         bd = n, lower.tail = FALSE))]
dt[, pVal    := pmin(1, 2 * pmin(rtp, ltp))]
dt[, isHypo  := NA]
dt[, c("ltp", "rtp") := NULL]
dt[, effSize := (x / n) - mu]
dt[effSize < -0.1, isHypo := TRUE]
dt[effSize >  0.1, isHypo := FALSE]

setorder(dt, pVal)
# Output format matches original borealis writeResults: chr pos pos x n mu theta pVal isHypo effSize
out <- as.data.frame(dt)[, c(1, 2, 2, 3:ncol(dt))]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.table(out, file = out_file, quote = FALSE, row.names = FALSE, sep = "\t")
message("Written: ", out_file)
