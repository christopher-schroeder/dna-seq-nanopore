suppressPackageStartupMessages({
    library(borealis)
    library(foreach)
    library(doParallel)
    library(data.table)
})

wd        <- getwd()
cov_dir   <- file.path(wd, snakemake@params[["cov_dir"]])
out_dir   <- file.path(wd, snakemake@params[["out_dir"]])
n_threads <- snakemake@threads
done_file <- file.path(wd, snakemake@output[["done"]])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chrs       <- c(as.character(seq_len(22)), "X", "Y")
chr_suffix <- paste0("_", paste(chrs, collapse = "_"), "_DMLs.tsv")
minObsDepth <- 10L

cov_files    <- list.files(cov_dir, pattern = "\\.bismark\\.cov\\.gz$", full.names = TRUE)
sample_names <- sub("\\.bismark\\.cov\\.gz$", "", basename(cov_files))
unique_idx   <- !duplicated(sample_names)
cov_files    <- cov_files[unique_idx]
sample_names <- sample_names[unique_idx]
n_samp       <- length(cov_files)
message("Found ", n_samp, " unique samples.")

# Phase 1: build M and Cov matrices without bsseq (avoids ~150 GB peak of makeBSseqData).
# loadBismarkData collects all 130 data.frames before calling makeBSseqData, which
# requires holding ~73 GB of input + ~30 GB output simultaneously.
# Here we: (a) collect the CpG position universe, then (b) fill matrices one file at a time.

message("Phase 1a: collecting CpG position universe from ", n_samp, " samples...")
# streaming anti-join avoids the 2^31 row limit of rbindlist on 130*28M rows
pos_set <- fread(cov_files[1], select = c(1L, 2L),
                 col.names = c("chr", "pos"))[chr %in% chrs]
setkeyv(pos_set, c("chr", "pos"))
for (i in seq_along(cov_files)[-1]) {
    tmp <- fread(cov_files[i], select = c(1L, 2L),
                 col.names = c("chr", "pos"))[chr %in% chrs]
    setkeyv(tmp, c("chr", "pos"))
    new_pos <- tmp[!pos_set]
    if (nrow(new_pos) > 0) {
        pos_set <- rbindlist(list(pos_set, new_pos))
        setkeyv(pos_set, c("chr", "pos"))
    }
    rm(tmp, new_pos)
}
gc()

chr_order <- c(as.character(seq_len(22)), "X", "Y")
pos_set[, chr_int := match(chr, chr_order)]
setorder(pos_set, chr_int, pos)
pos_set[, chr_int := NULL]
pos_set[, row_idx := .I]
setkeyv(pos_set, c("chr", "pos"))
npos <- nrow(pos_set)
message("Phase 1a complete: ", npos, " CpG sites in union.")

message("Phase 1b: loading methylation values into matrices...")
M   <- matrix(0L, nrow = npos, ncol = n_samp, dimnames = list(NULL, sample_names))
Cov <- matrix(0L, nrow = npos, ncol = n_samp, dimnames = list(NULL, sample_names))

for (i in seq_along(cov_files)) {
    # Bismark cov.gz: chr, start, end, pct_meth, count_meth, count_unmeth
    tmp <- fread(cov_files[i], select = c(1L, 2L, 5L, 6L),
                 col.names = c("chr", "pos", "meth", "unmeth"))
    tmp <- tmp[chr %in% chrs]
    tmp[, N := meth + unmeth]
    setkeyv(tmp, c("chr", "pos"))
    m <- pos_set[tmp, nomatch = 0]
    M[m$row_idx, i]   <- as.integer(m$meth)
    Cov[m$row_idx, i] <- as.integer(m$N)
    rm(tmp, m)
    if (i %% 20 == 0) { gc(); message("  loaded ", i, "/", n_samp) }
}
gc()
message("Phase 1b complete.")

# Phase 2: build beta-binomial models with SOCK workers + pre-split chunks.
chr_v <- pos_set$chr
pos_v <- pos_set$pos

keepInd <- rowSums(Cov >= 4L) >= 5L
x_f   <- M[keepInd, , drop = FALSE]
n_f   <- Cov[keepInd, , drop = FALSE]
chr_f <- chr_v[keepInd]
pos_f <- pos_v[keepInd]

niter <- nrow(x_f)
message("Looping over ", niter, " positions in model building. It will take some time.")

chunks <- lapply(seq_len(n_threads), function(i) {
    idx <- seq(i, niter, by = n_threads)
    list(chr = chr_f[idx], pos = pos_f[idx],
         x   = x_f[idx, , drop = FALSE],
         n   = n_f[idx, , drop = FALSE])
})
rm(x_f, n_f, chr_f, pos_f); gc()
message("Chunks built, starting SOCK cluster with ", n_threads, " workers.")

fitGamlss_fn <- borealis:::fitGamlss

workers <- parallel::makeCluster(n_threads, type = "SOCK")
registerDoParallel(workers)

modelDF <- foreach(
    chunk    = chunks,
    .combine = rbind,
    .export  = "fitGamlss_fn"
) %dopar% {
    nc    <- length(chunk$chr)
    dfInd <- data.frame(chr = chunk$chr, pos = chunk$pos,
                        mu = NA_real_, theta = NA_real_)
    for (i in seq_len(nc)) {
        pred <- tryCatch(
            fitGamlss_fn(chunk$x[i, ], chunk$n[i, ], 4L, TRUE, 10),
            error = function(e) NULL
        )
        if (!is.null(pred)) {
            dfInd$mu[i]    <- pred$mu[1]
            dfInd$theta[i] <- pred$sigma[1]
        }
    }
    dfInd
}

stopCluster(workers)
rm(chunks); gc()
message("Modeling complete. Writing model CSV.")

modelDF <- na.omit(modelDF)
modelDF <- modelDF[order(modelDF$chr, modelDF$pos), ]

write.csv(
    modelDF,
    file      = file.path(out_dir, paste0("CpG_model_", paste(chrs, collapse = "_"), ".csv")),
    row.names = FALSE, quote = FALSE
)

# Phase 3: write per-sample results.
# Replicates borealis:::writeResults without needing a BSseq object.
# Output format (identical to original): chr pos pos x n mu theta pVal isHypo effSize
# (the duplicate pos column is borealis's own convention).

message("Phase 3: writing per-sample results...")
setDT(modelDF); setkeyv(modelDF, c("chr", "pos"))
chr_dt <- as.character(chr_v)
pos_dt <- as.integer(pos_v)

for (i in seq_along(cov_files)) {
    samp <- sample_names[i]
    dt <- data.table(
        chr = chr_dt,
        pos = pos_dt,
        x   = as.integer(M[, i]),
        n   = as.integer(Cov[, i])
    )
    dt <- dt[n >= minObsDepth]
    setkeyv(dt, c("chr", "pos"))
    # inner join: modelDF[dt] returns modelDF cols then dt's non-key cols
    dt <- modelDF[dt, nomatch = 0]
    # column order after join: chr, pos, mu, theta, x, n
    # reorder to match original writeResults: chr(1), pos(2), x(3), n(4), mu(5), theta(6)
    setcolorder(dt, c("chr", "pos", "x", "n", "mu", "theta"))

    # compute p-values (replicates borealis:::computePvalsAndEffSize)
    dt[, leftTailProb  := pmax(0, gamlss.dist::pBB(x, mu = mu, sigma = theta, bd = n))]
    dt[, rightTailProb := pmax(0, gamlss.dist::pBB(pmax(x - 1L, 0L), mu = mu, sigma = theta,
                                                    bd = n, lower.tail = FALSE))]
    dt[, pVal    := pmin(1, 2 * pmin(rightTailProb, leftTailProb))]
    dt[, isHypo  := NA]
    dt[, leftTailProb  := NULL]
    dt[, rightTailProb := NULL]
    dt[, effSize := (x / n) - mu]
    dt[effSize < -0.1, isHypo := TRUE]
    dt[effSize >  0.1, isHypo := FALSE]

    setorder(dt, pVal)
    # convert to data.frame for duplicate-column indexing (original format has pos twice)
    out <- as.data.frame(dt)[, c(1, 2, 2, 3:ncol(dt))]

    write.table(
        out,
        file      = file.path(out_dir, paste0(samp, chr_suffix)),
        quote     = FALSE,
        row.names = FALSE,
        sep       = "\t"
    )
}
rm(M, Cov); gc()
message("Phase 3 complete.")

for (samp in sample_names) {
    src <- file.path(out_dir, paste0(samp, chr_suffix))
    dst <- file.path(out_dir, paste0(samp, ".tsv"))
    if (file.exists(src)) file.rename(src, dst)
}

file.create(done_file)
