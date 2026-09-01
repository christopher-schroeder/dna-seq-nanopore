suppressPackageStartupMessages({
    library(borealis)
    library(foreach)
    library(doParallel)
    library(data.table)
})

wd        <- getwd()
cov_dir   <- file.path(wd, snakemake@params[["cov_dir"]])
out_csv   <- file.path(wd, snakemake@output[["model"]])
n_threads <- snakemake@threads

chrs <- c(as.character(seq_len(22)), "X", "Y")

cov_files    <- list.files(cov_dir, pattern = "\\.bismark\\.cov\\.gz$", full.names = TRUE)
sample_names <- sub("\\.bismark\\.cov\\.gz$", "", basename(cov_files))
unique_idx   <- !duplicated(sample_names)
cov_files    <- cov_files[unique_idx]
sample_names <- sample_names[unique_idx]
n_samp       <- length(cov_files)
message("Found ", n_samp, " unique samples.")

# Spawn SOCK workers NOW while the parent process is tiny (~200 MB).
# Loading M/Cov later grows the parent but never forks again — no OOM risk.
message("Starting SOCK cluster with ", n_threads, " workers (parent still tiny)...")
workers <- parallel::makeCluster(n_threads, type = "SOCK")
clusterEvalQ(workers, suppressPackageStartupMessages(library(borealis)))
registerDoParallel(workers)
fitGamlss_fn <- borealis:::fitGamlss

# Phase 1: load M and Cov matrices (one file at a time, no bsseq).
message("Phase 1a: collecting CpG position universe...")
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
message("Phase 1a complete: ", npos, " CpG sites.")

message("Phase 1b: loading methylation matrices...")
M   <- matrix(0L, nrow = npos, ncol = n_samp, dimnames = list(NULL, sample_names))
Cov <- matrix(0L, nrow = npos, ncol = n_samp, dimnames = list(NULL, sample_names))
for (i in seq_along(cov_files)) {
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

chr_v <- pos_set$chr
pos_v <- pos_set$pos

# Phase 2: model fitting chromosome by chromosome.
# Each chr slice is ~1–3 GB, so chunks sent to workers are small.
message("Phase 2: fitting models per chromosome...")
all_models <- vector("list", length(chrs))

for (ci in seq_along(chrs)) {
    chr <- chrs[ci]
    chr_idx <- which(chr_v == chr)
    if (length(chr_idx) == 0L) next

    M_chr   <- M[chr_idx,   , drop = FALSE]
    Cov_chr <- Cov[chr_idx, , drop = FALSE]
    pos_chr <- pos_v[chr_idx]

    keepInd <- rowSums(Cov_chr >= 4L) >= 5L
    x_f <- M_chr[keepInd,   , drop = FALSE]
    n_f <- Cov_chr[keepInd, , drop = FALSE]
    pf  <- pos_chr[keepInd]
    rm(M_chr, Cov_chr); gc()

    niter <- nrow(x_f)
    message("  chr", chr, ": ", niter, " positions")
    if (niter == 0L) {
        all_models[[ci]] <- data.frame(chr = character(0), pos = integer(0),
                                       mu  = numeric(0),   theta = numeric(0))
        next
    }

    nt <- min(n_threads, niter)
    chunks <- lapply(seq_len(nt), function(i) {
        idx <- seq(i, niter, by = nt)
        list(chr = rep(chr, length(idx)), pos = pf[idx],
             x   = x_f[idx, , drop = FALSE],
             n   = n_f[idx, , drop = FALSE])
    })
    rm(x_f, n_f, pf); gc()

    chr_model <- foreach(
        chunk    = chunks,
        .combine = rbind,
        .export  = "fitGamlss_fn"
    ) %dopar% {
        nc    <- length(chunk$chr)
        dfInd <- data.frame(chr = chunk$chr, pos = chunk$pos,
                            mu  = NA_real_,  theta = NA_real_)
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
    rm(chunks); gc()
    all_models[[ci]] <- na.omit(chr_model)
}

stopCluster(workers)
rm(M, Cov); gc()

modelDF <- rbindlist(all_models)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
fwrite(modelDF, file = out_csv)
message("Model written to ", out_csv, " (", nrow(modelDF), " positions)")
