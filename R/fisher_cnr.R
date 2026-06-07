## Fisher distance and simulation
## Internalized from SCclust (Krasnitz et al., MIT License)
## https://github.com/KrasnitzLab/SCclust

#' Symmetric cell-by-cell matrix from upper-triangle vector
#' @param data numeric vector of upper-triangle values (length = n*(n-1)/2)
#' @param cell_names character vector of length n
#' @return symmetric n x n matrix with dimnames set to cell_names
#' @keywords internal
#' @noRd
.cell2cell_matrix <- function(data, cell_names) {
    n <- length(cell_names)
    m <- matrix(0, n, n, dimnames = list(cell_names, cell_names))
    m[upper.tri(m)] <- data
    pmin(m, t(m))
}

#' Metropolis-Hastings sampler for permuting binary pin vectors
#'
#' Performs \code{sweeps} Metropolis steps that swap 1s and 0s while
#' preserving the marginal pin probabilities \code{p}.  Used to generate
#' null-distributed pin matrices for the Fisher simulation.
#'
#' @param x binary vector (0/1) representing one cell's pin occupancy
#' @param p numeric vector of marginal probabilities (same length as x)
#' @param sweeps number of proposed swaps
#' @return permuted binary vector with same length as x
#' @keywords internal
#' @noRd
.metro <- function(x, p, sweeps) {
    for (i in seq_len(sweeps)) {
        occ <- which(x == 1 & p < 1); vac <- which(x == 0 & p > 0)
        if (!length(occ) || !length(vac)) return(x)
        if (length(occ) > length(vac)) { orig <- sample(occ, length(vac)); dest <- vac }
        else                           { orig <- occ; dest <- sample(vac, length(occ)) }
        mv <- p[dest] * (1 - p[orig]) / (p[orig] * (1 - p[dest])) > stats::runif(length(orig))
        x[dest[mv]] <- 1; x[orig[mv]] <- 0
    }
    x
}

#' Fisher exact test p-values for one row against all other rows
#'
#' Computes one-sided Fisher exact test p-values between cell \code{i} and
#' all cells \code{j >= i} from pre-computed contingency table components.
#'
#' @param i row index (cell) to test against all others
#' @param yy matrix of both-positive co-occurrences (cells x cells)
#' @param ny matrix of row-positive, col-negative co-occurrences
#' @param nn matrix of both-negative co-occurrences
#' @return numeric vector of p-values (length = nrow(yy))
#' @keywords internal
#' @noRd
.fast_fisher <- function(i, yy, ny, nn) {
    ftp <- rep(1, nrow(yy))
    for (j in i:nrow(yy))
        ftp[j] <- stats::fisher.test(
            matrix(c(yy[i,j], ny[i,j], ny[j,i], nn[i,j]), 2, 2),
            alternative = "greater")$p.value
    ftp[ftp == 0] <- 2 * .Machine$double.xmin
    ftp
}

#' Monte Carlo Fisher simulation for pairwise cell similarity
#'
#' Runs \code{nsim} permutation simulations of the Fisher combined test across
#' pin sign groups.  Each simulation permutes pin occupancy per cell using
#' Metropolis-Hastings (\code{.metro}) then computes pairwise Fisher p-values.
#' Returns a matrix of simulated p-values used as the null distribution by
#' \code{.fisher_fdr}.
#'
#' @param m list of binary matrices (one per sign group), each pins x cells
#' @param nsim number of simulation iterations
#' @param nsweep number of Metropolis sweeps per permutation step
#' @param seedme integer random seed
#' @param njobs number of parallel workers (passed to \code{parallel::makeCluster})
#' @param combo combination method: \code{"fisher"} (default) or \code{"stouffer"}
#' @return numeric matrix (cell-pairs x nsim) of simulated combined p-values
#' @keywords internal
#' @noRd
.sim_fisher <- function(m, nsim, nsweep, seedme, njobs = 1, combo = "fisher") {
    ncores <- min(njobs, parallel::detectCores())
    cl <- parallel::makeCluster(ncores, setup_strategy = "sequential")
    on.exit(parallel::stopCluster(cl))
    parallel::clusterSetRNGStream(cl)
    RNGkind("L'Ecuyer-CMRG"); set.seed(seedme)
    rf   <- lapply(m, rowMeans)
    npairs <- ncol(m[[1]]) * (ncol(m[[1]]) - 1) / 2
    tp   <- matrix(ncol = nsim, nrow = npairs)
    xmat <- matrix(ncol = length(m), nrow = npairs)
    for (i in seq_len(nsim)) {
        if (length(m) <= 1) next
        for (j in seq_along(m)) {
            if (nsweep > 0 && length(rf[[j]]) > 1)
                m[[j]] <- parallel::parApply(cl, m[[j]], 2, .metro,
                                             p = rf[[j]], sweeps = nsweep)
            yy <- t(m[[j]]) %*% m[[j]]
            ny <- t(1 - m[[j]]) %*% m[[j]]
            nn <- t(1 - m[[j]]) %*% (1 - m[[j]])
            lbi <- seq_len(nrow(yy))
            lbi[lbi %% 2 == 0] <- nrow(yy) + 2 - nrow(yy) %% 2 - lbi[lbi %% 2 == 0]
            x <- parallel::parSapply(cl, lbi, .fast_fisher,
                                     yy = yy, ny = ny, nn = nn)[, order(lbi)]
            x <- pmin(x, t(x))
            xmat[, j] <- x[upper.tri(x)]
        }
        if (length(m) == 1) {
            tp[, i] <- xmat[, 1]
        } else if (combo == "fisher") {
            Xsq      <- -2 * rowSums(log(xmat))
            tp[, i]  <- sapply(Xsq, stats::pchisq, df = 2 * ncol(xmat),
                                lower.tail = FALSE)
        } else {
            Z       <- colSums(apply(1 - xmat, 1, stats::qnorm)) / sqrt(ncol(xmat))
            tp[, i] <- 1 - sapply(Z, stats::pnorm)
        }
    }
    tp[tp == 0] <- 2 * .Machine$double.xmin
    tp
}

#' Wrapper: run true and simulated Fisher tests across pin sign groups
#'
#' Splits the pin matrix by breakpoint sign group, runs one real and
#' \code{nsim} permuted Fisher tests, and returns a list with the true and
#' simulated p-value vectors for use by \code{.fisher_fdr}.
#'
#' @param pinmat_df data.frame of pin matrix (pins x cells, binary)
#' @param pins_df data.frame with columns \code{bin} and \code{sign}
#' @param njobs number of parallel workers (NULL = detectCores() - 4)
#' @param nsim number of simulation iterations (default 150)
#' @param nsweep number of Metropolis sweeps per permutation step (default 200)
#' @param seedme integer random seed (default 123)
#' @return list with \code{true} and \code{sim} p-value matrices, or NULL if
#'   fewer than 2 sign groups are present
#' @keywords internal
#' @noRd
.sim_fisher_wrapper <- function(pinmat_df, pins_df, njobs = NULL,
                                nsim = 150, nsweep = 200, seedme = 123) {
    if (is.null(njobs)) njobs <- max(1L, parallel::detectCores() - 4L)
    signs <- unique(pins_df[, "sign"])
    if (length(signs) <= 1) return(NULL)
    m     <- lapply(signs, function(s)
                    as.matrix(pinmat_df[pins_df[, "sign"] == s, , drop = FALSE]))
    vtrue <- .sim_fisher(m, nsim = 1,    nsweep = 0,      seedme = seedme, njobs = njobs)
    msim  <- .sim_fisher(m, nsim = nsim, nsweep = nsweep, seedme = seedme, njobs = njobs)
    list(true = vtrue, sim = msim)
}

#' Compute pairwise log10(FDR) matrix from Fisher p-values and their null
#'
#' Estimates FDR for each observed pairwise p-value by comparing its
#' empirical CDF to the simulated null CDF, using piecewise linear
#' interpolation in log space.
#'
#' @param true_pv numeric vector of observed pairwise p-values (upper triangle)
#' @param sim_pv numeric vector of all simulated null p-values
#' @param cell_names character vector of cell IDs (length n, where n*(n-1)/2 == length(true_pv))
#' @param lmmax lower CDF quantile below which a linear model is used for extrapolation
#' @return symmetric n x n matrix of log10(FDR) values (≤ 0; more negative = more significant)
#' @keywords internal
#' @noRd
.fisher_fdr <- function(true_pv, sim_pv, cell_names, lmmax = 0.001) {
    assertthat::assert_that(length(cell_names) == (1 + sqrt(1 + 8 * length(true_pv))) / 2)
    sim_sort   <- sort(sim_pv);  sim_unique  <- unique(sim_sort)
    true_sort  <- sort(true_pv); true_unique <- unique(true_sort)
    sim_count  <- tapply(match(sim_sort,  sim_unique),  match(sim_sort,  sim_unique),  length)
    true_count <- tapply(match(true_sort, true_unique), match(true_sort, true_unique), length)
    z <- cbind(c(log(true_unique), log(sim_unique)),
               c(log(cumsum(true_count) / sum(true_count)),
                 log(cumsum(sim_count)  / sum(sim_count))),
               c(rep(0, length(true_unique)), rep(1, length(sim_unique))))
    z <- z[order(z[, 1]), ]
    simlow  <- cumsum(z[, 3])[z[, 3] == 0]
    valid   <- (simlow > 0) & (simlow < sum(z[, 3]))
    x1 <- z[match(simlow,     cumsum(z[, 3]))[valid], 1]
    x2 <- z[match(simlow + 1, cumsum(z[, 3]))[valid], 1]
    y1 <- z[match(simlow,     cumsum(z[, 3]))[valid], 2]
    y2 <- z[match(simlow + 1, cumsum(z[, 3]))[valid], 2]
    logfdr <- rep(0, length(true_unique))
    logfdr[valid] <- (y2 - y1) * log(true_unique)[valid] / (x2 - x1) +
                     (y1 * x2 - y2 * x1) / (x2 - x1) -
                     log(cumsum(true_count) / sum(true_count))[valid]
    lowp  <- (cumsum(sim_count) / sum(sim_count)) < lmmax
    if (sum(lowp) > 1) {
        lmu <- max(sim_unique[lowp])
        if (is.finite(lmu) && min(true_unique) < min(sim_unique)) {
            lmfit <- stats::lm(
                log(cumsum(sim_count[lowp]) / sum(sim_count)) ~ log(sim_unique[lowp]))
            lf_fit <- lmfit$coefficients[2] * log(true_unique) + lmfit$coefficients[1] -
                      log(cumsum(true_count) / sum(true_count))
            interp  <- true_unique < lmu & true_unique > min(sim_unique)
            logfdr[interp]  <- (lf_fit[interp]   * (log(lmu) - log(true_unique[interp])) -
                                 logfdr[interp]   * (log(min(sim_unique)) - log(true_unique[interp]))) /
                                (log(lmu) - log(min(sim_unique)))
            nointerp <- true_unique < min(sim_unique)
            logfdr[nointerp] <- lf_fit[nointerp]
        }
    }
    logfdr <- cummax(logfdr); logfdr[logfdr > 0] <- 0
    .cell2cell_matrix(logfdr[match(true_pv, true_unique)] / log(10), cell_names)
}

#' Build an annotated hierarchical clustering tree from pin matrix and Fisher distance
#'
#' Wraps \code{stats::hclust} and annotates each merge node with:
#' the maximum and mean pairwise log10(FDR) within the node (\code{mergefdr},
#' \code{meanfdr}), the fraction of cells sharing each pin (\code{sharing}),
#' node size (\code{nodesize}), and mean pin complexity (\code{complexity}).
#'
#' @param pinmat pins x cells binary matrix (or data.frame)
#' @param mat_fdr symmetric cell x cell log10(FDR) matrix
#' @param mat_dist symmetric cell x cell distance matrix.  Uses log10(p-value)
#'   convention: more negative = more similar; this is consumed by
#'   \code{stats::hclust} where numerically smaller = merged first.
#' @param hcmethod linkage method passed to \code{stats::hclust} (default \code{"average"})
#' @return \code{hclust} object augmented with \code{mergefdr}, \code{meanfdr},
#'   \code{nodesize}, \code{complexity}, \code{leaflist}, \code{labellist},
#'   \code{sharing}, and \code{featuremat} fields
#' @keywords internal
#' @noRd
.hclust_tree <- function(pinmat, mat_fdr, mat_dist, hcmethod = "average") {
    hc  <- stats::hclust(stats::as.dist(mat_dist), method = hcmethod)
    n   <- nrow(hc$merge)
    leaflist   <- vector("list", n); labellist  <- vector("list", n)
    mergefdr   <- numeric(n);        meanfdr    <- numeric(n)
    nodesize   <- integer(n);        complexity <- numeric(n)
    sharing    <- matrix(NA_real_, nrow = nrow(pinmat), ncol = n)
    for (i in seq_len(n)) {
        li <- if (hc$merge[i, 1] < 0) -hc$merge[i, 1] else leaflist[[hc$merge[i, 1]]]
        ri <- if (hc$merge[i, 2] < 0) -hc$merge[i, 2] else leaflist[[hc$merge[i, 2]]]
        leaflist[[i]]  <- c(li, ri)
        labellist[[i]] <- hc$labels[leaflist[[i]]]
        sharing[, i]   <- rowMeans(as.matrix(pinmat)[, labellist[[i]], drop = FALSE])
        complexity[i]  <- mean(colSums(as.matrix(pinmat)[, labellist[[i]], drop = FALSE]))
        nodesize[i]    <- length(leaflist[[i]])
        sub <- mat_fdr[leaflist[[i]], leaflist[[i]]]
        mergefdr[i] <- max(sub[upper.tri(sub)])
        meanfdr[i]  <- mean(sub[upper.tri(sub)])
    }
    hc$mergefdr  <- mergefdr; hc$meanfdr   <- meanfdr
    hc$nodesize  <- nodesize; hc$complexity <- complexity
    hc$leaflist  <- leaflist; hc$labellist  <- labellist
    hc$sharing   <- sharing;  hc$featuremat <- abs(as.matrix(pinmat))
    hc
}
