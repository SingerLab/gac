## Pin matrix construction
## Internalized from SCclust (Krasnitz et al., MIT License)
## https://github.com/KrasnitzLab/SCclust

#' Convert chromosome labels to numeric indices
#' @param chrom character or numeric chromosome vector (e.g. "chr1" or 1)
#' @return numeric vector
#' @keywords internal
#' @noRd
.chrom_numeric <- function(chrom) {
    if (is.numeric(chrom)) return(as.numeric(chrom))
    n <- length(unique(chrom))
    x <- substring(chrom, 4)
    x[chrom == "chrX"] <- n - 1
    x[chrom == "chrY"] <- n
    as.numeric(x)
}

#' Mode of modes across chromosome-stratified copy number values
#' @param v numeric copy number values
#' @param otherlabel chromosome labels (same length as v)
#' @param tiebreaker value used to break ties
#' @param tiebreakerside "greater" or "lesser"
#' @return scalar numeric mode
#' @keywords internal
#' @noRd
.modeofmodes <- function(v, otherlabel, tiebreaker, tiebreakerside) {
    countmat  <- as.matrix(tapply(v, list(as.factor(v), as.factor(otherlabel)), length))
    modes     <- dimnames(countmat)[[1]][apply(countmat, 2, which.max)]
    modecounts <- tapply(modes, as.factor(modes), length)
    topmodes  <- as.numeric(names(modecounts))[modecounts == max(modecounts)]
    untie     <- topmodes[abs(topmodes - tiebreaker) == min(abs(topmodes - tiebreaker))]
    if (tiebreakerside == "greater") return(max(untie))
    return(min(untie))
}

#' Estimate per-cell ploidy from rounded copy number segments
#' @param gc_df genomic coordinate data.frame with a numeric \code{chrom} column
#' @param segment_df copy number segment matrix (bins x cells), first 3 cols are coords
#' @param rounded_df optional pre-rounded integer matrix; computed from segment_df if NULL
#' @param chromrange integer vector of chromosome indices; last two are X and Y
#' @return matrix of ploidy estimates per cell (ploidymed, ploidymod, ploidychromod, homoloss)
#' @keywords internal
#' @noRd
.calc_ploidies_internal <- function(gc_df, segment_df, rounded_df = NULL,
                                    chromrange = 1:24) {
    last_auto <- utils::tail(chromrange, 1) - 2
    if (is.null(rounded_df)) rounded_df <- round(as.matrix(segment_df[, -(1:3)]))
    getmode   <- function(x) as.numeric(names(which.max(tapply(x, as.factor(x), length))))
    auto_idx  <- gc_df[, "chrom"] <= last_auto
    ploidymod     <- apply(rounded_df[auto_idx, ], 2, getmode)
    ploidymed     <- apply(rounded_df[auto_idx, ], 2, stats::median)
    ploidychromod <- apply(rounded_df[auto_idx, ], 2, .modeofmodes,
                           otherlabel = gc_df[auto_idx, "chrom"],
                           tiebreaker = 2, tiebreakerside = "greater")
    homoloss <- colSums(!rounded_df[gc_df[, "chrom"] < last_auto, ]) /
                sum(gc_df[, "chrom"] < last_auto)
    cbind(ploidymed, ploidymod, ploidychromod, homoloss)
}

#' Augment genomic coordinate data.frame with absolute positions
#' @param gc_df chromInfo with chrom.numeric, bin.start, bin.end columns
#' @param df segment data.frame with abspos column (same row order as gc_df)
#' @return data.frame with chrom, chromstart, chromend, absstart, absend
#' @keywords internal
#' @noRd
.augment_gc <- function(gc_df, df) {
    aug <- cbind(gc_df[, c("chrom.numeric", "bin.start", "bin.end")], df[, "abspos"])
    colnames(aug) <- c("chrom", "chromstart", "chromend", "absstart")
    aug <- cbind(aug, aug[, "absstart"] + (aug[, "chromend"] - aug[, "chromstart"]))
    colnames(aug)[5] <- "absend"
    aug
}

#' Compress copy number matrix into run-length encoded segments
#' @param gc_df augmented genomic coordinate data.frame (from \code{.augment_gc})
#' @param segment_df bins x cells copy number matrix with 3 coordinate columns prepended
#' @param homoloss maximum allowed fraction of homozygous-loss bins per cell
#' @param chromrange integer chromosome index range
#' @return list: short (segment table), ploidies, cells (passing cell IDs)
#' @keywords internal
#' @noRd
.calc_segments_short <- function(gc_df, segment_df, homoloss = 1, chromrange = 1:24) {
    a  <- round(as.matrix(segment_df[, -(1:3)]))
    b  <- a - rbind(matrix(0, 1, ncol(a)), a[-nrow(a), ])
    bb <- b != 0
    bb[match(unique(segment_df[, "chrom"]), segment_df[, "chrom"]), ] <- TRUE
    segstarts <- row(a)[bb]; segvals <- a[bb]
    profid    <- dimnames(a)[[2]][cumsum(segstarts == 1)]
    segends   <- c(segstarts[-1] - 1, nrow(a)); segends[segends == 0] <- nrow(a)
    segloc    <- cbind(gc_df[segstarts, c("chrom", "chromstart", "absstart")],
                       gc_df[segends,   c("chromend", "absend")])
    tshort    <- data.frame(I(profid), segloc, segstarts, segends,
                            segbins = segends - segstarts + 1, segvals)
    ploidies_df <- .calc_ploidies_internal(gc_df, segment_df, a, chromrange)
    tshort  <- cbind(tshort,
                     cvals = tshort[, "segvals"] - ploidies_df[tshort[, "profid"], "ploidychromod"])
    good    <- ploidies_df[, "homoloss"] <= homoloss
    list(short    = tshort[tshort[, "profid"] %in% rownames(ploidies_df)[good], ],
         ploidies = ploidies_df,
         cells    = rownames(ploidies_df[good, ]))
}

#' Compute centromeric segment indices to exclude from pin construction
#' @param short_df segment table from \code{.calc_segments_short}
#' @param dropareas data.frame with chrom, from, to columns (from \code{.pin_centroareas})
#' @return logical vector (TRUE = centromeric, to be dropped)
#' @keywords internal
#' @noRd
.calc_censored_index <- function(short_df, dropareas) {
    censored <- (short_df[, "chromstart"] >= dropareas[short_df[, "chrom"], "from"]) &
                (short_df[, "chromend"]   <= dropareas[short_df[, "chrom"], "to"])
    extra <- (which(censored) + 1); extra <- extra[extra <= nrow(short_df)]
    censored[extra] <- ((short_df[extra, "chrom"]  == short_df[extra - 1, "chrom"]) &
                        (short_df[extra, "profid"] == short_df[extra - 1, "profid"])) |
                       censored[extra]
    censored
}

#' Smear segment breakpoints into overlapping windows
#' @param short_df segment table
#' @param censored logical vector of centromeric segments to exclude (or NULL)
#' @param smear number of bins to extend each breakpoint in both directions
#' @param keepboundaries logical; whether to keep chromosome boundary breakpoints
#' @param chromrange integer chromosome index range; last two are X and Y
#' @return data.frame of smeared breakpoints with profid, chrom, bpstart, bpend, bpsign
#' @keywords internal
#' @noRd
.calc_smear_breakpoints <- function(short_df, censored = NULL, smear = 1,
                                    keepboundaries = FALSE, chromrange = 1:24) {
    last_auto <- utils::tail(chromrange, 1) - 2
    dtshort <- cbind(short_df[, c("profid", "chrom")],
                     bpstart = short_df[, "segstarts"] - smear,
                     bpend   = short_df[, "segstarts"] + smear,
                     bpsign  = sign(short_df[, "cvals"] -
                                    c(0, short_df[-nrow(short_df), "cvals"])))
    ustart <- short_df[match(unique(short_df[, "chrom"]), short_df[, "chrom"]), "segstarts"]
    uend   <- c((ustart - 1)[-1], short_df[nrow(short_df), "segends"])
    if (keepboundaries) {
        hit <- (dtshort[, "bpstart"] + smear) == ustart[dtshort[, "chrom"]]
        dtshort[hit, "bpsign"] <- 2 * dtshort[hit, "bpsign"]
        if (!is.null(censored)) dtshort <- dtshort[!censored, ]
    } else {
        keep <- (dtshort[, "bpstart"] + smear) > ustart[dtshort[, "chrom"]]
        if (!is.null(censored)) keep <- keep & !censored
        dtshort <- dtshort[keep, ]
        lo <- dtshort[, "bpstart"] < ustart[dtshort[, "chrom"]]
        dtshort[lo, "bpstart"] <- ustart[dtshort[lo, "chrom"]]
        hi <- dtshort[, "bpend"] > uend[dtshort[, "chrom"]]
        dtshort[hi, "bpend"] <- uend[dtshort[hi, "chrom"]]
    }
    dtshort[dtshort[, "chrom"] %in% seq_len(last_auto), ]
}

#' Containment indicator matrix for breakpoint interval matching
#' @param vstart,vend start/end positions of reference intervals (pins)
#' @param wstart,wend start/end positions of query intervals (cell breakpoints)
#' @return two-column integer matrix (startpin, endpin) assigning each query to a pin
#' @keywords internal
#' @noRd
.containment_indicator <- function(vstart, vend, wstart, wend) {
    lv <- length(vstart); lw <- length(wstart)
    z <- cbind(c(vend, wend), c(seq_len(lv), rep(0, lw)), c(rep(0, lv), seq_len(lw)))
    z <- z[order(z[, 1]), ]
    endbeforeend <- cummax(z[, 2])[order(z[, 3])][sort(z[, 3]) != 0]
    z <- cbind(c(wstart, vstart), c(rep(lv + 1, lw), seq_len(lv)), c(seq_len(lw), rep(0, lv)))
    z <- z[order(z[, 1]), ]
    startafterstart <- rev(cummin(rev(z[, 2])))[order(z[, 3])][sort(z[, 3]) != 0]
    cbind(startafterstart, endbeforeend)
}

#' Build binary pin matrix from smeared breakpoints
#' @param short_df segment table
#' @param smear_df smeared breakpoints (from \code{.calc_smear_breakpoints})
#' @return list: pinmat (data.frame, pins x cells), pins (data.frame, bin + sign per pin)
#' @keywords internal
#' @noRd
.calc_pinmat_short <- function(short_df, smear_df) {
    allsigns <- sort(unique(smear_df[, "bpsign"]))
    pinmat <- NULL; pins <- NULL
    for (vsign in allsigns) {
        a <- smear_df[smear_df[, "bpsign"] == vsign, ]
        a <- a[order(a[, "bpend"]), ]
        apins <- NULL
        while (nrow(a) > 0) {
            apins <- c(apins, a[1, "bpend"])
            a <- a[!(a[, "bpstart"] <= apins[length(apins)] &
                     a[, "bpend"]   >= apins[length(apins)]), , drop = FALSE]
        }
        a  <- smear_df[smear_df[, "bpsign"] == vsign, ]
        ci <- .containment_indicator(apins, apins,
                                     a[order(a[, "bpend"]), "bpstart"],
                                     a[order(a[, "bpend"]), "bpend"])
        a <- cbind(a[order(a[, "bpend"]), ], ci)
        colnames(a)[(ncol(a) - 1):ncol(a)] <- c("startpin", "endpin")
        apm <- matrix(0, nrow = length(apins) + 2,
                      ncol = length(unique(short_df[, "profid"])),
                      dimnames = list(NULL, unique(short_df[, "profid"])))
        for (id in unique(a[, "profid"])) {
            apm[a[a[, "profid"] == id, "startpin"] + 1, id] <- 1
            apm[a[a[, "profid"] == id, "endpin"]   + 2, id] <-
                apm[a[a[, "profid"] == id, "endpin"] + 2, id] - 1
        }
        apm   <- apply(apm, 2, cumsum)[-c(1, length(apins) + 2), , drop = FALSE]
        apins_df <- data.frame(bin = apins, sign = rep(vsign, length(apins)))
        pinmat <- rbind(pinmat, apm)
        pins   <- rbind(pins, apins_df)
    }
    keep <- rowSums(pinmat) < ncol(pinmat)
    list(pinmat = data.frame(pinmat[keep, , drop = FALSE]),
         pins   = pins[keep, , drop = FALSE])
}

#' Build pin matrix from integer copy number profiles
#'
#' Converts a copy number matrix into a binary pins-by-cells matrix where
#' each row is a consensus breakpoint interval (pin) and each cell is 1 if
#' it carries that breakpoint. Centromeric regions are excluded via
#' \code{dropareas}.
#'
#' @param gc_df chromInfo with chrom.numeric, bin.start, bin.end, abspos columns
#' @param segment_df copy number data.frame: first 3 cols are chrom coords, remaining are cells
#' @param homoloss maximum fraction of homozygous-loss bins allowed per cell (default 1 = no filter)
#' @param dropareas centromere exclusion table from \code{.pin_centroareas} (or NULL)
#' @param smear breakpoint smear window in bins (default 1)
#' @param chromrange integer chromosome index vector; last two assumed X and Y
#' @param keepboundaries logical; keep chromosome boundary breakpoints
#' @return list: pinmat, pins, cells, ploidies
#' @keywords internal
#' @noRd
.calc_pinmat <- function(gc_df, segment_df, homoloss = 1, dropareas = NULL,
                         smear = 1, chromrange = 1:24, keepboundaries = FALSE) {
    aug      <- .augment_gc(gc_df, segment_df)
    res      <- .calc_segments_short(aug, segment_df, homoloss = homoloss,
                                     chromrange = chromrange)
    censored <- if (!is.null(dropareas))
                    .calc_censored_index(res$short, dropareas) else NULL
    smear_df <- .calc_smear_breakpoints(res$short, censored = censored, smear = smear,
                                        keepboundaries = keepboundaries,
                                        chromrange = chromrange)
    res2 <- .calc_pinmat_short(res$short, smear_df)
    list(pinmat = res2$pinmat, pins = res2$pins,
         cells  = res$cells,  ploidies = res$ploidies)
}

#' Identify centromere regions from a UCSC cytoBand table
#' @param cyto UCSC cytoBand data.frame (chrom, chromStart, chromEnd, name, gieStain)
#' @param centromere character(2) patterns matching centromere cytoband names
#'   (default \code{c("p11", "q11")})
#' @return data.frame with columns chrom (numeric), from, to
#' @keywords internal
#' @noRd
.pin_centroareas <- function(cyto, centromere = c("p11", "q11")) {
    cyto[, 1] <- .chrom_numeric(cyto[, 1])
    cyto <- cyto[order(cyto[, 1]), ]
    cl   <- cyto[grep(centromere[1], cyto[, 4]), ]
    cr   <- cyto[grep(centromere[2], cyto[, 4]), ]
    cl   <- cl[match(unique(cl[, 1]), cl[, 1]), ]
    cr   <- cr[nrow(cr):1, ]
    cr   <- cr[match(unique(cr[, 1]), cr[, 1]), ]
    cr   <- cr[nrow(cr):1, ]
    out  <- cbind(cl[, c(1, 2)], cr[, 3])
    colnames(out) <- c("chrom", "from", "to")
    out
}

#' Map genomic regions to bin row indices in chromInfo
#' @param gc_df chromInfo data.frame; must have \code{chrom.numeric}, \code{bin.start},
#'   \code{bin.end} columns
#' @param regions data.frame with chrom, from, to columns (from \code{.pin_centroareas})
#' @return integer vector of row indices in gc_df falling within any region
#' @keywords internal
#' @noRd
.pin_regions2bins <- function(gc_df, regions) {
    assertthat::assert_that(!is.null(gc_df$chrom.numeric))
    bins <- vector("list", nrow(regions))
    for (i in seq_len(nrow(regions))) {
        r  <- regions[i, ]
        df <- gc_df[gc_df$chrom.numeric == r$chrom, ]
        df <- df[((df$bin.start >= r$from) & (df$bin.start <= r$to)) |
                 ((df$bin.end   >= r$from) & (df$bin.end   <= r$to)), ]
        bins[[i]] <- rownames(df)
    }
    as.numeric(unlist(bins))
}
