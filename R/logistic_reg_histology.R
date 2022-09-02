#' Perform logistic regression for a binary/binarized trait
#'
#' These functions were developed for the genetic analysis of cancer subtypes,
#' where "histology" is considered the phenotype, and copy number the genotype.
#' The association is carried out using logistic regresion by `glm` using
#' `family = "binary"`, as this recreates the standard model of quantitative
#' genetics P = G + E.  Alternate arguments of `family` have no been tested.
#'
#' @param cnr a cnr bundle
#'
#' @param trait character, name of the trait of interest to analyze.  Must
#'  be a column in the phenotype matrix (Y). e.g. "binary1"
#'
#' @param pheno0  character, phenotype(s) to use as baseline, e.g. "0"
#'
#' @param pheno1  character, phenotype(s) to use as alternate, e.g. "1"
#'
#' @param exclude.cluster character, list of clusters to exclude, e.g.
#'  hypersegmented, Stroma, etc. Default "HC"
#'
#' @param family  character, description of the error distribution and
#' link function to be used in the model. See \link[stats]{glm} for details.
#'  Default "binomial"
#'  
#' @param na.action character, handling of NA. default is "na.exclude"
#'
#' @param ... additional arguments passed to glm
#'
#' @return
#' a CNR object with results from a logistic regression analysis
#'  (family = "binomial") with effect estimates, and p-values attached to
#'  the chromInfo and gene.index matrices.
#'
#' Results columns are "Estimate", "Std.Error", "z.value", "p.value" ; with
#'  the phenotype comparison pre-apended as <pheno0>.vs.<pheno1>.lr.<value>.
#'  Using grade as an example  these would be:
#' 
#' 0.vs.1.lr.Estimate
#' 0.vs.1.lr.Std.Error
#' 0.vs.1.lr.z.value
#' 0.vs.1.lr.p.value
#'
#' 
#' @examples
#'
#' data(cnr)
#'
#' cnr <- histo_logit(cnr, trait = "binary1",
#'    pheno0 = 0, pheno1 = 1)
#'
#' cnr <- histo_logit(cnr, trait = "category1",
#'    pheno0 = "A", pheno1 = c("B", "C"))
#'
#' cnr <- histo_logit(cnr, trait = "category2",
#'    pheno0 = c("X", "Y"), pheno1 = "Z")
#' 
#' @export
histo_logit <- function(cnr, trait, pheno0, pheno1,
                        exclude.cluster = "HC",
                        family = "binomial", na.action = "na.exclude", ...) {

    ## bin level
    cnr <- histo_logit_bin(cnr, trait = trait,
                           pheno0 = pheno0, pheno1 = pheno1,
                           exclude.cluster = exclude.cluster,
                           family = family, na.action = na.action, ...)

    ## gene level
    cnr <- histo_logit_gene(cnr, trait = trait,
                            pheno0 = pheno0, pheno1 = pheno1,
                            exclude.cluster = exclude.cluster,
                            family = family, na.action = na.action, ...)
    
    return(cnr)
}

#' perform logistic regression for a binary/binarized trait with a covaritate
#'
#' These functions were developed for the genetic analysis of cancer subtypes,
#' where "histology" is considered the phenotype, and copy number the genotype.
#' The association is carried out using logistic regresion by `glm` using
#' `family = "binary"`, as this recreates the standard model of quantitative
#' genetics P = G + E.  Alternate arguments of `family` have no been tested.
#'
#' @param cnr a cnr bundle
#'
#' @param trait character, name of the trait of interest to analyze.  Must
#'  be a column in the phenotype matrix (Y). e.g. "binary1"
#'
#' @param pheno0  character, phenotype(s) to use as baseline, e.g. "0"
#'
#' @param pheno1  character, phenotype(s) to use as alternate, e.g. "1"
#'
#' @param covar  character, model covariates to include in the model,
#'  e.g. "category1"
#'
#' @param exclude.cluster character, list of clusters to exclude, e.g.
#'  hypersegmented, Stroma, etc. Default "HC"
#'
#' @param family  character, description of the error distribution and
#'  link function to be used in the model. See \link[stats]{glm} for details.
#'  Default "binomial"
#'  
#' @param na.action character, handling of NA. default is "na.exclude"
#'
#' @param ... additional arguments passed to glm
#' 
#'
#' @return
#' a CNR object with results from a logistic regression analysis
#' (family = "binomial") with effect estimates, and p-values attached to
#' the chromInfo and gene.index matrices.
#'
#' Results columns are "Estimate", "Std.Error", "z.value", "p.value" ; with
#'  the phenotype comparison pre-apended as
#'  <pheno0>.vs.<pheno1>.cv<covar>.lr.<value>.
#'  Using grade as an example  these would be:
#' 
#' 0.vs.1..lr.Estimate
#' 0.vs.1.quantitative1.lr.Std.Error
#' 0.vs.1.quantitative1.lr.z.value
#' 0.vs.1.quantitative1.lr.p.value
#'
#' 
#' @examples
#'
#' data(cnr)
#'
#' cnr <- histo_logit_cov(cnr, trait = "binary1",
#'    pheno0 = 0, pheno1 = 1, covar = "category1")
#'
#' cnr <- histo_logit_cov(cnr, trait = "category1",
#'    pheno0 = "A", pheno1 = c("B", "C"), covar = "category2")
#'
#' cnr <- histo_logit_cov(cnr, trait = "category2",
#'    pheno0 = c("X", "Y"), pheno1 = "Z", covar = "category1")
#' 
#' @export
histo_logit_cov <- function(cnr, trait, pheno0, pheno1, covar, 
                            exclude.cluster = "HC",
                            family = "binomial", na.action = "na.exclude", ...) {
    
    cnr <- histo_logit_bin_cov(cnr, trait = trait,
                           pheno0 = pheno0, pheno1 = pheno1, covar = covar,
                           exclude.cluster = exclude.cluster,
                           family = family, na.action = na.action, ...)
    
    cnr <- histo_logit_gene_cov(cnr, trait = trait,
                            pheno0 = pheno0, pheno1 = pheno1, covar = covar,
                            exclude.cluster = exclude.cluster,
                            family = family, na.action = na.action, ...)
    
    return(cnr)
}



#' perform logistic regression for a binary/binarized trait
#' 
#' Internal to histo_logit for gene level analysis
#' 
#' @param cnr a cnr bundle
#'
#' @param trait character, name of the trait of interest to analyze.  Must
#'  be a column in the phenotype matrix (Y). e.g. "binary1"
#'
#' @param pheno0  character, phenotype(s) to use as baseline, e.g. "0"
#'
#' @param pheno1  character, phenotype(s) to use as alternate, e.g. "1"
#'
#' @param exclude.cluster character, list of clusters to exclude, e.g.
#'  hypersegmented, Stroma, etc. Default "HC"
#'
#' @param family  character, description of the error distribution and
#' link function to be used in the model. See \link[stats]{glm} for details.
#'  Default "binomial"
#'  
#' @param na.action character, handling of NA. default is "na.exclude"
#'
#' @param ... additional arguments passed to glm
#'
#' @return
#' a CNR object with results from a logistic regression analysis
#'  (family = "binomial") with effect estimates, and p-values attached to
#'  the chromInfo matrices.
#'
#' Results columns are "Estimate", "Std.Error", "z.value", "p.value" ; with
#'  the phenotype comparison pre-apended as <pheno0>.vs.<pheno1>.lr.<value>.
#'  Using grade as an example  these would be:
#' 
#' 0.vs.1.lr.Estimate
#' 0.vs.1.lr.Std.Error
#' 0.vs.1.lr.z.value
#' 0.vs.1.lr.p.value
#'
#' 
#' @examples \dontrun{
#'
#' data(cnr)
#'
#' cnr <- histo_logit(cnr, trait = "binary1",
#'    pheno0 = 0, pheno1 = 1)
#'
#' cnr <- histo_logit(cnr, trait = "category1",
#'    pheno0 = "A", pheno1 = c("B", "C"))
#'
#' cnr <- histo_logit(cnr, trait = "category2",
#'    pheno0 = c("X", "Y"), pheno1 = "Z")
#' }
#' 
#' @importFrom stats glm coef
#' 
#' @keywords internal
#' @noRd
histo_logit_gene <- function(cnr, trait, pheno0, pheno1,
                             exclude.cluster = "HC",
                             family = "binomial",
                             na.action = "na.exclude", ...) {
    
    ## ----set.phenotype----
    y <- matrix(NA, nrow = nrow(cnr$Y))
    table(cnr$Y[, trait])
    y[cnr$Y[, trait] %in% pheno0] <- 0
    y[cnr$Y[, trait] %in% pheno1] <- 1
    y[!cnr$Y[, trait] %in% c(pheno0, pheno1)] <- NA
    y[cnr$Y$final_cluster %in% exclude.cluster] <- NA
    
    ## perform genome scan with LR
    reg.eff <- t(apply(cnr$genes, 2, function(x) {
        ## perform glm
        out <- stats::glm(y ~  x, family = family, na.action = na.action, ...)
        ## exctract genotype effects
        if(nrow(stats::coef(summary(out))) == 2) {
            eff <- stats::coef(summary(out))["x", ]
        } else {
            ## if no variation on x, set effects to 0
            ## and stats to NA
            eff <- cbind("Estimate" = 0,
                         "Std. Error" = 0,
                         "z value" = NA,
                         "Pr(>|z|)" = NA)
        }
        ## return effects tables
        return(eff)
    }))
    
    ## copy to gene.index
    colnames(reg.eff) <- c("Estimate", "Std.Error", "z.value", "p.value")
    ## create new column names containing the phenotype comparison
    colnames(reg.eff) <- gsub(" ", ".",
                              paste0(paste(pheno0, collapse = "."),
                                     ".vs.",
                                     paste(pheno1, collapse = "."),
                                     ".lr.", colnames(reg.eff)))
    ## add hgnc.symbol
    reg.eff <- cbind(reg.eff, hgnc.symbol = gsub("\\.", "-", rownames(reg.eff)))
    
    ## merge with gene.index
    cnr[["gene.index"]] <- merge(cnr$gene.index, reg.eff,
                                 by = "hgnc.symbol", sort = FALSE)
    rownames(cnr[["gene.index"]]) <- cnr[["gene.index"]]$hgnc.symbol
    
    ## output cnr
    return(cnr)
} ## end histo_logit_gene


#' perform logistic regresion for a binarized trait
#'
#' Internal to histo_logit for bin level analysis
#' 
#' @param cnr a cnr bundle
#'
#' @param trait character, name of the trait of interest to analyze.  Must
#'  be a column in the phenotype matrix (Y). e.g. "binary1"
#'
#' @param pheno0  character, phenotype(s) to use as baseline, e.g. "0"
#'
#' @param pheno1  character, phenotype(s) to use as alternate, e.g. "1"
#'
#' @param exclude.cluster character, list of clusters to exclude, e.g.
#'  hypersegmented, Stroma, etc. Default "HC"
#'
#' @param family  character, description of the error distribution and
#' link function to be used in the model. See \link[stats]{glm} for details.
#'  Default "binomial"
#'  
#' @param na.action character, handling of NA. default is "na.exclude"
#'
#' @param ... additional arguments passed to glm
#'
#' @return
#' a CNR object with results from a logistic regression analysis
#'  (family = "binomial") with effect estimates, and p-values attached to
#'  the chromInfo matrices.
#'
#' Results columns are "Estimate", "Std.Error", "z.value", "p.value" ; with
#'  the phenotype comparison pre-apended as <pheno0>.vs.<pheno1>.lr.<value>.
#'  Using grade as an example  these would be:
#' 
#' 0.vs.1.lr.Estimate
#' 0.vs.1.lr.Std.Error
#' 0.vs.1.lr.z.value
#' 0.vs.1.lr.p.value
#'
#' 
#' @examples \dontrun{
#'
#' data(cnr)
#'
#' cnr <- histo_logit(cnr, trait = "binary1",
#'    pheno0 = 0, pheno1 = 1)
#'
#' cnr <- histo_logit(cnr, trait = "category1",
#'    pheno0 = "A", pheno1 = c("B", "C"))
#'
#' cnr <- histo_logit(cnr, trait = "category2",
#'    pheno0 = c("X", "Y"), pheno1 = "Z")
#'
#' }
#' 
#' @importFrom stats glm coef
#' 
#' @keywords internal
#' @noRd
histo_logit_bin <- function(cnr, trait, pheno0, pheno1, 
                            exclude.cluster = "HC",
                            family = "binomial",
                            na.action = "na.exclude", ...) {
    
    ## ----set.phenotype----
    y <- matrix(NA, nrow = nrow(cnr$Y))
    table(cnr$Y[, trait])
    y[cnr$Y[, trait] %in% pheno0] <- 0
    y[cnr$Y[, trait] %in% pheno1] <- 1
    y[!cnr$Y[, trait] %in% c(pheno0, pheno1)] <- NA
    y[cnr$Y$final_cluster %in% exclude.cluster] <- NA

    ## perform genome scan with LR
    reg.eff <- t(apply(cnr$X, 1, function(x) {
        ## perform glm
        out <- stats::glm(y ~  x, family = family, na.action = na.action, ...)
        ## exctract genotype effects
        if(nrow(stats::coef(summary(out))) == 2) {
            eff <- stats::coef(summary(out))["x", ]
        } else {
            ## if no variation on x, set effects to 0
            ## and stats to NA
            eff <- cbind("Estimate" = 0,
                         "Std. Error" = 0,
                         "z value" = NA,
                         "Pr(>|z|)" = NA)
        }
        ## return effects tables
        return(eff)
    }))
        
    ## copy to gene.index
    colnames(reg.eff) <- c("Estimate", "Std.Error", "z.value", "p.value")
    ## create new column names containing the phenotype comparison
    colnames(reg.eff) <- gsub(" ", ".",
                              paste0(paste(pheno0, collapse = "."),
                                     ".vs.",
                                     paste(pheno1, collapse = "."),
                                     ".lr.", colnames(reg.eff)))
    ## add hgnc.symbol
    reg.eff <- cbind(reg.eff, hgnc.symbol = gsub("\\.", "-", rownames(reg.eff)))

    
    ## merge with chromInfo
    cnr <- addInfo(cnr, reg.eff)    
        
    ## output cnr
    return(cnr)
} ## histo_logit_bin


#' perform logistic regresion for a binarized trait with one covariate
#'
#' Internal for histo_logit_cov for gene level analysis
#' 
#' @param cnr a cnr bundle
#'
#' @param trait character, name of the trait of interest to analyze.  Must
#'  be a column in the phenotype matrix (Y). e.g. "binary1"
#'
#' @param pheno0  character, phenotype(s) to use as baseline, e.g. "0"
#'
#' @param pheno1  character, phenotype(s) to use as alternate, e.g. "1"
#'
#' @param covar  character, model covariates to include in the model,
#'  e.g. "category1"
#'
#' @param exclude.cluster character, list of clusters to exclude, e.g.
#'  hypersegmented, Stroma, etc. Default "HC"
#'
#' @param family  character, description of the error distribution and
#'  link function to be used in the model. See \link[stats]{glm} for details.
#'  Default "binomial"
#'  
#' @param na.action character, handling of NA. default is "na.exclude"
#'
#' @param ... additional arguments passed to glm
#' 
#'
#' @return
#' a CNR object with results from a logistic regression analysis
#' (family = "binomial") with effect estimates, and p-values attached to
#' the chromInfo and gene.index matrices.
#'
#' Results columns are "Estimate", "Std.Error", "z.value", "p.value" ; with
#'  the phenotype comparison pre-apended as
#'  <pheno0>.vs.<pheno1>.cv<covar>.lr.<value>.
#'  Using grade as an example  these would be:
#' 
#' 0.vs.1..lr.Estimate
#' 0.vs.1.quantitative1.lr.Std.Error
#' 0.vs.1.quantitative1.lr.z.value
#' 0.vs.1.quantitative1.lr.p.value
#'
#' 
#' @examples \dontrun{
#'
#' data(cnr)
#'
#' cnr <- histo_logit_cov(cnr, trait = "binary1",
#'    pheno0 = 0, pheno1 = 1, covar = "category1")
#'
#' cnr <- histo_logit_cov(cnr, trait = "category1",
#'    pheno0 = "A", pheno1 = c("B", "C"), covar = "category2")
#'
#' cnr <- histo_logit_cov(cnr, trait = "category2",
#'    pheno0 = c("X", "Y"), pheno1 = "Z", covar = "category1")
#'
#' }
#' @importFrom stats glm coef
#' 
#' @keywords internal
#' @noRd
histo_logit_gene_cov <- function(cnr, trait, pheno0, pheno1, covar,
                                 exclude.cluster = "HC",
                                 family = "binomial",
                                 na.action = "na.exclude", ...) {
    ## ----set.phenotype----
    y <- matrix(NA, nrow = nrow(cnr$Y))
    table(cnr$Y[, trait])
    y[cnr$Y[, trait] %in% pheno0] <- 0
    y[cnr$Y[, trait] %in% pheno1] <- 1
    y[!cnr$Y[, trait] %in% c(pheno0, pheno1)] <- NA
    y[cnr$Y$final_cluster %in% exclude.cluster] <- NA

    cc <- cnr$Y[, covar]

    ## perform genome scan with LR
    reg.eff <- t(apply(cnr$genes, 2, function(x) {
        ## perform glm
        out <- stats::glm(y ~  x + cc, family = family,
                   na.action = na.action, ...)
        ## exctract genotype effects
        if(nrow(stats::coef(summary(out))) == 2) {
            eff <- stats::coef(summary(out))["x", ]
        } else {
            ## if no variation on x, set effects to 0
            ## and stats to NA
            eff <- cbind("Estimate" = 0,
                         "Std. Error" = 0,
                         "z value" = NA,
                         "Pr(>|z|)" = NA)
        }
        ## return effects tables
        return(eff)
    }))
        
    ## copy to gene.index
    colnames(reg.eff) <- c("Estimate", "Std.Error", "z.value", "p.value")
    ## create new column names containing the phenotype comparison
    colnames(reg.eff) <- gsub(" ", ".",
                              paste0(paste(pheno0, collapse = "."),
                                     ".vs.",
                                     paste(pheno1, collapse = "."),
                                     ".cv.", covar,
                                     ".lr.", colnames(reg.eff)))
    ## add hgnc.symbol
    reg.eff <- cbind(reg.eff, hgnc.symbol = gsub("\\.", "-", rownames(reg.eff)))

    ## merge with gene.index
    cnr[["gene.index"]] <- merge(cnr$gene.index, reg.eff,
                                 by = "hgnc.symbol", sort = FALSE)
    rownames(cnr[["gene.index"]]) <- cnr[["gene.index"]]$hgnc.symbol

    ## output cnr
    return(cnr)
} ## end histo_logit_gene_cov



#' perform logistic regresion for a binarized trait with one covariate
#'
#' Internal to histo_logit_cov for bin level analysis
#' 
#' @param cnr a cnr bundle
#'
#' @param trait character, name of the trait of interest to analyze.  Must
#'  be a column in the phenotype matrix (Y). e.g. "binary1"
#'
#' @param pheno0  character, phenotype(s) to use as baseline, e.g. "0"
#'
#' @param pheno1  character, phenotype(s) to use as alternate, e.g. "1"
#'
#' @param covar  character, model covariates to include in the model,
#'  e.g. "category1"
#'
#' @param exclude.cluster character, list of clusters to exclude, e.g.
#'  hypersegmented, Stroma, etc. Default "HC"
#'
#' @param family  character, description of the error distribution and
#'  link function to be used in the model. See \link[stats]{glm} for details.
#'  Default "binomial"
#'  
#' @param na.action character, handling of NA. default is "na.exclude"
#'
#' @param ... additional arguments passed to glm
#' 
#'
#' @return
#' a CNR object with results from a logistic regression analysis
#' (family = "binomial") with effect estimates, and p-values attached to
#' the chromInfo and gene.index matrices.
#'
#' Results columns are "Estimate", "Std.Error", "z.value", "p.value" ; with
#'  the phenotype comparison pre-apended as
#'  <pheno0>.vs.<pheno1>.cv<covar>.lr.<value>.
#'  Using grade as an example  these would be:
#' 
#' 0.vs.1..lr.Estimate
#' 0.vs.1.quantitative1.lr.Std.Error
#' 0.vs.1.quantitative1.lr.z.value
#' 0.vs.1.quantitative1.lr.p.value
#'
#' 
#' @examples \dontrun{
#'
#' data(cnr)
#'
#' cnr <- histo_logit_cov(cnr, trait = "binary1",
#'    pheno0 = 0, pheno1 = 1, covar = "category1")
#'
#' cnr <- histo_logit_cov(cnr, trait = "category1",
#'    pheno0 = "A", pheno1 = c("B", "C"), covar = "category2")
#'
#' cnr <- histo_logit_cov(cnr, trait = "category2",
#'    pheno0 = c("X", "Y"), pheno1 = "Z", covar = "category1")
#'
#' }
#' 
#' @importFrom stats glm coef
#' 
#' @keywords internal
#' @noRd
histo_logit_bin_cov <- function(cnr, trait, pheno0, pheno1, covar,
                                exclude.cluster = "HC",
                                family = "binomial",
                                na.action = "na.exclude", ...) {

    ## ----set.phenotype----
    y <- matrix(NA, nrow = nrow(cnr$Y))
    table(cnr$Y[, trait])
    y[cnr$Y[, trait] %in% pheno0] <- 0
    y[cnr$Y[, trait] %in% pheno1] <- 1
    y[!cnr$Y[, trait] %in% c(pheno0, pheno1)] <- NA
    y[cnr$Y$final_cluster %in% exclude.cluster] <- NA

    ## get covariate
    cc <- cnr$Y[, covar]
    
    ## perform genome scan with LR
    reg.eff <- t(apply(cnr$X, 1, function(x) {
        ## perform glm
        out <- stats::glm(y ~  x + cc, family = family, na.action = na.action,
                   ...)
        ## exctract genotype effects
        if(nrow(stats::coef(summary(out))) == 2) {
            eff <- stats::coef(summary(out))["x", ]
        } else {
            ## if no variation on x, set effects to 0
            ## and stats to NA
            eff <- cbind("Estimate" = 0,
                         "Std. Error" = 0,
                         "z value" = NA,
                         "Pr(>|z|)" = NA)
        }
        ## return effects tables
        return(eff)
    }))
        
    ## copy to gene.index
    colnames(reg.eff) <- c("Estimate", "Std.Error", "z.value", "p.value")
    ## create new column names containing the phenotype comparison
    colnames(reg.eff) <- gsub(" ", ".",
                              paste0(paste(pheno0, collapse = "."),
                                     ".vs.",
                                     paste(pheno1, collapse = "."),
                                     ".cv.", covar,
                                     ".lr.", colnames(reg.eff)))
    ## add hgnc.symbol
    reg.eff <- cbind(reg.eff, hgnc.symbol = gsub("\\.", "-", rownames(reg.eff)))
    
    
    ## merge with chromInfo
    cnr <- addInfo(cnr, reg.eff)    
    
    ## output cnr
    return(cnr)
} ## histo_logit_bin_cov


#' Manhattan plot of logistic regression output
#'
#' @param cnr a cnr bundle
#'
#' @param pval character, name of the column in the phenotype matrix (Y)
#' containing p-values
#'
#' @param sig.threshold numeric, significance threshold to use
#'
#' @param pch plot character, default 21
#'
#' @param cex character, expansion, default 0.64
#'
#' @param ylab character, y-axis label
#'
#' @param xlab character, x-axis label
#'
#' @param ... additional arguments passed to `plot`
#'
#' @return
#' A manhattan plot with p-values for e.g. logistic regresion or other comparisons
#'
#' @export
plot_lr <- function(cnr, pval = "LL.vs.SCL.lr.p.value",
                    sig.threshold = 10^-8,
                    pch = 21, cex = 0.64,
                    ylab = expression(-Log[10]~P-value), xlab = "HSA Genome",
                    ...) {
    ## set chr breaks and michr
    chrBreaks <- cumsum(table(cnr$chromInfo$bin.chr))
    midChr <- chrBreaks - floor((chrBreaks - c(1, chrBreaks[1:(length(chrBreaks) - 1)]))/2)
    ## manhattan plot by bin
    plot(-log10(cnr$chromInfo[, pval]),
         pch = pch, cex = cex, , xaxt = "n", xaxs = "i",
         ylab = ylab, xlab = xlab,
         ...)
    abline(h = -log10(sig.threshold), col = "#D40000")
    abline(v = chrBreaks, lty = 2, lwd = 0.8, col = "gray40")
    axis(1, at = midChr, labels = names(midChr))
}  ## plot_lr



#' Export the p-values for comparisons of interest to visualize in IGV
#' 
#' @param cnr a cnr bundle
#'
#' @param outdir output directory, default is "."
#'
#' @param pval.column column to export
#'
#' @param outfile.prefix basename of file
#'
#' @param extension extesion compatible with IGV. Any of ".linear",
#'  ".logistic", ".assoc", ".qassoc", ".gwas".  Default ".logistic"
#'
#' @param na.value value to give to NA, default 1 (i.e. p-value = 1)
#'
#' @param print.output print output on console as an object list with
#'  bin.pvalues and gene.pvalues
#'
#' @return
#' 
#' Exports p-values from an analysis into two files, one containing "bin" level
#' p-values, and another containing "gene" level p-values in IGV's GWAS format.
#'
#' Bin-level p-values are expanded to all genes within each bin, such that
#' the user can visualize the bin blocks with added resolution.
#'
#' By default, any p-values with  NA are converted to 1.  These can arise in
#' monomorphic regions of the genome and don't have an association to
#' the phenotype.  
#' 
#' 
#' @references https://software.broadinstitute.org/software/igv/GWAS
#'
#' @importFrom dplyr %>% mutate select
#' 
#' @export
export_pval_igv <- function(cnr, pval.column, outdir = ".",
                            outfile.prefix, extension = ".logistic",
                            na.value = "1", print.output = FALSE ) {
    ## get bin pvalues
    pvals1 <- cnr$chromInfo[, c("bin.id", "bin.chrom", pval.column)]
    names(pvals1) <- c("bin.id", "chr", "pval")
    ## expand to gene ID
    column.keep <- c("bin.id", "start", "end", "hgnc.symbol")
    pvals1 <- merge(cnr$gene.index[, column.keep],
                    pvals1, by = "bin.id")
    ## wrangle & write
    pvals1 %>%
        dplyr::mutate(bp = floor((.data$start + .data$end)/2),
               snp = .data$hgnc.symbol) %>%
        dplyr::select(.data$chr, .data$bp, .data$snp, .data$pval) %>%
        write.table(
            file = file.path(outdir, paste0(outfile.prefix, "_bin_pval",
                                            extension)),
            sep = "\t", quote = FALSE, row.names = FALSE, na = na.value)
    ## get gene pvalues
    gene.column.keep <- c("seqnames", "start", "end", "hgnc.symbol", pval.column)
    genePvals <- cnr$gene.index[, gene.column.keep]
    names(genePvals) <- c("chr", "start", "end", "hgnc.symbol", "pval")
    ## wrangle & write
    genePvals %>%
        dplyr::mutate(bp = floor((.data$start + .data$end)/2),
               snp = .data$hgnc.symbol) %>%
        dplyr::select(.data$chr, .data$bp, .data$snp, .data$pval) %>%
        write.table(
            file = file.path(outdir,
                             paste0(outfile.prefix, "_gene_pval",
                                    extension)),
            sep = "\t", quote = FALSE, row.names = FALSE, na = na.value)
    if(!print.output) {
        out <- NULL
    } else {
        out <- list(interpolated.bin.pvalues = pvals1, gene.pvalues = genePvals)
    }

    return(out)
}

#' Wrapper to joint CNA effect estimation using GLMNet
#'
#' Function is still in development
#'
#' @param cnr a cnr bundle
#'
#' @param trait character, name of the trait of interest to analyze.  Must
#'  be a column in the phenotype matrix (Y). e.g. "binary1"
#'
#' @param pheno0  character, phenotype(s) to use as baseline, e.g. "0"
#'
#' @param pheno1  character, phenotype(s) to use as alternate, e.g. "1"
#'
#' @param exclude.cluster character, list of clusters to exclude, e.g.
#'  hypersegmented, Stroma, etc. Default "HC"
#'
#' @return
#' A CNR object containing joint CNA effect estimates using GLMNet.  Analysis
#' is carried out in both `X` and `genes` copy number matrices.
#'
#' Effects are estimated at two Lambda; min, 1se.  Lambda values are appended
#' as attributes to the chromInfo, and gene.index, for the bin, respectively.
#' 
#' @examples
#' data(cnr)
#'
#' cnr <- estimate_joint_effects(cnr, trait = "binary1", pheno0 = 0, pheno1 = 1)
#' 
#' @export
estimate_joint_effects <- function(cnr, trait, pheno0, pheno1,
                                   exclude.cluster = "HC") {

    cnr <- estimate_joint_effects_bin(cnr = cnr, trait = trait,
                                      pheno0 = pheno0, pheno1 = pheno1,
                                      exclude.cluster = exclude.cluster)

    cnr <- estimate_joint_effects_gene(cnr = cnr, trait = trait,
                                      pheno0 = pheno0, pheno1 = pheno1,
                                      exclude.cluster = exclude.cluster)

    return(cnr)

}

#' Joint CNA effect estimation using GLMNet at the bin level
#'
#' Internal function to \link{estimate_joint_effects}
#'
#' @param cnr a cnr bundle
#'
#' @param trait character, name of the trait of interest to analyze.  Must
#'  be a column in the phenotype matrix (Y). e.g. "binary1"
#'
#' @param pheno0  character, phenotype(s) to use as baseline, e.g. "0"
#'
#' @param pheno1  character, phenotype(s) to use as alternate, e.g. "1"
#'
#' @param exclude.cluster character, list of clusters to exclude, e.g.
#'  hypersegmented, Stroma, etc. Default "HC"
#'
#' @return
#' A CNR object containing joint CNA effect estimates using GLMNet at the
#' `BIN` level.  Hence, this function runs the analysis on the `X` matrix.
#'
#' Effects are estimated at two Lambda; min, 1se.  Lambda values are appended
#' as attributes to the chromInfo, and gene.index, for the bin, respectively.
#' 
#' @examples \dontrun {
#' data(cnr)
#'
#' cnr <- estimate_joint_effects(cnr, trait = "binary1",
#'    pheno0 = 0, pheno1 = 1)
#'
#' }
#' 
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats complete.cases predict
#' @keywords internal
#' @noRd
estimate_joint_effects_bin <- function(cnr, trait, pheno0, pheno1,
                                   exclude.cluster = "HC") {

    ## set phenotype
    y <- matrix(NA, nrow = nrow(cnr$Y))
    table(cnr$Y[, trait])
    y[cnr$Y[, trait] %in% pheno0] <- 0
    y[cnr$Y[, trait] %in% pheno1] <- 1
    y[!cnr$Y[, trait] %in% c(pheno0, pheno1)] <- NA
    y[cnr$Y$final_cluster %in% exclude.cluster] <- NA

    ## remove missing values and create genotype
    ## and phenotype matrices
    P <- y[!is.na(y)]
    G <- cnr$X[, !is.na(y)]

    ## effect estimation
    fit = glmnet::glmnet(t(G), P, family = "binomial")
    cvfit = glmnet::cv.glmnet(t(G), P, family = "binomial",
                              type.measure = "class")
    beta_glmnet = as.matrix(stats::predict(fit, type = "coefficients")[-1,])
    
    ## plot(fit, xvar = "dev", label = TRUE)
    ## plot(cvfit)
    ## print(cvfit)
    beta_coeff <- as.matrix(cbind(coef(fit, s = cvfit$lambda.min),
                        coef(fit, s = cvfit$lambda.1se)))

    prefix <- paste0(paste(pheno0, collapse = "."), ".vs.",
                     paste(pheno1, collapse = "."))
    
    colnames(beta_coeff) <- gsub(" ", ".",
                                 paste0(prefix,
                                        ".lr.Effect.at_Lambda.",
                                        c("min", "1se")))
    
    ## cbind but remove intercept
    cnr <- addInfo(cnr, df = beta_coeff[-1,])
    
    attr(cnr$chromInfo, paste0(prefix, ".cnvfit")) <-
                            c("lambda.min" = cvfit$lambda.min,
                              "lambda.1se" = cvfit$lambda.1se)
    return(cnr)
}



#' Joint CNA effect estimation using GLMNet at the gene level
#' 
#' Internal function to \link{estimate_joint_effects}
#'
#' @param cnr a cnr bundle
#'
#' @param trait character, name of the trait of interest to analyze.  Must
#'  be a column in the phenotype matrix (Y). e.g. "binary1"
#'
#' @param pheno0  character, phenotype(s) to use as baseline, e.g. "0"
#'
#' @param pheno1  character, phenotype(s) to use as alternate, e.g. "1"
#'
#' @param exclude.cluster character, list of clusters to exclude, e.g.
#'  hypersegmented, Stroma, etc. Default "HC"
#'
#' @return
#' A CNR object containing joint CNA effect estimates using GLMNet at the
#' `GENE` level.  Hence, this function runs the analysis on the `genes`
#' matrix.
#'
#' Effects are estimated at two Lambda; min, 1se.  Lambda values are appended
#' as attributes to the chromInfo, and gene.index, for the bin, respectively.
#' 
#' @examples \dontrun{
#' data(cnr)
#'
#' cnr <- estimate_joint_effects(cnr, trait = "binary1",
#'    pheno0 = 0, pheno1 = 1)
#'
#' }
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats complete.cases predict
#' 
#' @keywords internal
#' @noRd
estimate_joint_effects_gene <- function(cnr, trait, pheno0, pheno1,
                                        exclude.cluster = "HC") {

    ## set phenotype
    y <- matrix(NA, nrow = nrow(cnr$Y))
    table(cnr$Y[, trait])
    y[cnr$Y[, trait] %in% pheno0] <- 0
    y[cnr$Y[, trait] %in% pheno1] <- 1
    y[!cnr$Y[, trait] %in% c(pheno0, pheno1)] <- NA
    y[cnr$Y$final_cluster %in% exclude.cluster] <- NA

    ## remove missing values and create genotype
    ## and phenotype matrices
    P <- y[!is.na(y)]
    genes.all.ok <- stats::complete.cases(t(cnr$genes))
    G <- as.matrix(cnr$genes[!is.na(y), genes.all.ok])

    ## effect estimation
    fit = glmnet::glmnet(G, P, family = "binomial")
    cvfit = glmnet::cv.glmnet(G, P, family = "binomial",
                              type.measure = "class")
    beta_glmnet = as.matrix(stats::predict(fit, type = "coefficients")[-1,])
    
    ## plot(fit, xvar = "dev", label = TRUE)
    ## plot(cvfit)
    ## print(cvfit)
    beta_coeff <- as.matrix(cbind(coef(fit, s = cvfit$lambda.min),
                        coef(fit, s = cvfit$lambda.1se)))
    beta_coeff <- as.data.frame(beta_coeff)

    prefix <- paste0(paste(pheno0, collapse = "."),
                     ".vs.",
                     paste(pheno1, collapse = "."))
        
        
    colnames(beta_coeff) <- gsub(" ", ".",
                                 paste0(prefix,
                                        ".lr.Effect.at_Lambda.",
                                        c("min", "1se")))
                                        
    
    beta_coeff$hgnc.symbol <- gsub("\\.", "-", rownames(beta_coeff))
    
    cnr[["gene.index"]] <- merge(cnr$gene.index, beta_coeff[-1,],
                                 by = "hgnc.symbol",
                                 all.x = TRUE, sort = FALSE)
    
    cnr[["gene.index"]] <- cnr[["gene.index"]][cnr$gene.index$hgnc.symbol != "",]
    rownames(cnr[["gene.index"]]) <- cnr[["gene.index"]]$hgnc.symbol
    
    attr(cnr$gene.index, paste0(prefix, ".cnvfit")) <-
                             c("lambda.min" = cvfit$lambda.min,
                               "lambda.1se" = cvfit$lambda.1se)
    
    return(cnr)
}


#' Genome-wide effects plot
#'
#' Plot estimated effects genome wide.
#'
#' @param cnr a cnr bundle
#'
#' @param effect.column character, name of the effect column to plot
#'
#' @param type character, plot type, defaults to "h". See \link[graphics]{plot}
#'
#' @param ylab character, y-axis label
#'
#' @param xlab character, x-axis lable
#'
#' @param ... additional arguments passed to plot
#'
#' @return
#'
#' base graphics histogram like plot containing the estimated genome-side
#' genotypic effects.
#'
#' @examples
#' data(cnr)
#'
#' cnr <- histo_logit(cnr, trait = "binary1",
#'    pheno0 = 0, pheno1 = 1)
#' 
#' plot_effect(cnr, effect.column = "0.vs.1.lr.Estimate")  ## effect from logistic reg
#'
#' @export
plot_effect <- function(cnr,
                        effect.column = grep("Estimate", names(cnr$Y), value = TRUE)[1],
                        type = "h", 
                        ylab = expression(hat(italic(Beta))),
                        xlab = "HSA Genome",
                        ...) {

    ## set chr breaks and michr
    chrBreaks <- cumsum(table(cnr$chromInfo$bin.chr))
    midChr <- chrBreaks - floor((chrBreaks - c(1, chrBreaks[1:(length(chrBreaks) - 1)]))/2)

    ## histogram like effect plot
    plot(cnr$chromInfo[, effect.column],
         type = type,
         xaxt = "n", xaxs = "i",
         ylab = ylab, xlab = xlab,
         ...)
    abline(v = chrBreaks, lty = 2, lwd = 0.8, col = "gray40")
    axis(1, at = midChr, labels = names(midChr))
}  ## plot_effects

