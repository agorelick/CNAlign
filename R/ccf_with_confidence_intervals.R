# Calculcate cancer cell fractions and 95% confidence intervals (CCFs) for tumor mutations with given:
# - purity (cellularity)
# - number of ref, and alt reads
# - major and minor copy number (larger and smaller allele specific integer copy numbers) 
# - normal copy number (number of copies at given locus in normal human genome)
#
# code is based on this unpublished MSKCC github repo:  https://github.com/mskcc/facets-suite/blob/master/R/ccf-annotate-maf.R
# concepts are based on the referenced PMIDs (28270531, 25877892)
#
# input data: tab-delimited table with fields mutation_id, sample_id, purity, t_ref_count, t_alt_count, major_cn, minor_cn, normal_cn
# output data: same as above with added fields: t_depth, total_copies, mutant_copies, ccf, ccf_lwr95, ccf_upr95


##' expected_mutant_copies
##'
##' Estimate mutant copy number, given observed VAF, purity, and local ploidy. Based on PMID 28270531
##'
##' @export
expected_mutant_copies = function(t_var_freq, total_copies, purity, round.mcn=T) {
    if (is.na(total_copies)) {
        as.integer(NA)
    } else {
        if (total_copies == 0) total_copies = 1
        mu = t_var_freq * (1 / purity) * (purity * total_copies + (1 - purity) * 2)
        alt_copies = ifelse(mu < 1, 1, abs(mu)) # mu < 1 ~ 1, mu >= 1 ~ abs(mu)
        alt_copies = ifelse(alt_copies > total_copies, total_copies, alt_copies) # don't allow expected number greater than the total
       ifelse(round.mcn, round(alt_copies), alt_copies)
    }
}


##' estimate_ccf
##'
##' Estimate most likely CCF given observed VAF, purity and local ploidy. Based on PMID 25877892
##'
##' @export
estimate_ccf = function(purity, total_copies, mutant_copies, t_alt_count, t_depth, resolution=0.001) {

    expected_vaf = function(ccf, purity, total_copies) {
        purity * ccf * mutant_copies / (2 * (1 - purity) + purity * total_copies)
    }

    ccfs = seq(0.000, 1, resolution)
    probs = sapply(ccfs, function(c) {
                   stats::dbinom(t_alt_count, t_depth, expected_vaf(c, purity, total_copies)) })
    probs = probs / sum(probs)
    ccf_mostlikely_index = which.max(probs)
    ccf_value <- ccfs[ccf_mostlikely_index]

    ## empirical CDF to determine 95% CIs for the CCF
    CDF <- cumsum(probs) 

    # closest CCF outside the lower bound of the 95% CI
    below95 <- which(CDF <= 0.025)
    ccf_lwr95_index <- max(below95) 
    ccf_lwr95 <- ccfs[ccf_lwr95_index]

    # closest CCF outside the upper bound of the 95% CI
    above95 <- which(CDF >= 0.975)
    ccf_upr95_index <- min(above95)
    ccf_upr95 <- ccfs[ccf_upr95_index]

    list(ccf=ccf_value, ccf_lwr95=ccf_lwr95, ccf_upr95=ccf_upr95)
}



