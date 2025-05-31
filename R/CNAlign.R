##' CNAlign
##'
##' Determine purity and ploidy values for multiple tumor samples with shared ancestry. This function uses the GuRoBi solver to determine purity/ploidy values for each sample that will *maximize* the number of segments with the same (allele-specific) integer copy numbers in at least rho% of samples. 
##'
##' @export
CNAlign <- function(dat, license_path, min_ploidy=1.7, max_ploidy=6.5, min_purity=0.05, max_purity=0.95, min_aligned_seg_mb=5.0, max_homdel_mb=100.0, 
                    delta_tcn_to_int=0.2, delta_tcn_to_avg=0.1, delta_tcnavg_to_int=0.1, delta_mcn_to_int=0.2, delta_mcn_to_avg=0.1, delta_mcnavg_to_int=0.1, 
                    mcn_weight=0.5, rho=1.0, timeout=120, min_cna_segments_per_sample=5, obj2_clonalonly=0) {

    # ensure dat is a data.frame and available in the inner call
    stopifnot(is.data.frame(dat))

    # recode NA BAFs to -9
	dat$BAF[is.na(dat$BAF)] <- -9
	
    start_time <- now()
	message('Started CNAlign at ',as.character(start_time))

    out <- basilisk::basiliskRun(env = cnalign_env, fun = function(dat, license_path, ...) {
                              library(reticulate)

                              # Source the Python file
                              script_path <- system.file("python", "CNAlign.py", package = "CNAlign")
                              source_python(script_path)

                              # Run the Python function
                              result <- do_CNAlign(dat = dat, gurobi_license = license_path, ...)
                              
                              # collect all the input arguments
                              params <- list(...)
                                
                              # return the result and the params as a named-list
                              return(list(result=result, params=params))

                    }, dat = dat, license_path = license_path, 
                    min_ploidy=min_ploidy, max_ploidy=max_ploidy, min_purity=min_purity, max_purity=max_purity, 
                    min_aligned_seg_mb=min_aligned_seg_mb, max_homdel_mb=max_homdel_mb, 
                    delta_tcn_to_int=delta_tcn_to_int, delta_tcn_to_avg=delta_tcn_to_avg, delta_tcnavg_to_int=delta_tcnavg_to_int, 
                    delta_mcn_to_int=delta_mcn_to_int, delta_mcn_to_avg=delta_mcn_to_avg, delta_mcnavg_to_int=delta_mcnavg_to_int, 
                    mcn_weight=mcn_weight, rho=rho, timeout=timeout, min_cna_segments_per_sample=min_cna_segments_per_sample, 
                    obj2_clonalonly=obj2_clonalonly)

    #def do_CNAlign(dat, gurobi_license, min_ploidy=1.7, max_ploidy=6.5, min_purity=0.05, max_purity=0.95, min_aligned_seg_mb=5.0, max_homdel_mb=50.0, delta_tcn_to_int=0.2, delta_tcn_to_avg=0.1, delta_tcnavg_to_int=0.05, delta_mcn_to_int=0.2, delta_mcn_to_avg=0.1, delta_mcnavg_to_int=0.05, mcn_weight=0.5, rho=1.0, timeout=10*60, min_cna_segments_per_sample=3, obj2_clonalonly=0):

	end_time <- now()
	run_date <- as.character(format(end_time,format='%Y-%m-%d %H:%M'))
    message('Ended CNAlign at ',as.character(end_time))
    sec_elapsed <- round(as.numeric(end_time) - as.numeric(start_time))
    message('Time elapsed: ',sec_elapsed,'s')

    # return the model object and a list of all the parameters
    out
}


