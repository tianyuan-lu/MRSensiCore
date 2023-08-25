#' Permutation test
#'
#' Obtain distribution of MR estimates by randomly sampling instruments
#'
#' @param harmonized_data output of TwoSampleMR::harmonise_data()
#' @param mr_sens_result output of \link{mr_sens_h2}, \link{mr_sens_z}, or \link{mr_sens_beta}
#' @param core_group vector containing groups that will be considered core instrument groups, default = 1
#' @param exclude_group vector containing groups that will be excluded when sampling instruments, default = 1
#' @param Nperm number of permutations
#' @param method_list MR methods, selected from TwoSampleMR::mr_method_list()
#' 
#' @return data frame containing all MR estimates from permutations
#' 
#' @export
#' 

permute_core_instrument <- function(harmonized_data, mr_sens_result, core_group = 1, exclude_group = 1, Nperm = 1000, method_list = c("mr_ivw","mr_weighted_mode","mr_weighted_median","mr_penalised_weighted_median","mr_egger_regression")) {
  harmonized_data = harmonized_data[harmonized_data$mr_keep,]
  harmonized_data = harmonized_data[!is.na(harmonized_data$se.exposure),]
  if ("eaf.exposure" %in% colnames(mr_sens_result)) {
    harmonized_data = harmonized_data[!is.na(harmonized_data$eaf.exposure),]
  }
  harmonized_data = harmonized_data[!harmonized_data$SNP %in% mr_sens_result$instrument$instrument[mr_sens_result$instrument$instrument_group %in% exclude_group],]
  Nins = sum(mr_sens_result$instrument$instrument_group %in% core_group)
  all_sample_result = data.frame(id.exposure = character(),
                                 id.outcome = character(),
                                 outcome = character(),
                                 exposure = character(),
                                 method = character(),
                                 nsnp = numeric(),
                                 b = numeric(),
                                 se = numeric(),
                                 pval = numeric())
  for (ite in 1:Nperm) {
    set.seed(ite)
    sample_harmonized_data = harmonized_data[sample(1:nrow(harmonized_data), Nins),]
    sample_result = mr(sample_harmonized_data,method_list = method_list)
    all_sample_result = rbind.data.frame(all_sample_result, sample_result)
  }
  return(all_sample_result)
}
