#' Sensitivity analysis
#'
#' Rank core instruments by magnitude of effect (beta)
#'
#' @param harmonized_data output of TwoSampleMR::harmonise_data()
#' @param K number of instrument groups, default = 5
#' @param method_list MR methods, selected from TwoSampleMR::mr_method_list()
#' 
#' @return list containing \itemize{
#' \item \code{result} results of MR sensitivity analysis
#' \item \code{instrument} summary of instruments and group membership
#' \item \code{Q} Cochran’s Q statistic
#' \item \code{I2} Higgins & Thompson’s I2 statistic
#' \item \code{het_P} heterogeneity p-value
#' }
#' 
#' @export
#' 

mr_sens_beta <- function(harmonized_data, K = 5, method_list = c("mr_ivw","mr_simple_mode","mr_weighted_mode","mr_weighted_median","mr_egger_regression")) {
  harmonized_data = harmonized_data[harmonized_data$mr_keep,]
  harmonized_data = harmonized_data[!is.na(harmonized_data$se.exposure),]
  harmonized_data$instrument_group = (K + 1 - as.numeric(cut_number(abs(harmonized_data$beta.exposure),K)))
  if ("chr.exposure" %in% colnames(harmonized_data) & "pos.exposure" %in% colnames(harmonized_data)) {
    instrument_group = data.frame(instrument = harmonized_data$SNP,
                                  beta.exposure = harmonized_data$beta.exposure,
                                  instrument_group = harmonized_data$instrument_group,
                                  chr = harmonized_data$chr.exposure,
                                  pos = harmonized_data$pos.exposure)
  } else {
    instrument_group = data.frame(instrument = harmonized_data$SNP,
                                  beta.exposure = harmonized_data$beta.exposure,
                                  instrument_group = harmonized_data$instrument_group)
  }
  all_result = data.frame(id.exposure = character(),
                          id.outcome = character(),
                          outcome = character(),
                          exposure = character(),
                          method = character(),
                          nsnp = numeric(),
                          b = numeric(),
                          se = numeric(),
                          pval = numeric(),
                          Fstat = numeric(),
                          group = numeric(),
                          range = character())
  for (group in 1:K) {
    cat("obtaining discrete estimates for instrument group",group,"/",K,"...",collapse="\n")
    result_group = mr(harmonized_data[harmonized_data$instrument_group==group,],method_list = method_list)
    if ("mr_egger_regression" %in% method_list) {
      result_group_egger_intercept = mr_egger_regression(harmonized_data[harmonized_data$instrument_group==group,]$beta.exposure,
                                                         harmonized_data[harmonized_data$instrument_group==group,]$beta.outcome,
                                                         harmonized_data[harmonized_data$instrument_group==group,]$se.exposure,
                                                         harmonized_data[harmonized_data$instrument_group==group,]$se.outcome)
      result_group = rbind.data.frame(result_group,
                                      data.frame(id.exposure = result_group$id.exposure[1],
                                                 id.outcome = result_group$id.outcome[1],
                                                 outcome = result_group$outcome[1],
                                                 exposure = result_group$exposure[1],
                                                 method = "MR Egger intercept",
                                                 nsnp = result_group_egger_intercept$nsnp,
                                                 b = result_group_egger_intercept$b_i,
                                                 se = result_group_egger_intercept$se_i,
                                                 pval = result_group_egger_intercept$pval_i))
    }
    R2 = sum((harmonized_data[harmonized_data$instrument_group==group,]$beta.exposure / harmonized_data[harmonized_data$instrument_group==group,]$se.exposure)**2 / ((harmonized_data[harmonized_data$instrument_group==group,]$beta.exposure / harmonized_data[harmonized_data$instrument_group==group,]$se.exposure)**2 + harmonized_data[harmonized_data$instrument_group==group,]$samplesize.exposure - 2))
    result_group$Fstat = (mean(harmonized_data[harmonized_data$instrument_group==group,]$samplesize.exposure) - nrow(harmonized_data[harmonized_data$instrument_group==group,]) - 1) / nrow(harmonized_data[harmonized_data$instrument_group==group,]) * R2 / (1 - R2)
    result_group$group = group
    result_group$range = "Discrete"
    all_result = rbind.data.frame(all_result, result_group)
  }
  for (group in 1:K) {
    cat("obtaining cumulative estimates for instrument group",group,"/",K,"...",collapse="\n")
    result_group = mr(harmonized_data[harmonized_data$instrument_group<=group,],method_list = method_list)
    if ("mr_egger_regression" %in% method_list) {
      result_group_egger_intercept = mr_egger_regression(harmonized_data[harmonized_data$instrument_group<=group,]$beta.exposure,
                                                         harmonized_data[harmonized_data$instrument_group<=group,]$beta.outcome,
                                                         harmonized_data[harmonized_data$instrument_group<=group,]$se.exposure,
                                                         harmonized_data[harmonized_data$instrument_group<=group,]$se.outcome)
      result_group = rbind.data.frame(result_group,
                                      data.frame(id.exposure = result_group$id.exposure[1],
                                                 id.outcome = result_group$id.outcome[1],
                                                 outcome = result_group$outcome[1],
                                                 exposure = result_group$exposure[1],
                                                 method = "MR Egger intercept",
                                                 nsnp = result_group_egger_intercept$nsnp,
                                                 b = result_group_egger_intercept$b_i,
                                                 se = result_group_egger_intercept$se_i,
                                                 pval = result_group_egger_intercept$pval_i))
    }
    R2 = sum((harmonized_data[harmonized_data$instrument_group<=group,]$beta.exposure / harmonized_data[harmonized_data$instrument_group<=group,]$se.exposure)**2 / ((harmonized_data[harmonized_data$instrument_group<=group,]$beta.exposure / harmonized_data[harmonized_data$instrument_group<=group,]$se.exposure)**2 + harmonized_data[harmonized_data$instrument_group<=group,]$samplesize.exposure - 2))
    result_group$Fstat = (mean(harmonized_data[harmonized_data$instrument_group<=group,]$samplesize.exposure) - nrow(harmonized_data[harmonized_data$instrument_group<=group,]) - 1) / nrow(harmonized_data[harmonized_data$instrument_group<=group,]) * R2 / (1 - R2)
    result_group$group = group
    result_group$range = "Cumulative"
    all_result = rbind.data.frame(all_result, result_group)
  }
  Q = sapply(unique(all_result$method[all_result$method!="MR Egger intercept"]), function(m) Qstatistic(all_result[all_result$method==m & all_result$range=="Discrete",]$b,
                                                                                                        all_result[all_result$method==m & all_result$range=="Discrete",]$se,
                                                                                                        all_result[all_result$method==m & all_result$range=="Cumulative" & all_result$group==K,]$b))
  I2 = (Q - (K - 1)) / Q
  I2 = ifelse(I2 > 0, I2, 0)
  het_P = pchisq(Q, df = K - 1, lower.tail = F)
  return(list(result = all_result,
              instrument = instrument_group,
              Q = Q,
              I2 = I2,
              het_P = het_P))
}
