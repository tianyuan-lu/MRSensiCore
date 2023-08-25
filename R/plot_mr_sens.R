#' Plot results of sensitivity analysis
#'
#' Illustrate MR estimates by instrument group
#'
#' @param mr_sens_result output of \link{mr_sens_h2}, \link{mr_sens_z}, or \link{mr_sens_beta}
#' @param scale "linear": linear scale, or "exp": exponential scale
#' 
#' @return ggplot object
#' 
#' @export
#' 

plot_mr_sens <- function(mr_sens_result, scale = "linear") {
  plotdat = mr_sens_result$result
  plotdat$group = factor(plotdat$group, levels = c(1:max(plotdat$group)))
  if (scale == "linear") {
    ggplot(plotdat[plotdat$method!="MR Egger intercept",], aes(x = group, y = b, color = range)) +
      geom_hline(yintercept = 0, lty = 2, col = "darkgrey") +
      facet_wrap(method~.) +
      geom_point(position = position_dodge(width = 0.3)) +
      geom_errorbar(aes(ymin = b - 1.96 * se,
                        ymax = b + 1.96 * se),width = 0,
                    position = position_dodge(width = 0.3)) +
      theme_bw() +
      xlab("Instrument group (core instrument >>> --- >>> peripheral instrument)") +
      ylab(expression(hat(beta)[MR])) +
      ggtitle(paste0("Exposure: ",plotdat$exposure[1],"\n","Outcome: ",plotdat$outcome[1])) +
      labs(color = "") +
      scale_color_manual(values = c("red","royalblue")) +
      theme(axis.title = element_text(size = 15),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.position = "bottom") -> plt
  }
  if (scale == "exp") {
    ggplot(plotdat[plotdat$method!="MR Egger intercept",], aes(x = group, y = exp(b), color = range)) +
      geom_hline(yintercept = 1, lty = 2, col = "darkgrey") +
      facet_wrap(method~.) +
      geom_point(position = position_dodge(width = 0.3)) +
      geom_errorbar(aes(ymin = exp(b - 1.96 * se),
                        ymax = exp(b + 1.96 * se)),width = 0,
                    position = position_dodge(width = 0.3)) +
      theme_bw() +
      xlab("Instrument group (core instrument >>> --- >>> peripheral instrument)") +
      ylab(expression(hat(OR)[MR])) +
      ggtitle(paste0("Exposure: ",plotdat$exposure[1],"\n","Outcome: ",plotdat$outcome[1])) +
      labs(color = "") +
      scale_color_manual(values = c("red","royalblue")) +
      theme(axis.title = element_text(size = 15),
            axis.text = element_text(size = 14),
            strip.text = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.position = "bottom") -> plt
  }
  return(plt)
}
