R0_scatterplot <- function(pred_R01, pred_R02, filename, xLabel, yLabel, labelPoints){
  filepath <- paste0('figures/', filename, '.pdf')
  maxLim <- round(max(unlist(pred_R01), unlist(pred_R02), na.rm = T))
  pdf(filepath, height = 6, width = 8)
  plot(
    unlist(pred_R01)
    , unlist(pred_R02)
    , pch = 16
    , xlab = xLabel
    , ylab = yLabel
    , xlim = c(0, maxLim)
    , ylim = c(0, maxLim)
  )
  abline(h = 1, lty = 2)
  abline(v = 1, lty = 2)
  if(labelPoints == TRUE){
    pts_to_label <- which(unlist(pred_R02) > 1)
    text(pred_R01[pts_to_label], pred_R02[pts_to_label]-0.08, labels = survey_data$location_code[pts_to_label])
  }
  dev.off()
  
}