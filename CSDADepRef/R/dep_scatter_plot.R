# label => if label == 1: label="A"
#          else: label is used
#' @export
dep_scatter_plot <- function(x, y, label_id, xlabel, ylabel, color, corr_var_name="r_{0}", alpha=0.5) {
  if(as.numeric(label_id)) {
    label_id <- intToUtf8(utf8ToInt("A") + label_id - 1)
  }
  #theme <- theme_linedraw()
  theme <- csda_theme()
  #lags_to_plot <- floor(seq(from=0, to=max_lags_plot, length=max_plots+1))[1 + (1:max_plots)]
  base_plot <- ggplot2::qplot(
    x, y,
    geom="point",
    main=TeX(paste0(label_id, ") $", corr_var_name, " = ", sprintf("%.4g", cor(
      x, y
    ) ^ 2 ), "$")),
    xlab=TeX(xlabel),
    ylab=TeX(ylabel),
    alpha=I(alpha),
    color=I(color),
  ) + theme
}

