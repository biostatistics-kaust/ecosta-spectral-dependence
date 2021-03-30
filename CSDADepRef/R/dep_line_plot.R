# label => if label == 1: label="A"
#          else: label is used
#' @export
dep_line_plot <- function(x, fs, label_id, xlabel="Relative time \\[sec\\]", ylabel="Amplitude \\[AU\\]", title="", color="black", alpha=0.9) {
  if(as.numeric(label_id)) {
    label_id <- intToUtf8(utf8ToInt("A") + label_id - 1)
  }
  theme <- csda_theme()
  t <- (1: length(x)) / fs
  base_plot <- ggplot2::qplot(
    t, x,
    geom="line",
    main=TeX(paste0(label_id, ") ", title)),
    xlab=TeX(xlabel),
    ylab=TeX(ylabel),
    alpha=I(alpha),
    color=I(color),
  ) + theme
}

