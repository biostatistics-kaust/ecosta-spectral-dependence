#' @export
dep_plot_spectrum <- function(fs, Sx0, Sx=NULL, flim=c(0, 0.5*fs), color="red", label_id=1, channel_label="CH") {
  include_filter <- !is.null(Sx)
  f <- ((1:length(Sx0)) - 1) / length(Sx0) * (0.5 * fs)
  if(as.numeric(label_id)) {
    label_id <- intToUtf8(utf8ToInt("A") + label_id - 1)
  }
  theme <- csda_theme()
  
  ggdata <- data.frame(f=f, Sx0=Sx0)
  maxSx <- max(Sx0) * 1.3
  if(include_filter){
    ggdata$Sx <- Sx
    maxSx <- max(Sx) * 1.3
  }
  psx <- ggplot(ggdata) +
    geom_line(aes(x=f, y=Sx0), color=(if(include_filter) "black" else color), alpha=0.3) +
    (if(include_filter) geom_line(aes(x=f, y=Sx), color=color) else labs(title="")) +
    geom_ribbon(aes(x=f, ymin=0, ymax=(if(include_filter) Sx else Sx0 )), fill=color, alpha=0.3) +  #IDK
    labs(title=TeX(sprintf(
      "%s) $C_{%s}}(\\omega)$",
      label_id,
      channel_label
    )), x="", y="") +
    #coord_cartesian(xlim=flim, ylim=c(0, maxSx)) +
    theme
}

