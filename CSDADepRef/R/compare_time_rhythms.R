#' Plot and compare EEG brain rhythms from two time series
#' 
#' @param spectrum_smoothness      Cut-frequency (0-1) of a 3rd Butterworth LPF
#'
#'
#' @keywords EEG
#' @export
#' @examples 
#' \dontrun{
#'  compare_time_rhythms()
#' }
compare_time_rhythms <- function(x=NULL, y=NULL, target_rhythm="alpha", fs=100, N=600, flim=c(0, 0.5*fs), channel_label1="p", channel_label2="q", spectrum_smoothness=0.4) {
  if(is.null(x)){
    x <- simulate_EEG(N, fs=fs)
  }else{
    x <- x[1: N]
  }
  if(is.null(y)){
    y <- simulate_EEG(N, fs=fs)
  }else{
    y <- y[1: N]
  }
  include_filter <- !is.null(target_rhythm)
  
  plot_labels <- list(
    f="Frequency [Hz]",
    t="Time [s]",
    S="Spectrum [AU]",
    x="Amplitude [AU]",
    fx=paste0("Amplitude [AU] ", if(!is.null(target_rhythm)) sprintf("(%s band)", target_rhythm) else "")
  )
  label_id <- function(i) intToUtf8(utf8ToInt("A") + i - 1)
  
  Sx0 <- estimate_smooth_spectrum(x, spectrum_smoothness)
  Sy0 <- estimate_smooth_spectrum(y, spectrum_smoothness)
  t <- (1:length(x)) / fs
  f <- (1:length(Sx0)) / length(Sx0) * 0.5 * fs
  f_indices <- (f >= flim[1]) & (f <= flim[2])
  f <- f[f_indices]
  Sx0 <- Sx0[f_indices]
  Sy0 <- Sy0[f_indices]
  Sx <- NULL
  Sy <- NULL
  
  if(include_filter){
    r.x <- extract_rhythms(x, fs, rhythms=c(target_rhythm))
    r.y <- extract_rhythms(y, fs, rhythms=c(target_rhythm))
    Sx <- estimate_smooth_spectrum(r.x[[target_rhythm]], 0.4)
    Sy <- estimate_smooth_spectrum(r.y[[target_rhythm]], 0.4)
    Sx <- Sx[f_indices]
    Sy <- Sy[f_indices]
  }
  
  color1 <- "#E69F00"
  color2 <- "#56B4E9"
  theme <- csda_theme()

  # https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
  #psx <- ggplot2::qplot(f, Sx, geom="line", main=TeX("$S_p$"), xlab="", ylab="", xlim=c(0, 0.5 *fs)) + theme
  #psy <- ggplot2::qplot(f, Sy, geom="line", main=TeX("$S_q$"), xlab=plot_labels$f, ylab=plot_labels$S, xlim=c(0, 0.5 *fs)) + theme
  ncols <- if(include_filter) 3 else 2
  
  psx <- dep_plot_spectrum(fs=fs, Sx0=Sx0, Sx=Sx, flim=flim, color=color1, label_id=1, channel_label=channel_label1)
  psy <- dep_plot_spectrum(fs=fs, Sx0=Sy0, Sx=Sy, flim=flim, color=color2, label_id=2, channel_label=channel_label2)
  
  
  px <- ggplot2::qplot(t, x, color=I(color1), geom="line", main=TeX(paste0(label_id(3), ") ", "$X_{", channel_label1, "}$")), xlab="", ylab="") + theme
  py <- ggplot2::qplot(t, y, color=I(color2), geom="line", main=TeX(paste0(label_id(4), ") ", "$X_{", channel_label2, "}$")), xlab=plot_labels$t, ylab=plot_labels$x) + theme
  
  if(include_filter) {
    pax <- ggplot2::qplot(t, r.x[[target_rhythm]], color=I(color1), geom="line", main=TeX(paste0(label_id(5), ") ", "$X_{\\", target_rhythm, ",", channel_label1, "}$")), xlab="", ylab="") + theme
    pay <- ggplot2::qplot(t, r.y[[target_rhythm]], color=I(color2), geom="line", main=TeX(paste0(label_id(6), ") ", "$X_{\\", target_rhythm, ",", channel_label2, "}$")), xlab=plot_labels$t, ylab=plot_labels$fx) + theme
    
    gridExtra::grid.arrange(psx, px, pax, psy, py, pay, nrow=2, widths=c(2, 3, 3))
  } else {
    gridExtra::grid.arrange(psx, px, psy, py, nrow=2, widths=c(2, 3))
  }
}
