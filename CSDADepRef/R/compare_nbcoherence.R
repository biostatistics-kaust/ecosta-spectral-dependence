#' Plot and compare narrow-band coherence between two channels
#' 
#' @keywords EEG
#' @export
#' @examples 
#' \dontrun{
#'  compare_nbcoherence(max_lags_plot=200, max_lags_calculate=200)
#' }

compare_nbcoherence <- function(x=NULL, y=NULL, target_rhythm="alpha", max_lags_plot=40, max_lags_calculate=40, max_plots=3, fs=100, N=600) {
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
  
  r.x <- extract_rhythms(x, fs, rhythms=c(target_rhythm))
  r.y <- extract_rhythms(y, fs, rhythms=c(target_rhythm))
  t <- (1:length(x)) / fs
  
  Sx <- estimate_smooth_spectrum(x)
  Sy <- estimate_smooth_spectrum(y)
  f <- (1:length(Sx)) / length(Sx) * 0.5 * fs
  
  plot_labels <- list(
    f="Frequency [Hz]",
    t="Time [s]",
    S="Spectrum [AU]",
    x="Amplitude [AU]",
    fx=paste0("Filtered amplitude [AU] (", target_rhythm, " band)")
  )
  label_id <- function(i) intToUtf8(utf8ToInt("A") + i)
  n <- 0
  theme <- csda_theme()
  
  lags_to_plot <- floor(seq(from=0, to=max_lags_plot, length=max_plots+1))[1 + (1:max_plots)]
  
  #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  color1 <- "#E69F00"
  color2 <- "#56B4E9"
  color_red <- "#db5067"
  color_blue <- "#2b76b2"
  color_gold <- "#e4a027"
  
  base_plot <- ggplot2::qplot(
    x, y,
    geom="point",
    main=TeX(paste0(label_id(n), ") $r_{0} = ", round(cor(
      x, y
    ) ^ 2, 3), "$")),
    xlab=TeX(paste0("$X_{\\p}$")),
    ylab=TeX(paste0("$X_{\\q}$")),
    alpha=I(0.1),
    color=I("black"),
  ) + theme
  #n <- n + 1
    
  # https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
  lag_plots_neg <- lapply(c(rev(lags_to_plot), 0), function(lag)
    {
      n <<- n + 1
      ggplot2::qplot(
        r.x[[target_rhythm]][1:(length(y)-lag)], r.y[[target_rhythm]][(lag+1):length(x)],
        geom="point",
        #main=TeX(paste0("$l=", -lag, "$")),
        main=TeX(paste0(
          label_id(n), ") ",
            #"Lag ", -lag, ": ",
          "$r_{p,q}(", -lag, ") = ", round(cor(
          r.x[[target_rhythm]][1:(length(y)-lag)], r.y[[target_rhythm]][(lag+1):length(x)]
        ) ^ 2, 3), "$")),
        xlab=TeX(paste0("$X_{\\p,", target_rhythm, "}$")),
        ylab=if(lag==max_lags_plot) TeX(paste0("$X_{\\q,", target_rhythm, "}$")) else "",
        alpha=I(0.1),
        color=I(if(lag!=0) color_red else color_gold),
      ) + theme
    }
  )
  lag_plots_pos <- lapply(lags_to_plot, function(lag)
    {
      n <<- n + 1
      ggplot2::qplot(
        r.x[[target_rhythm]][(lag+1):length(x)], r.y[[target_rhythm]][1:(length(y)-lag)], 
        geom="point",
        #main=TeX(paste0("$l=", lag, "$")),
        main=TeX(paste0(
          label_id(n), ") ",
            #"Lag ", lag, ": ",
          "$r_{p,q}(", lag, ") = ", round(cor(
          r.x[[target_rhythm]][(lag+1):length(x)], r.y[[target_rhythm]][1:(length(y)-lag)]
        ) ^ 2, 3), "$")),
        xlab=TeX(paste0("$X_{\\", target_rhythm, ",p}$")),
        ylab="", #""TeX(paste0("$X_{\\", target_rhythm, ",q}$")),
        alpha=I(0.1),
        color=I(color_blue)
      ) + theme
    }
  )

  corr_plot <- dep_coherence_lagged_plot(
    max_lags=max_lags_calculate,
    x=r.x[[target_rhythm]],
    y=r.y[[target_rhythm]],
    corr_var_name="r_{p,q}",
    x_var_name=paste0("X_{\\p,", target_rhythm, "}"),
    y_var_name=paste0("X_{\\q,", target_rhythm, "}")
  )

  M <- length(lag_plots_pos) * 2 + 2
  lag_plots <- c(list(base_plot), lag_plots_neg, lag_plots_pos, list(corr_plot))
  

  gridExtra::grid.arrange(grobs=lag_plots, nrow=2, layout_matrix=rbind(1:M, rep(M+2, M)))
}
