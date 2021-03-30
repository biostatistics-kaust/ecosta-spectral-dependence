separate_low_high_freq <- function(x, freq_low, freq_high, fs) {
  coeffs_low <- signal::fir1(100, pmax(1e-4, freq_low/(0.5 * fs)), type="pass")
  coeffs_high <- signal::fir1(100,  pmax(1e-4, freq_high/(0.5 * fs)), type="pass")
  L <- length(x)
  y_l <- x
  y_h <- x
  for(i in 1:1){
    y_l <- as.numeric(signal::filtfilt(coeffs_low, c(y_l, y_l))[(L+1):(2 * L)])
    y_h <- as.numeric(signal::filtfilt(coeffs_high, c(y_h, y_h))[(L+1):(2 * L)])
  }
  y_l <- y_l - mean(y_l)
  y_h <- y_h - mean(y_h)
  list(low=y_l, high=y_h, fs=fs)
}

modulation_index <- function(phase, amplitude) {
  N <- length(phase)
  amplitude <- amplitude / sum(amplitude)
  reference <- 1 / N
  sum(amplitude * (log10(amplitude) - log10(reference))) / log10(N)
}

plot_modulation_index <- function(data, levels=20, low_wrt_high=TRUE, label_id="") {
  if(low_wrt_high) {
    phase_h <- Arg(seewave::hilbert(data$low, data$fs))
    amp_h <- abs(seewave::hilbert(data$high, data$fs))
  } else {
    phase_h <- Arg(seewave::hilbert(data$high, data$fs))
    amp_h <- abs(seewave::hilbert(data$low, data$fs))
  }
  if(is.numeric(label_id)) {
    label_id <- intToUtf8(utf8ToInt("A") + label_id - 1)
  }

  breaks <- seq(-pi, pi, length.out=levels)
  lev_phases <- cut(phase_h, breaks=breaks)
  amplitude <- sapply(levels(lev_phases), function(x) mean(amp_h[x == lev_phases]))
  phase <- 0.5 * (breaks[2:length(breaks)] + breaks[1:length(breaks)-1])
  
  title <- sprintf("%s) Amplitude-phase plot. MI: %.2fe-3",
    label_id,
    modulation_index(phase, amplitude) * 1e3
  )
  
  amplitude_noise <- mean(amplitude)#1 / levels
  pseudo_hist <- data.frame(phase=phase, amplitude=amplitude)
  #data_collect <- data.frame(phase0=phase_h, amplitude0=amp_h )
  q <- ggplot(pseudo_hist) +
  #geom_line(aes(x=phase, y=amplitude)) + 
  #geom_ribbon(aes(x=phase, ymin=0, ymax=amplitude), color="black", alpha=0.1) + 
  geom_col(aes(x=phase, y=amplitude)) + 
  #geom_point(aes(x=phase0, y=amplitude0), alpha=0.05, data=data_collect) + 
  geom_line(aes(x=phase, y=amplitude_noise), color="blue", linetype="dashed") + 
  labs(title=title) +
  csda_theme()
  q
}
