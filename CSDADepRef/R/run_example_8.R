#' @export
run_example_8 <- function() {
  band_pass_filter_helper <- function(x, fs, freqs=c(1, 3), causal=TRUE){
    #
    coeffs_low <- signal::butter(4, freqs/(0.5 * fs), type="pass")
    coeffs_low <- signal::fir1(20, freqs/(0.5 * fs), type="pass")
    coeffs_low <- signal::fir1(100, freqs/(0.5 * fs), type="pass")
    if(causal)
      y_c <- as.numeric(signal::filter(coeffs_low, x))
    else
      y_c <- as.numeric(signal::filtfilt(coeffs_low, x))
    y_c <- y_c - mean(y_c)
    y_c
  }
  set.seed(1008)
  fs <- 128*10
  N = fs * 10
  #N <- 200
  mN = 500
  t <- (1:N) / fs
  x1 <- simulate_unstable_oscillator(N=N+mN, f=2, tau=4, fs=fs, sigma=1)
  x2 <- simulate_unstable_oscillator(N=N+mN, f=15, tau=4, fs=fs, sigma=1)
  x2 <- x2 + simulate_unstable_oscillator(N=N+mN, f=30, tau=4, fs=fs, sigma=1)
  
  lag = 10
  x <- rnorm(N) * 0.5 + x1[1:N+lag]
  y <- rnorm(N) * 0.5 + x1[1:N+0] + 1 * x2[1:N+0]
  
  L <- length(x)
  
  Tf <- (fs * 2)
  T0 <- (fs * 5)
  y_c <- band_pass_filter_helper(x, fs, c(0.1,4), causal=TRUE)[(1:Tf)+T0]
  y_nc <- band_pass_filter_helper(x, fs, c(0.1,4), causal=FALSE)[(1:Tf)+T0]
  y_c0 <- band_pass_filter_helper(y, fs, c(0.1,4), causal=TRUE)[(1:Tf)+T0]
  y_nc0 <- band_pass_filter_helper(y, fs, c(0.1,4), causal=FALSE)[(1:Tf)+T0]
  
  df_plot <- data.frame(
    t=t[(1:Tf)], 
    #y_causal=y_c, 
    #y_noncausal=y_nc, 
    y_causal=y_c0, 
    y_noncausal=y_nc0,
    y_causal_ref=y[(1:Tf)+T0], 
    y_noncausal_ref=y[(1:Tf)+T0]
    #y_causal_ref=x[(1:Tf)+T0], 
    #y_noncausal_ref=x[(1:Tf)+T0]
  )
  
  q1 <- ggplot(df_plot) +
    geom_line(aes(x=t, y=y_causal_ref), color="gray", alpha=0.4) +
    geom_line(aes(x=t, y=y_causal), color="darkred") +
    csda_theme() +
    labs(title="A) One-sided FIR(100) filtered signal", x="time [sec]", y="amplitude")
    
  q2 <- ggplot(df_plot) +
    geom_line(aes(x=t, y=y_noncausal_ref), color="gray", alpha=0.4) +
    geom_line(aes(x=t, y=y_noncausal), color="darkblue") +
    csda_theme() +
    labs(title="B) Two-sided FIR(100) filtered signal", x="time [sec]", y="amplitude")
  
  corr_helper <- function(a, b) {
    corr_info <- coherence_corr(a, b, max_lags=fs, sign_lags="-")
    data.frame(
      lags=corr_info$correlation_indices,
      xcorr=corr_info$correlation_values
    )
  }
  
  corr1 <- corr_helper(y_c, y_c0)
  corr_max <- max(corr1$xcorr)
  corr_lag_max <- corr1$lags[which.max(corr1$xcorr)]
  q_corr1 <- ggplot(corr1) +
    geom_line(aes(x=lags, y=xcorr), color="darkblue") +
    geom_vline(xintercept=corr_lag_max, color="darkblue", linetype="dashed") +
    csda_theme() +
    coord_cartesian(ylim=c(0, 1)) +
    labs(y="Cross-correlation", title=TeX(sprintf(
      "C) Cross-corr. $%s$ and $%s$. $\\\\textbf{C}$oh=%.3f (lag=%d)",
      #"\\textbf{C}^{causal}_{delta}\\[  Z_{delta}(t)  \\]",
      #"\\textbf{C}^{causal}_{delta}\\[  Z_{delta}(t-10) + Z_{alpha}(t)  \\]",
      #"\\textbf{C}^{causal}_{delta}\\[  X_1(t)  \\]",
      #"\\textbf{C}^{causal}_{delta}\\[  X_2(t)  \\]",
      "X^{(1-sided)}_{1, delta}(t)",
      "X^{(1-sided)}_{2, delta}(t)",
      corr_max,
      corr_lag_max
    )))
  
  corr2 <- corr_helper(y_nc, y_nc0)
  corr_max <- max(corr2$xcorr)
  corr_lag_max <- corr2$lags[which.max(corr2$xcorr)]
  q_corr2 <- ggplot(corr2) +
    geom_line(aes(x=lags, y=xcorr), color="darkred") +
    geom_vline(xintercept=corr_lag_max, color="darkred", linetype="dashed") +
    csda_theme() +
    coord_cartesian(ylim=c(0, 1)) +
    labs(y="Cross-correlation", title=TeX(sprintf(
      "D) Cross-corr. $%s$ and $%s$. $\\\\textbf{C}$oh=%.3f (lag=%d)",
      #"\\textbf{C}^{non causal}_{delta}\\[  X_1(t)  \\]",
      #"\\textbf{C}^{non causal}_{delta}\\[  X_2(t)  \\]",
      "X^{(2-sided)}_{1, delta}(t)",
      "X^{(2-sided)}_{2, delta}(t)",
      corr_max,
      corr_lag_max
    )))
    
  q <- gridExtra::grid.arrange(grobs=list(
      q1,
      q2,
      
      q_corr1,
      q_corr2
      
    ), ncol=2,
    layout_matrix=rbind(
      rep(c(1, 2), each=1),
      rep(c(3, 4), each=1)
    )
  )
  q
}

