#' @export
run_example_6 <- function () {
  set.seed(1006)
  fs <- 128*4*2
  N <- fs * 2.5
  tau <- 0.5
  Z_theta <- simulate_unstable_oscillator(N=N, f=5, tau=6, fs=fs, sigma=1)
  Z_gamma <- simulate_unstable_oscillator(N=N, f=32, tau=7, fs=fs, sigma=1)
  freq_low <- c(0, 30)
  freq_high <- c(30, 45)
  am_modulation <- 1*(Z_gamma+1.0) * (0 + Z_theta) + 2 * Z_gamma + rnorm(N)*0.1
  extreme_mod <- 4 * Z_theta + (Z_gamma) * (Z_theta) + rnorm(N)*0.1
  t <- (1:length(Z_theta)) / fs
  q_src1 <- qplot(t, Z_theta, geom="line") +
    labs(title=TeX("A) Latent source $Z_{delta}$"), x="time [sec]", y="amplitude") +
    csda_theme()
  q_mod <- qplot(t, am_modulation, geom="line") +
    labs(title=TeX("G) Component $X_1$"), x="time [sec]", y="amplitude") +
    csda_theme()
  
  q_src2 <- qplot(t, Z_gamma, geom="line") +
    labs(title=TeX("B) Latent source $Z_{gamma}$"), x="time [sec]", y="amplitude") +
    csda_theme()
  q_extr <- qplot(t, extreme_mod, geom="line") +
    labs(title=TeX("H) Component $X_2$"), x="time [sec]", y="amplitude") +
    csda_theme()
  
  snoise <- rnorm(N)
  q_snoise <- qplot(t, snoise, geom="line") +
    labs(title=TeX("C) White noise"), x="time [sec]", y="amplitude") +
    csda_theme()
  
  q_split_src1 <- plot_modulation_index(separate_low_high_freq(Z_theta, freq_low=freq_low, freq_high=freq_high, fs=fs), label_id=4)
  q_split_src2 <- plot_modulation_index(separate_low_high_freq(Z_gamma, freq_low=freq_low, freq_high=freq_high, fs=fs), label_id=5)
  q_split_snoise <- plot_modulation_index(separate_low_high_freq(snoise, freq_low=freq_low, freq_high=freq_high, fs=fs), label_id=6)
  
  q_split_mod <- plot_modulation_index(separate_low_high_freq(am_modulation, freq_low=freq_low, freq_high=freq_high, fs=fs), label_id=9)
  q_split_extr <- plot_modulation_index(separate_low_high_freq(extreme_mod, freq_low=freq_low, freq_high=freq_high, fs=fs), label_id=10)
    
  gridExtra::grid.arrange(grobs=list(
      q_src1,
      q_src2,
      
      q_split_src1,
      q_split_src2,
      
      q_mod,
      q_extr,
      
      q_split_mod,
      q_split_extr,
      
      
      q_snoise,
      q_split_snoise
      
    ), ncol=6,
    layout_matrix=rbind(
      rep(c(1, 2, 9), each=2),
      rep(c(3, 4, 10), each=2),
      rep(c(5, 6), each=3),
      rep(c(7, 8), each=3)
    )
  )# + csda_theme()
}


