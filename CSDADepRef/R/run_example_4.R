#' @export
run_example_4 <- function(){
  set.seed(1004)
  fs <- 50
  x1 <- simulate_unstable_oscillator(N=300, f=10, tau=-log(0.05), fs=fs, sigma=1)
  x2 <- simulate_unstable_oscillator(N=300, f=10, tau=-log(0.50), fs=fs, sigma=1)  
  compare_time_rhythms(x1, x2, N=length(x1), target_rhythm=NULL, flim=c(0, 0.5*fs), channel_label1="1", channel_label2="2", fs=fs, spectrum_smoothness=0.9)
}

#' @export
run_example_4b <- function(){
  set.seed(1004)
  fs <- 50
  f <- 10
  tau1 <- -log(0.05)
  tau2 <- -log(0.50)
  #tau1 <- -log(0.00001)
  r1 <- simulate_unstable_oscillator(N=NULL, f=f, tau=tau1, fs=fs, sigma=1)
  r2 <- simulate_unstable_oscillator(N=NULL, f=f, tau=tau2, fs=fs, sigma=1)
  S1 <- estimate_true_ar2_spectrum(M=1000, fs=fs, r1$phi1, r1$phi2)
  S2 <- estimate_true_ar2_spectrum(M=1000, fs=fs, r2$phi1, r2$phi2)
  
  ggdata <- data.frame(
    S1=S1$S,
    S2=S2$S,
    f=S1$f
  )
  color1 <- "#db5067"
  color2 <- "#56B4E9"
  #cat(sprintf("$\\phi_1=\\frac{2}{%.3g}cos(2\\pi\\frac{%.3g}{%.3g})$", 1+exp(tau1), f, fs))
  eq1 <- (sprintf("$\\phi_1=\\frac{2}{%.3g}cos(2\\pi\\frac{%.3g}{%.3g})$   .     $\\phi_2=-\\frac{1}{(%.3g)^2}$", 1+exp(-tau1), f, fs, 1+exp(-tau1)))
  eq2 <- (sprintf("$\\phi_1=\\frac{2}{%.3g}cos(2\\pi\\frac{%.3g}{%.3g})$   .     $\\phi_2=-\\frac{1}{(%.3g)^2}$", 1+exp(-tau2), f, fs, 1+exp(-tau2)))
  #print(eq1)
  q <- ggplot(ggdata) +
    geom_ribbon(aes(x=f, ymin=0, ymax=S1, color="A"), fill=color1, alpha=0.2) +
    geom_line(aes(x=f, y=S1), color=color1) +
    geom_ribbon(aes(x=f, ymin=0, ymax=S2, color="B"), fill=color2, alpha=0.2) +
    geom_line(aes(x=f, y=S2), color=color2) +
    #guides(color=guide_legend(title=NULL)) +
    scale_color_discrete(labels=lapply(c(eq1, eq2), TeX)) + 
    labs(colour="AR(2) model", x="Frequency [Hertz]", y="Normalized PSD") +
    csda_theme()
  q
}
