#' @export
run_example_1 <- function() {
  fs = 128
  N = fs * 60
  T0 = fs * 10
  
  #"Fp1" "Fp2" "F3"  "F4"  "C3"  "C4"  "P3"  "P4"  "O1"  "O2"  "F7"  "F8" "T7"  "T8"  "P7"  "P8"  "Fz"  "Cz"  "Pz"
  x <- read.csv("../dataset/10.21227--rzfh-zn36/decompress-data/ADHD/v1p-cleaned.csv")[1:N+T0, "F3"] * 1e-9
  y <- read.csv("../dataset/10.21227--rzfh-zn36/decompress-data/ADHD/v1p-cleaned.csv")[1:N+T0, "F4"] * 1e-9
  x <- read.csv("../dataset/10.21227--rzfh-zn36/decompress-data/Control/v41p-cleaned.csv")[1:N+T0, "F3"] * 1e-9
  y <- read.csv("../dataset/10.21227--rzfh-zn36/decompress-data/Control/v41p-cleaned.csv")[1:N+T0, "F4"] * 1e-9
  
  #x <- simulate_EEG(N, fs=fs)
  #y <- simulate_EEG(N, fs=fs)
  q <- compare_time_rhythms(x, y, target_rhythm="alpha", flim=c(3, 20), channel_label1="F3", channel_label2="F4")
  csda_save_image(q, "time-comparison-alpha")
  q <- compare_time_rhythms(x, y, target_rhythm="all", flim=c(3, 20))
  csda_save_image(q, "time-comparison-all")
  q <- compare_nbcoherence(x, y, target_rhythm="alpha", max_plots=2, max_lags_calculate=100)
  csda_save_image(q, "coh-comparison-alpha")
  q <- compare_nbcoherence(x, y, target_rhythm="all", max_plots=2, max_lags_calculate=100)
  csda_save_image(q, "coh-comparison-all")
}
