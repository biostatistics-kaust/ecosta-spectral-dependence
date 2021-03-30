#' @export
estimate_true_ar2_spectrum <- function(M, fs, phi1, phi2, normalized=TRUE) {
  f <- ((1:M) - 1) / M * (0.5)
  S <- 1 / abs(
    1
    - phi1 * exp(-2i * pi * f)
    - phi2 * exp(-4i * pi * f)
  ) ^ 2
  if(normalized){
    S <- S / max(S)
  }
  list(f=f * fs, S=S)
}
