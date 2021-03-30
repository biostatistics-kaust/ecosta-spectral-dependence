#' Provides a smooth estimator of the spectrum of a signal
#' @export
estimate_smooth_spectrum <- function(x, fcut=0.05) {
  Sx <- spec.pgram(x, taper=0, log="no", plot=FALSE)$spec
  coeffs <- signal::butter(3, fcut, type="low")
  pmax(0, signal::filtfilt(coeffs, Sx))
}


