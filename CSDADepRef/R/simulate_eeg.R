#' Simulate EEG
#' 
#' @export
#' @examples 
#' \dontrun{
#'   y <- simulate_EEG(1000) # Sample 1000 points at 200Hz
#'   spectrum(simulate_EEG(1000, electrical_line=50,), log="no")
#' }
simulate_EEG <- function(N, fs=200, f=c(0,2,4,10,20,34,45), tau=c(1,1,1,1,1,1,1), sigma=c(3,1,1,1,1,1,1), electrical_line=NULL) {
  stopifnot(all(f < 0.5 * fs))
  stopifnot(all(tau >= 0))
  stopifnot(all(sigma > 0))
  stopifnot(is.null(electrical_line) || electrical_line > 40)
  stopifnot(length(f) == length(tau) && length(tau) == length(sigma))
  y <- 0
  for(i in 1:length(f)){
    y <- y + simulate_unstable_oscillator(N, f=f[i], tau=tau[i], sigma=sigma[i], fs=fs)
  }
  if(!is.null(electrical_line)){
    y <- y + simulate_unstable_oscillator(N, f=electrical_line, tau=10, sigma=0.5*mean(sigma))
  }
  y
}

