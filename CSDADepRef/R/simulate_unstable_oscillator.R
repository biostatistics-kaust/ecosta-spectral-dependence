
#' Simulate unstable oscillators
#' 
#' For the purposes of this package, we use the definition of [A](doi:dasar-link):
#' an unstable oscillator is a oscillator which frequency cauld vary around
#' a central frequency.
#' Here, we simulate it through an AR(2) model [A](www.doi.org/10.1109/ACCESS.2020.3019243):
#' \loadmathjax
#' \mjdeqn{
#'    \phi_{1} = \frac {2}{1+e^{-\tau }}\cos({2\pi \omega^{*}})
#'  }{phi1 = 2 / (1+e^tau) cos(2*pi*w) }
#' \mjdeqn{
#'    \phi_{2} = -\frac {1}{\left ({1+e^{-\tau }}\right)^{2}}
#'  }{phi1 = -1 /(1+e^tau)^2 }
#' The power spectral density  of the AR(2) is
#' \mjdeqn{
#' S_{x}\left ({\omega }\right)
#' = \frac
#'    {\sigma _{\varepsilon }^{2}}
#'    {\left |{1-\phi_{1}e^{-2j\pi\omega }-\phi_{2}e^{-4j\pi \omega }}\right |^{2}}
#' }{S(w) = s^2/abs(1-phi1*e^(-2j*pi*w)-phi2*e^(-4j*pi*w))}
#' @keywords oscillator, AR2
#' @export
#' @examples 
#' \dontrun{
#'   simulate_unstable_oscillator(100, 10, 4, 30) # Sample 100 points of an 10Hz oscillator with $tau=4$ sampled at 30Hz
#' }
simulate_unstable_oscillator <- function(N, f, tau=1, fs=200, sigma=1) {
  stopifnot(f < 0.5 * fs)
  stopifnot(tau >= 0)
  stopifnot(sigma > 0)
  phi_1 <- 2 / (1 + exp(-tau)) * cos(2 * pi * f / fs)
  phi_2 <- -1 / (1 + exp(-tau)) ^ 2
  if(is.null(N)){
    return(list(phi1=phi_1, phi2=phi_2))
  }
  y <- as.numeric(arima.sim(model=list(ar=c(phi_1, phi_2)), n=N, sd=1))
  (y - mean(y)) / sd(y) * sigma
}

