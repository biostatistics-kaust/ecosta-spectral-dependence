#' Extract EEG brain rhythms from a time series
#' 
#' @keywords EEG
#' @export
#' @examples 
#' \dontrun{
#'   simulate_unstable_oscillator(100, 10, 4, 30) # Sample 100 points of an 10Hz oscillator with $tau=4$ sampled at 30Hz
#'  x <- simulate_unstable_oscillator(1000, 1, 1, 200);
#'  y <- extract_rhythms(x, 200)
#'  plot(y$delta); lines(y$theta, col="red"); lines(x, col="blue")
#' }
extract_rhythms <- function(x, fs, rhythms=NULL, print_filters=TRUE) {
  rhythms_info <- eeg_rhythms()
  valid_rhythms <- names(rhythms_info)
  if(is.null(rhythms)){
    rhythms <- valid_rhythms
  }
  #
  result <- list()
  for(rhythm in rhythms){
    stopifnot(rhythm %in% valid_rhythms)
    frequencies <- rhythms_info[[rhythm]]
    if(any(frequencies >= (0.5 * fs))){
      warning(paste("Rhythm", rhythm, "requires at least fs =", frequencies[2] * 2, "Hz. Provided: ", fs))
    }
    coeffs <- signal::fir1(10, frequencies/(0.5 * fs),  type="pass")
    #coeffs <- signal::butter(3, frequencies/(0.5 * fs),  type="pass")
    if(print_filters){
      print(sprintf(":: RHYTHM: %s", rhythm))
      print(coeffs)
    }
    y <- x
    for(i in 1:5){
      L <- length(x)
      y <- as.numeric(signal::filter(coeffs, c(y, y))[(L+1):(2 * L)])
    }
    result[[rhythm]] <- y
  }
  result
}


