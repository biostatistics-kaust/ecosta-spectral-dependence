#' @param sign_lags    We have two interpretations of x(t-14)
#'                   lag is defined as x(t-h), lag is therefore 14 (sign_lags="-")
#'                   or lag is defined as x(t+h), and lag then is -14 (sign_lags="+")
#'                   By default, and following the paper's notation, sign_lags="-"
#' @export
coherence_corr <- function(x, y, max_lags=10, sign_lags="-") {
  stopifnot(length(x) == length(y))
  stopifnot(max_lags < 0.7 * length(x))
  stopifnot(sign_lags %in% c("+", "-"))
  
  # How lag algorithm works (1-indices): # Revise for (0-indices langs) ### Assuming x(t+h)
  # Vector indices (at t=0): | length=15
  #    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
  # Time indices (at t=0) that it represents:
  #   -7    -6    -5    -4    -3    -2    -1     0     1     2     3     4     5     6     7
  # Time indices (at t=+3):
  #   -4    -3    -2    -1     0     1     2     3     4     5     6     7
  # Time indices (at t=-2):
  # ??-9  ??-8    -7    -6    -5    -4    -3    -2    -1     0     1     2     3     4     5     6     7
  lagged_positive_correlation <- sapply(1:max_lags, function(lag)  # x(t) <-> y(t + l)
    cor(x[1:(length(y)-lag)], y[(lag+1):length(x)])
  )
  lagged_negative_correlation <- sapply(max_lags:0, function(lag)  # x(t + l) <-> y(t)
    cor(x[(lag+1):length(x)], y[1:(length(y)-lag)] )
  )
  #
  corr_indices <- -max_lags:max_lags
  corr_lag <- c(lagged_negative_correlation, lagged_positive_correlation) ^ 2
  if(sign_lags == "-") {
    #NO need to invert corr_indices, it's a symmetric array
    corr_lag <- rev(corr_lag)
  }
  index_maximum <- which.max(corr_lag)
  corr_idx_maximum <- corr_indices[index_maximum]
  corr_maximum <- corr_lag[index_maximum]
  #
  list(
    correlation_values=corr_lag,
    correlation_indices=corr_indices,
    coherence=corr_maximum,
    max_corr_lag=corr_idx_maximum
  )
}
