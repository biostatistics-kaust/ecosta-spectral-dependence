RHYTHMS_INFO <- list(
  delta=c(0.5, 4),
  theta=c(4, 8),
  alpha=c(8, 12),
  beta=c(12, 30),
  gamma=c(30, 45),
  highgamma=c(55, 95),
  all=c(0.5, 45)
)

VALID_RHYTHMS <- names(RHYTHMS_INFO)

#' @export
eeg_rhythms <- function() RHYTHMS_INFO
