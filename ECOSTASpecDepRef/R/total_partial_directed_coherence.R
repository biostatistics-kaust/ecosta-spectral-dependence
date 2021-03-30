complex_cov_to_cor <- function (A){
  p <- ncol(A)
  a <- sqrt(1/diag(A))
  Pi <- Conj(a) * A * rep(a, each = p)
  Pi
}

normalize_PDC_matrix <- function (A){
    p <- ncol(A)
    A <- (diag(p) - A)
    inv_a <- 1 / sqrt(colSums(abs(A) ^ 2))
    Pi <- A * rep(inv_a, each = p)
    Pi
}

normalize_GPDC_matrix = function (A, Sigma){
    p <- ncol(A)
    inv_sigmas <- sqrt(1 / diag(Sigma))
    A <- (diag(p) - A) * inv_sigmas
    inv_a <- 1 / sqrt(colSums(abs(A) ^ 2))
    #Pi <- A * rep(inv_a, each = p) * Conj(rep(inv_a, each = p))
    Pi <- A * rep(inv_a, each = p)
    Pi
}

phi_spectrum <- function(Phi, w) {
    phi_rows <- nrow(Phi)
    phi_cols <- ncol(Phi)
    lags <- (1:(phi_cols * phi_rows) - 1) %/% (phi_rows * phi_rows)
    lags <- matrix(lags, nrow=phi_rows) + 1
    complex_exponentials_lags <- exp(-2i * pi * lags * w) * Phi
    complex_exponentials_lags <- matrix(complex_exponentials_lags, nrow=phi_rows * phi_rows)
    #print(sprintf("phi_cols > phi_rows:: %s %s", phi_cols, phi_rows))
    if(phi_cols > phi_rows){
        A <- t(matrix(rowSums(complex_exponentials_lags),nro=phi_rows))
    }else{
        A <- t(matrix(complex_exponentials_lags, nrow=phi_rows))
    }
    A
}

estimate_partial_directed_coherence <- function(Phi, Sigma, w, type="PDC") {
    eps_Phi <- + 1e-10 * diag(nrow(Phi))
    A <- phi_spectrum(Phi, w)
    Pi <- if(type=="PDC") {normalize_PDC_matrix(A)}
          else if(type=="GPDC") {normalize_GPDC_matrix(A, Sigma)}
          else {NULL}
    Pi
}

abs_partial_directed_coherence_freqs <- function(Phi, Sigma, freqs=(0:100)/200, type="PDC") {
  Sx <- array(NA, dim=c(length(freqs), nrow(Phi), ncol(Phi)))
  for(i in 1:length(freqs)){
    Sx[i,,] <- abs(estimate_partial_directed_coherence(Phi, Sigma, freqs[i], type))
  }
  Sx
}

total_abs_partial_directed_coherence_freqs <- function(Phi, Sigma, freqs=(0:100)/200, type="PDC") {
  Sx <- 0
  Sx_max <- Phi * 0
  for(i in 1:length(freqs)){
    Px <- abs(estimate_partial_directed_coherence(Phi, Sigma, freqs[i], type))
    Sx <- Sx + Px
    Sx_max <- pmax(Sx_max, Px)
    P <- abs(estimate_partial_directed_coherence(Phi, Sigma, freqs[i], type))
    if(any(abs(P) > 1)){
      print(c("===>", max(P), min(P)))
      print(Phi)
      print(Sigma)
      print(freqs[i])
      print(type)
    }
  }
  Sx / length(freqs)
}


partial_directed_coherence <- function(data, var_order, fs, number_freq_points=200, type="PDC"){
	model <- var_coefficients(data, var_order=var_order)
  freqs <- (0:(number_freq_points-1)) / (2 * number_freq_points)
  pdcoh <- abs_partial_directed_coherence_freqs(model$Phi, model$Sigma, freqs, type)
  list(
    freqs=freqs * fs,
    part_direct_coherence=pdcoh
  )
}

#' @export
total_partial_directed_coherence <- function(data, var_order, fs, number_freq_points=200, type="PDC", freq_region=c(0, 0.5*fs), var_type="olsVARS"){
  stopifnot(all(freq_region >= 0))
  stopifnot(all(freq_region <= 0.5 * fs))
  stopifnot(freq_region[2] > freq_region[1])
	model <- var_coefficients(data, var_order=var_order, var_type=var_type)
  freqs <- (0:(number_freq_points-1)) / (2 * number_freq_points) * (freq_region[2] - freq_region[1]) / fs + freq_region[1] / fs
  pdcoh <- total_abs_partial_directed_coherence_freqs(model$Phi, model$Sigma, freqs, type)
  pdcoh
}











