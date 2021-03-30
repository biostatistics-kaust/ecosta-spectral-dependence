
#' \donotrun{
#'  mixture_matrix <- rbind(
#'     c(2, 1),
#'     c(0, 0.5)
#'   )
#'   Y <- simulate_sparse_contemporaneous_mixture(N, mixture_matrix=mixture_matrix, fs=fs, f=c(2, 10), tau=c(5, 5), sigma=c(2, 2), sigma_epsilon=c(1, 1))
#'   plot(Y[1,], type="l")
#'   lines(Y[2,])  
#' }
#' @export
simulate_sparse_contemporaneous_mixture <- function(N, mixture_matrix, fs, lag_matrix=NULL, f=c(2, 10), tau=c(1, 1), sigma=c(2, 2), sigma_epsilon=c(1, 1)){
  stopifnot(N > 100)
  stopifnot(fs >= 1)
  stopifnot(all(f < 0.5 * fs))
  stopifnot(all(tau >= 0))
  stopifnot(all(sigma > 0))
  stopifnot(length(f) == length(tau) && length(tau) == length(sigma))
  #stopifnot(nrow(mixture_matrix) == ncol(mixture_matrix) && nrow(mixture_matrix) == length(sigma))
  #stopifnot(is.null(lag_matrix) || (nrow(lag_matrix) == ncol(lag_matrix) && nrow(lag_matrix) == length(sigma)))
  stopifnot(ncol(mixture_matrix) == length(sigma))
  stopifnot(is.null(lag_matrix) || (ncol(mixture_matrix) == ncol(lag_matrix) && nrow(mixture_matrix) == nrow(lag_matrix)))
  stopifnot(nrow(mixture_matrix) == length(sigma_epsilon))
  #E <- matrix(NA, ncol=N, nrow=length(f))
  #G <- matrix(NA, ncol=N, nrow=length(f))
  n_col <- nrow(mixture_matrix)
  n_row <- ncol(mixture_matrix)
  #print(c(n_row, n_col))
  E <- matrix(NA, nrow=n_col, ncol=N)
  G <- matrix(NA, nrow=n_row, ncol=N)
  #print(c(length(f), -1, n_col))
  #print(N)
  #print((mixture_matrix))
  #print(dim(G))
  #print("===")
  Nf <- length(f)
  for(i in 1:n_row){
    #print(c(i, f[i], tau[i], sigma[i]))
    G[i, ] <- simulate_unstable_oscillator(N, f=f[i], tau=tau[i], sigma=sigma[i], fs=fs)# + rnorm(N, sd=sigma_epsilon[i])
  }
  for(i in 1:n_col){
    E[i, ] <- rnorm(N, sd=1) * sigma_epsilon[i] # sigma_epsilon[i] is outside, because it could be zero.
    #print(i)
    #print(E[i,])
  }
  #print("=============================>"); print(sigma_epsilon); print(lag_matrix)
  if(is.null(lag_matrix)){
    Y <- (mixture_matrix %*% G) + E
    ## Result in the shape (time x channels)
    #t(Y)
  }else{
    # How lag algorithm works (1-indices): # Revise for (0-indices langs)
    # Vector indices (at t=0): | length=15
    #    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
    # Time indices (at t=0) that it represents:
    #   -7    -6    -5    -4    -3    -2    -1     0     1     2     3     4     5     6     7
    # Time indices (at t=+3):
    #   -4    -3    -2    -1     0     1     2     3     4     5     6     7
    # Time indices (at t=-2):
    # ??-9  ??-8    -7    -6    -5    -4    -3    -2    -1     0     1     2     3     4     5     6     7
    # |+----|+------ Do not exist in real vector
    # After selecting::  lower_bound_idx=2   upper_bound_idx=15-3=12 
    # Vector indices (at t=0): | length=15
    #    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
    # Time indices (at t=0):     lower_bound_idx + indices
    #    x     x    -5    -4    -3    -2    -1     0     1     2     3     4     x     x     x
    # Time indices (at t=+3):    lower_bound_idx + 3 + indices
    #    x     x    -2    -1     0     1     2     3     4     5     6     7
    # Time indices (at t=-2):    lower_bound_idx - 2 + indices
    #    x     x    -7    -6    -5    -4    -3    -2    -1     0     1     2     x     x     x     x     x
    # |+----|+------ Do not exist in real vector
    # Length: 
    lower_bound_idx <- min(lag_matrix)
    lower_bound_idx <- if(lower_bound_idx < 0) (-lower_bound_idx) else 0
    upper_bound_idx <- min(N, N - max(lag_matrix))
    #TODO: check what happen if upper_bound_idx < 0
    N1 <- upper_bound_idx - lower_bound_idx
    Y <- matrix(NA, ncol=N1, nrow=Nf)
    for(i in 1:Nf){
      #print(i)
      y_i <- E[i,1:N1]
      for(j in 1:Nf){
        lag <- lag_matrix[i, j]
        idx <- (lower_bound_idx + lag) + (1:N1)
        #print(c( idx[1], idx[length(idx)-1] ))
        #print(c( lower_bound_idx, N1, length(y_i), length(G[i, idx]) ))
        y_i <- y_i + mixture_matrix[i, j] * G[j, idx]
      }
      Y[i, ] <- y_i
      #print(y_i)
    }
  }
  list(G, Y)
}
