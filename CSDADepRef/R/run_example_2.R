#' @export
run_example_2 <- function() {
  set.seed(1002)  
  fs = 128
  #N = 130
  N = fs * 60
  #N = fs * 2
  N = fs * 10
  #T0 = fs * 10
  mixture_matrix <- rbind(
    c(0.1, 0.5),
    c(0, 2)
  )
  
  R <- simulate_sparse_contemporaneous_mixture(
    N, mixture_matrix=mixture_matrix, fs=fs,
    f=c(2, 20), tau=c(3, 3), sigma=c(2, 2),
    sigma_epsilon=c(1, 1)
    #sigma_epsilon=c(0, 0)
  )
  Z <- R[[1]]
  X <- R[[2]]
  
  compare_two_latent_process(
    Z1=Z[1, ],
    Z2=Z[2, ],
    X1=X[1, ],
    X2=X[2, ],
    fs=fs,
    band_A="delta",
    band_B="beta",
    mixture_matrix=mixture_matrix
  )
}
