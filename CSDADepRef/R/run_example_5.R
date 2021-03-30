#' @export
base_data_example_5 <- function(add_lag=FALSE) {
  set.seed(1005)
  if(add_lag){
    set.seed(10055)
  }
  fs <- 128
  #N = 130
  N = fs * 60
  #N = fs * 2
  N = fs * 10
  #T0 = fs * 10
  mixture_matrix <- rbind(
    c(2.0, 1.0, 0.5),
    c(0.0, 1.0, 0.0),
    c(0.0, 0.0, 1.0)
  )
  mixture_matrix <- rbind(
    c(0.0, 1.0, 0.5),
    c(1.0, 0.0, 1.0),
    c(0.0, 0.0, 1.0)
  )
  lag_matrix <- NULL
  if(add_lag){
    # X3 -> X1, 2
    # X2 -> X1
    lag_matrix <- -1 * rbind(
      c(1, 1, 1),
      c(0, 3, 0),
      c(0, 0, 2)
    )
    #lag_matrix=t(lag_matrix)
    print("lag_matrix")
    print(lag_matrix)
  }
  R <- simulate_sparse_contemporaneous_mixture(
    N, mixture_matrix=mixture_matrix, fs=fs,
    f=c(2, 10, 40), tau=c(3, 3, 3), sigma=c(2, 2, 2),
    lag_matrix=lag_matrix,
    sigma_epsilon=c(1, 1, 1)
    #sigma_epsilon=c(0, 0)
  )
  
  list(
    Z=R[[1]],
    X=R[[2]],
    fs=fs,
    mixture_matrix=mixture_matrix
  )
}


#' @export
base_data_example_5_second <- function(add_lag=FALSE) {
  set.seed(1005)
  if(add_lag){
    set.seed(10055)
  }
  fs <- 128
  #N = 130
  N = fs * 60
  #N = fs * 2
  N = fs * 10
  #T0 = fs * 10
  mixture_matrix <- rbind(
    c(0.0, 1.0, 0.5),
    c(1.0, 1.0, 1.0),
    c(0.0, 0.0, 1.0)
  )
  lag_matrix <- NULL
  if(add_lag){
    # X3 -> X1, 2
    # X2 -> X1
    lag_matrix <- -1 * rbind(
      c(1, 1, 1),
      c(0, 3, 0),
      c(0, 0, 2)
    )
    #lag_matrix=t(lag_matrix)
    print("lag_matrix")
    print(lag_matrix)
  }
  R <- simulate_sparse_contemporaneous_mixture(
    N, mixture_matrix=mixture_matrix, fs=fs,
    f=c(2, 10, 40), tau=c(3, 3, 3), sigma=c(2, 2, 2),
    lag_matrix=lag_matrix,
    sigma_epsilon=c(1, 1, 1)
    #sigma_epsilon=c(0, 0)
  )
  
  list(
    Z=R[[1]],
    X=R[[2]],
    fs=fs,
    mixture_matrix=mixture_matrix
  )
}


#' @export
run_example_5 <- function() {
  Q <- base_data_example_5()
  compare_three_latent_process(
    Q$Z[1,], Q$Z[2,], Q$Z[3,],
    Q$X[1,], Q$X[2,], Q$X[3,],
    fs=Q$fs,
    mixture_matrix=Q$mixture_matrix,
    lag_matrix=NULL,
    band_A="delta",
    band_B="alpha",
    band_C="gamma"
  )
}

#' @export
run_example_5b <- function() {
  Q <- base_data_example_5()
  print(("========== COH =========="))
  plot_coh <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=c("delta", "alpha", "gamma"),
    network_estimator=coherence_callback(max_lags=150, sign_lags="-"),
    metric_name="COH",
    is_symmetric=TRUE
  ) + labs(title="A) Coherence network")
  plot(plot_coh)
  
  print(("========== P-COH =========="))
  plot_pcoh <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=c("delta", "alpha", "gamma"),
    network_estimator=partial_coherence_callback(max_lags=150, sign_lags="-"),
    threshold_rel_max=0.20,
    #padding=0.2,
    metric_name="PCOH",
    is_symmetric=TRUE
  ) + labs(title="B) Partial coherence network")
  plot(plot_pcoh)
  
  q <- gridExtra::grid.arrange(grobs=list(plot_coh, plot_pcoh), ncol=2,
    layout_matrix=rbind(
      c(1, 2)
    )
  )
}

#' @export
run_example_5c <- function() {
  Q <- base_data_example_5(add_lag=TRUE)
  #
  #
  var_order <- 2
  Z <- Q$Z
  X <- Q$X
  
  print(("========== COH =========="))
  plot_coh <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=c("delta", "alpha", "gamma"),
    network_estimator=coherence_callback(max_lags=150, sign_lags="-"),
    threshold_rel_max=0.01,
    metric_name="COH",
    is_symmetric=TRUE
  ) + labs(title="A) Coherence network")
  plot(plot_coh)
  
  print(("========== P-COH =========="))
  plot_pcoh <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=c("delta", "alpha", "gamma"),
    network_estimator=partial_coherence_callback(max_lags=150, sign_lags="-"),
    threshold_rel_max=0.01,
    metric_name="PCOH",
    is_symmetric=TRUE
  ) + labs(title="B) Partial coherence network")
  plot(plot_pcoh)

  print(("========== PD-COH =========="))
  plot_pdcoh_var2 <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=c("delta", "alpha", "gamma"),
    network_estimator = function(X, fs, selected_band) total_partial_directed_coherence(
      X,
      fs=fs,
      freq_region=eeg_rhythms()[[selected_band]],
      var_order=3,
      number_freq_points=200,
      var_type="RawLSE",
      type="PDC"
    ),
    metric_name="|PDC|",
    is_symmetric=FALSE
  ) + labs(title="C) PDC:VAR(3) ")
  plot(plot_pdcoh_var2)

  plot_pdcoh_var7 <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=c("delta", "alpha", "gamma"),
    network_estimator = function(X, fs, selected_band) total_partial_directed_coherence(
      X,
      fs=fs,
      freq_region=eeg_rhythms()[[selected_band]],
      var_order=15,
      number_freq_points=200,
      var_type="RawLSE",
      type="PDC"
    ),
    metric_name="|PDC|",
    is_symmetric=FALSE
  ) + labs(title="D) PDC:VAR(15) ")
  plot(plot_pdcoh_var7)

  plot_pdcoh_lassle2 <-plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=c("delta", "alpha", "gamma"),
    network_estimator = function(X, fs, selected_band) total_partial_directed_coherence(
      X,
      fs=fs,
      freq_region=eeg_rhythms()[[selected_band]],
      var_order=3,
      number_freq_points=200,
      var_type="LASSLE",
      type="PDC"
    ),
    metric_name="|PDC|",
    is_symmetric=FALSE
  ) + labs(title="E) PDC:LASSLE(3) ")
  plot(plot_pdcoh_lassle2)

  plot_pdcoh_lassle7 <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=c("delta", "alpha", "gamma"),
    network_estimator = function(X, fs, selected_band) total_partial_directed_coherence(
      X,
      fs=fs,
      freq_region=eeg_rhythms()[[selected_band]],
      var_order=15,
      number_freq_points=200,
      var_type="LASSLE",
      type="PDC"
    ),
    is_symmetric=FALSE
  ) + labs(title="F) PDC:LASSLE(15) ")
  plot(plot_pdcoh_lassle7)

  gridExtra::grid.arrange(grobs=list(
    plot_coh, plot_pcoh,
    plot_pdcoh_var2, plot_pdcoh_var7,
    plot_pdcoh_lassle2, plot_pdcoh_lassle7
    ), ncol=2,
    layout_matrix=rbind(
      c(1, 2),
      c(3, 4),
      c(5, 6)
    )
  )
}


#' @export
run_example_5d <- function() {
  Q <- base_data_example_5_second()
  print(("========== COH =========="))
  plot_coh <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=c("delta", "alpha", "gamma"),
    network_estimator=coherence_callback(max_lags=150, sign_lags="-"),
    metric_name="COH",
    is_symmetric=TRUE
  ) + labs(title="A) Coherence network")
  plot(plot_coh)
  
  print(("========== P-COH =========="))
  plot_pcoh <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=c("delta", "alpha", "gamma"),
    network_estimator=partial_coherence_callback(max_lags=150, sign_lags="-"),
    threshold_rel_max=0.20,
    #padding=0.2,
    metric_name="PCOH",
    is_symmetric=TRUE
  ) + labs(title="B) Partial coherence network")
  plot(plot_pcoh)
  
  q <- gridExtra::grid.arrange(grobs=list(plot_coh, plot_pcoh), ncol=2,
    layout_matrix=rbind(
      c(1, 2)
    )
  )
}
