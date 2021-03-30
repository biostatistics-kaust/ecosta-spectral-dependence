AR2_coefficients <- function(f, fs, tau=3)
  list(
    phi_1 = (2 / (1 + exp(-tau)) * cos(2 * pi * f / fs)),
    phi_2 = (-1 / (1 + exp(-tau)) ^ 2)
  )

#'@export
base_data_example_7 <- function(add_lag=FALSE) {
  set.seed(1005)
  if(add_lag){
    set.seed(10055)
  }
  
  fs <- 128
  #N = 130
  N = fs * 60
  #N = fs * 2
  N = fs * 10
  #N = fs * 2
  #T0 = fs * 10
  #library(BigVAR)
  k=4
  p=2
  B=matrix(0,nrow=k,ncol=p*k)

  phi_delta <- AR2_coefficients(f=2, fs=fs)
  phi_gamma <- AR2_coefficients(f=40, fs=fs)
  phi_beta <- AR2_coefficients(f=20, fs=fs)
  
  transition_matrix1 <- rbind(
    c(phi_beta$phi_1, .5, 0, 0),
    c(0, 0, 1, 0),
    c(0, 0, phi_delta$phi_1, 0),
    c(0, 0, 0, phi_gamma$phi_1)
  )
  transition_matrix2 <- rbind(
    c(phi_beta$phi_2, 0, 0, 0),
    c(0, 0, 0, 1),
    c(0, 0, phi_delta$phi_2, 0),
    c(0, 0, 0, phi_gamma$phi_2)
  )

  B[, 1:k] = transition_matrix1
  B[, 1:k + k] = transition_matrix2
  A <- BigVAR::VarptoVar1MC(B,p,k)
  Y <- t(BigVAR::MultVarSim(k, A, p, diag(c(1e-7, 1e-7, 1e-3, .1)), N))
  Y <- Y - rowMeans(Y)
  #
  list(
    X=Y,
    Y=Y,
    fs=fs,
    mixture_matrix=B
  )
}

#'@export
run_example_7 <- function() {
  options(width=140)

  Q <- base_data_example_7()
  selected_bands <- c("delta", "beta", "gamma")
  selected_bands <- c("delta", "gamma")
  threshold_rel_max <- 0.35
  
  #X1 <- Q$Y[1,]; X2 <- Q$Y[2,]; X3 <- Q$Y[3,]; X4 <- Q$Y[4,]
  print(("========== COH =========="))
  print("A) Coherence network")
  plot_coh <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=selected_bands,
    network_estimator=coherence_callback(max_lags=150, sign_lags="-"),
    threshold_rel_max=threshold_rel_max,
    return_metric=TRUE,
    is_symmetric=TRUE
  )
  plot_coh$plot <- plot_coh$plot + labs(title="A) Coherence network")
  plot(plot_coh$plot)
  
  
  print("B) Partial coherence network")
  print(("========== P-COH =========="))
  plot_pcoh <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=selected_bands,
    network_estimator=partial_coherence_callback(max_lags=150, sign_lags="-"),
    threshold_rel_max=threshold_rel_max,
    return_metric=TRUE,
    is_symmetric=TRUE
  )
  plot_pcoh$plot <- plot_pcoh$plot + labs(title="B) Partial coherence network")
  plot(plot_pcoh$plot)

  
  print(("========== PD-COH =========="))
  print("C) PDC:VAR(2) ")
  plot_pdcoh_var2 <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=selected_bands,
    network_estimator = function(X, fs, selected_band) total_partial_directed_coherence(
      X,
      fs=fs,
      freq_region=eeg_rhythms()[[selected_band]],
      var_order=2,
      number_freq_points=200,
      var_type="RawLSE",
      type="PDC"
    ),
    threshold_rel_max=threshold_rel_max,
    return_metric=TRUE,
    is_symmetric=FALSE
  )
  plot_pdcoh_var2$plot <- plot_pdcoh_var2$plot + labs(title="C) PDC:VAR(2) ")
  plot(plot_pdcoh_var2$plot)

  
  print("D) PDC:VAR(15) ")
  plot_pdcoh_var7 <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=selected_bands,
    network_estimator = function(X, fs, selected_band) total_partial_directed_coherence(
      X,
      fs=fs,
      freq_region=eeg_rhythms()[[selected_band]],
      var_order=15,
      number_freq_points=200,
      var_type="RawLSE",
      type="PDC"
    ),
    threshold_rel_max=threshold_rel_max,
    return_metric=TRUE,
    is_symmetric=FALSE
  )
  plot_pdcoh_var7$plot <- plot_pdcoh_var7$plot + labs(title="D) PDC:VAR(15) ")
  plot(plot_pdcoh_var7$plot)

  
  print("E) PDC:LASSLE(2) ")
  plot_pdcoh_lassle2 <-plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=selected_bands,
    network_estimator = function(X, fs, selected_band) total_partial_directed_coherence(
      X,
      fs=fs,
      freq_region=eeg_rhythms()[[selected_band]],
      var_order=2,
      number_freq_points=200,
      var_type="LASSLE",
      type="PDC"
    ),
    threshold_rel_max=threshold_rel_max,
    return_metric=TRUE,
    is_symmetric=FALSE
  )
  plot_pdcoh_lassle2$plot <- plot_pdcoh_lassle2$plot + labs(title="E) PDC:LASSLE(2) ")
  plot(plot_pdcoh_lassle2$plot)

  
  print("F) PDC:LASSLE(15) ")
  plot_pdcoh_lassle7 <- plot_network_estimator(
    X=t(Q$X),
    fs=Q$fs,
    selected_bands=selected_bands,
    network_estimator = function(X, fs, selected_band) total_partial_directed_coherence(
      X,
      fs=fs,
      freq_region=eeg_rhythms()[[selected_band]],
      var_order=15,
      number_freq_points=200,
      var_type="LASSLE",
      type="PDC"
    ),
    threshold_rel_max=threshold_rel_max,
    return_metric=TRUE,
    is_symmetric=FALSE
  )
  plot_pdcoh_lassle7$plot <- plot_pdcoh_lassle7$plot + labs(title="F) PDC:LASSLE(15) ")
  plot(plot_pdcoh_lassle7$plot)

  remove_unnecessary_cols <- function(dataframe, varname) {
    print(varname)
    print(dataframe)
    new_dataframe <- data.frame(start=dataframe$start, end=dataframe$end, band=dataframe$band)
    new_dataframe[[varname]] <- dataframe$weight
    new_dataframe
  }
  general_metrics <- Reduce(
    function(x, y, ...) merge(x, y, by=c("start", "end", "band"), all=TRUE, ...),
    lapply(list(
      list(plot_coh$data, "COH"),
      list(plot_pcoh$data, "PCOH"),
      list(plot_pdcoh_var2$data, "PDC_VAR2"), 
      list(plot_pdcoh_var7$data, "PDC_VAR7"),
      list(plot_pdcoh_lassle2$data, "PDC_LASSLE2"), 
      list(plot_pdcoh_lassle7$data, "PDC_LASSLE7")
    ), function(x) remove_unnecessary_cols(x[[1]], x[[2]]))
  )
  print(general_metrics)
  
  gridExtra::grid.arrange(grobs=list(
    plot_coh$plot, plot_pcoh$plot,
    plot_pdcoh_var2$plot, plot_pdcoh_var7$plot,
    plot_pdcoh_lassle2$plot, plot_pdcoh_lassle7$plot
    ), ncol=2,
    layout_matrix=rbind(
      c(1, 2),
      c(3, 4),
      c(5, 6)
    )
  )
}

