#' @export
compare_three_latent_process <- function(Z1, Z2, Z3, X1, X2, X3, fs, mixture_matrix, lag_matrix=NULL, band_A="delta", band_B="beta", band_C="gamma") {  
  N <- length(X1)
  max_lags <- ceiling(N * 0.25)
  color_red <- "#db5067"
  color_blue <- "#2b76b2"
  color_gold <- "#e4a027"
  color_green <- "#2a9f78"
  
  X_1A <- extract_rhythms(X1, fs=fs, rhythms=c(band_A))[[band_A]]
  X_2A <- extract_rhythms(X2, fs=fs, rhythms=c(band_A))[[band_A]]
  X_3A <- extract_rhythms(X3, fs=fs, rhythms=c(band_A))[[band_A]]
  X_1B <- extract_rhythms(X1, fs=fs, rhythms=c(band_B))[[band_B]]
  X_2B <- extract_rhythms(X2, fs=fs, rhythms=c(band_B))[[band_B]]
  X_3B <- extract_rhythms(X3, fs=fs, rhythms=c(band_B))[[band_B]]
  X_1C <- extract_rhythms(X1, fs=fs, rhythms=c(band_C))[[band_C]]
  X_2C <- extract_rhythms(X2, fs=fs, rhythms=c(band_C))[[band_C]]
  X_3C <- extract_rhythms(X3, fs=fs, rhythms=c(band_C))[[band_C]]
  
  plot_net <- plot_network_estimator(
    X=cbind(X1, X2, X3),
    fs=fs,
    label_id=16,
    selected_bands=c(band_A, band_B, band_C),
    network_estimator=coherence_callback(max_lags=150, sign_lags="-"),
    padding=0.3,
    threshold_rel_max=0.1,
    metric_name="COH",
    is_symmetric=TRUE
  )
  
  coh_plot <- function(X_1A, X_2B, ch1, ch2, band, label_id, color="#e4a027", max_lags=200) dep_coherence_lagged_plot(
    X_1A, X_2B,
    max_lags=max_lags,
    label_id=label_id,
    xlabel=sprintf("$X_%s$", ch1),
    ylabel=sprintf("$X_%s$", ch2),
    color=color,
    corr_var_name="\\\\textbf{C}orr\\[\\cdot, \\cdot\\]",
    x_var_name=sprintf("X_{%s,%s}", ch1, band),
    y_var_name=sprintf("X_{%s,%s}", ch2, band),
    add_corr_definition_at_title=FALSE,
    use_dbl_line_title=FALSE
  )
  
  
  # mixture_matrix=mixture_matrix, lag_matrix=lag_matrix
  latent_label <- function(row, category="Latent source") sprintf(
    "%s: $Z_%s(t)$",
    category,
    row
  )
  print("mixture_matrix")
  print(mixture_matrix)
  
  component_label <- function(row, category="Component") sprintf(
    "%s: $X_%s(t)=%s\\cdot{}Z_1(t%s) + %s\\cdot{}Z_2(t%s) + %s\\cdot{}Z_3(t%s) + \\epsilon_{%s}(t)$ ",
    category,
    row,
    mixture_matrix[row, 1], 
    (if(is.null(lag_matrix) || (lag_matrix[row, 1] == 0)) "" else sprintf("%s", lag_matrix[row, 1])),
    mixture_matrix[row, 2], 
    (if(is.null(lag_matrix) || (lag_matrix[row, 2] == 0)) "" else sprintf("%s", lag_matrix[row, 2])),
    mixture_matrix[row, 3], 
    (if(is.null(lag_matrix) || (lag_matrix[row, 3] == 0)) "" else sprintf("%s", lag_matrix[row, 3])),
    row
  )
  plot_Z1 <- dep_line_plot(Z1, fs, label_id=1, title=latent_label(1), color="#2c9e78")
  plot_Z2 <- dep_line_plot(Z2, fs, label_id=2, title=latent_label(2), color="#2c9e78")
  plot_Z3 <- dep_line_plot(Z3, fs, label_id=3, title=latent_label(3), color="#2c9e78")
  plot_X1 <- dep_line_plot(X1, fs, label_id=4, title=component_label(1), color="#d6641e")
  plot_X2 <- dep_line_plot(X2, fs, label_id=5, title=component_label(2), color="#d6641e")
  plot_X3 <- dep_line_plot(X3, fs, label_id=6, title=component_label(3), color="#d6641e")
  #
  plot_X12_A <- coh_plot(X_1A, X_2A, ch1=1, ch2=2, band=band_A, label_id=7, color="#8c564b")
  plot_X13_A <- coh_plot(X_1A, X_3A, ch1=1, ch2=3, band=band_A, label_id=8, color="#8c564b")
  plot_X23_A <- coh_plot(X_2A, X_3A, ch1=2, ch2=3, band=band_A, label_id=9, color="#8c564b")
  #
  plot_X12_B <- coh_plot(X_1B, X_2B, ch1=1, ch2=2, band=band_B, label_id=10, color="#ff7f00")
  plot_X13_B <- coh_plot(X_1B, X_3B, ch1=1, ch2=3, band=band_B, label_id=11, color="#ff7f00")
  plot_X23_B <- coh_plot(X_2B, X_3B, ch1=2, ch2=3, band=band_B, label_id=12, color="#ff7f00")
  #
  plot_X12_C <- coh_plot(X_1C, X_2C, ch1=1, ch2=2, band=band_C, label_id=13, color="#66a61e")
  plot_X13_C <- coh_plot(X_1C, X_3C, ch1=1, ch2=3, band=band_C, label_id=14, color="#66a61e")
  plot_X23_C <- coh_plot(X_2C, X_3C, ch1=2, ch2=3, band=band_C, label_id=15, color="#66a61e")
  #
  gridExtra::grid.arrange(grobs=list(
      plot_Z1,
      plot_Z2,
      plot_Z3,
      
      plot_X1,
      plot_X2,
      plot_X3,
      
      plot_X12_A,
      plot_X13_A,
      plot_X23_A,
      
      plot_X12_B,
      plot_X13_B,
      plot_X23_B,
      
      plot_X12_C,
      plot_X13_C,
      plot_X23_C,
      
      plot_net
      
    ), ncol=8,
    layout_matrix=rbind(
      rep(c(1, 4), each=4),
      rep(c(2, 5), each=4),
      rep(c(3, 6), each=4),
      rep(c( 16,  7,  8,  9), each=2),
      rep(c( 16, 10, 11, 12), each=2),
      rep(c( 16, 13, 14, 15), each=2)
    )
  )
}
