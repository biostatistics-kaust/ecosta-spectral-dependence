#' @export
compare_two_latent_process <- function(Z1, Z2, X1, X2, fs, mixture_matrix, lag_matrix=NULL, band_A="delta", band_B="beta") {  
  N <- length(X1)
  max_lags <- ceiling(N * 0.25)
  color_red <- "#db5067"
  color_blue <- "#2b76b2"
  color_gold <- "#e4a027"
  color_green <- "#2a9f78"
  
  # mixture_matrix=mixture_matrix, lag_matrix=lag_matrix
  latent_label <- function(row, category="Latent source") sprintf(
    "%s: $Z_%s(t)$",
    category,
    row
  )
  
  component_label <- function(row, category="Component") sprintf(
    "%s: $X_%s(t)=%s\\cdot{}Z_1(t%s) + %s\\cdot{}Z_2(t%s) + \\epsilon_{%s}(t)$ ",
    category,
    row,
    mixture_matrix[row, 1], 
    (if(is.null(lag_matrix) || (lag_matrix[row, 1] == 0)) "" else sprintf("%s", lag_matrix[row, 1])),
    mixture_matrix[row, 2], 
    (if(is.null(lag_matrix) || (lag_matrix[row, 2] == 0)) "" else sprintf("%s", lag_matrix[row, 2])),
    row
  )
  
  plot_Z1 <- dep_line_plot(Z1, fs, label_id=1, title=latent_label(1), color="#2c9e78")
  plot_Z2 <- dep_line_plot(Z2, fs, label_id=2, title=latent_label(2), color="#2c9e78")
  plot_X1 <- dep_line_plot(X1, fs, label_id=3, title=component_label(1), color="#d6641e")
  plot_X2 <- dep_line_plot(X2, fs, label_id=4, title=component_label(2), color="#d6641e")

  plot_R_12 <- dep_scatter_plot(
    X1, X2, 
    corr_var_name="\\\\textbf{C}orr\\[X_1,X_2\\]",
    label_id=5, 
    xlabel="X", ylabel="Y", 
    color="#2b76b2")
  
  X_1A <- extract_rhythms(X1, fs=fs, rhythms=c(band_A))[[band_A]]
  X_2A <- extract_rhythms(X2, fs=fs, rhythms=c(band_A))[[band_A]]
  X_1B <- extract_rhythms(X1, fs=fs, rhythms=c(band_B))[[band_B]]
  X_2B <- extract_rhythms(X2, fs=fs, rhythms=c(band_B))[[band_B]]
  
  plot_R_1A_2A <- dep_scatter_plot(
    X_1A, X_2A,
    corr_var_name="\\textbf{C}orr\\[X_{1,\\Omega_0}, X_{2,\\Omega_0}\\]",
    label_id=6, 
    xlabel="X_1", ylabel="X_2", 
    color="#2c77b1")
  
  plot_R_1A_2B <- dep_scatter_plot(
    X_1A, X_2B,
    corr_var_name="\\\\textbf{C}orr\\[X_{1,\\Omega_0}, X_{2,\\Omega_1}\\]",
    label_id=7, 
    xlabel="X_1", ylabel="X_2", 
    color="#2c77b1")
  
  plot_R_1B_2B <- dep_scatter_plot(
    X_1B, X_2B,
    corr_var_name="\\\\textbf{C}orr\\[X_{1,\\Omega_1}, X_{2,\\Omega_1}\\]",
    label_id=8, 
    xlabel="X_1", ylabel="X_2", 
    color="#2c77b1") +
  theme(plot.background = element_rect(fill = "#d2e5ef", colour="transparent"))
  
  plot_R_1B_2A <- dep_scatter_plot(
    X_1B, X_2A,
    corr_var_name="\\\\textbf{C}orr\\[X_{1,\\Omega_1}, X_{2,\\Omega_0}\\]",
    label_id=9, 
    xlabel="X_1", ylabel="X_2", 
    color="#2c77b1")
  
  
  plot_C_12 <- dep_coherence_lagged_plot(
    X1, X2,
    max_lags=max_lags,
    label_id=10,
    xlabel="$X_1$", ylabel="$X_2$", 
    color="#e4a027",
    corr_var_name="\\\\textbf{C}orr\\[\\cdot, \\cdot\\]",
    #corr_var_name="\\\\textbf{C}orr\\[X_{1}, X_{2}\\]",
    x_var_name="X_{1}",
    y_var_name="X_{2}",
    add_corr_definition_at_title=FALSE,
    use_dbl_line_title=TRUE
  )
  plot_C_1A_2A <- dep_coherence_lagged_plot(
    X_1A, X_2A,
    max_lags=max_lags,
    label_id=11,
    xlabel="$X_1$", ylabel="$X_2$", 
    color="#e4a027",
    corr_var_name="\\\\textbf{C}orr\\[\\cdot, \\cdot\\]",
    #corr_var_name="\\\\textbf{C}orr\\[X_{1,\\Omega_0}, X_{2,\\Omega_0}\\]",
    x_var_name="X_{1,\\Omega_0}",
    y_var_name="X_{2,\\Omega_0}",
    add_corr_definition_at_title=FALSE,
    use_dbl_line_title=TRUE
  )
  plot_C_1A_2B <- dep_coherence_lagged_plot(
    X_1A, X_2B,
    max_lags=max_lags,
    label_id=12,
    xlabel="$X_1$", ylabel="$X_2$", 
    color="#e4a027",
    corr_var_name="\\\\textbf{C}orr\\[\\cdot, \\cdot\\]",
    #corr_var_name="\\\\textbf{C}orr\\[X_{1,\\Omega_0}, X_{2,\\Omega_1}\\]",
    x_var_name="X_{1,\\Omega_0}",
    y_var_name="X_{2,\\Omega_1}",
    add_corr_definition_at_title=FALSE,
    use_dbl_line_title=TRUE
  )
  plot_C_1B_2A <- dep_coherence_lagged_plot(
    X_1B, X_2A,
    max_lags=max_lags,
    label_id=14,
    xlabel="$X_1$", ylabel="$X_2$", 
    color="#e4a027",
    corr_var_name="\\\\textbf{C}orr\\[\\cdot, \\cdot\\]",
    #corr_var_name="\\\\textbf{C}orr\\[X_{1,\\Omega_1}, X_{2,\\Omega_0}\\]",
    x_var_name="X_{1,\\Omega_1}",
    y_var_name="X_{2,\\Omega_0}",
    add_corr_definition_at_title=FALSE,
    use_dbl_line_title=TRUE
  )
  plot_C_1B_2B <- dep_coherence_lagged_plot(
    X_1B, X_2B,
    max_lags=max_lags,
    label_id=13,
    xlabel="$X_1$", ylabel="$X_2$", 
    color="#e4a027",
    corr_var_name="\\\\textbf{C}orr\\[\\cdot, \\cdot\\]",
    #corr_var_name="\\\\textbf{C}orr\\[X_{1,\\Omega_1}, X_{2,\\Omega_1}\\]",
    x_var_name="X_{1,\\Omega_1}",
    y_var_name="X_{2,\\Omega_1}",
    add_corr_definition_at_title=FALSE,
    use_dbl_line_title=TRUE
  ) +
  theme(plot.background = element_rect(fill = "#f6deb4", colour="transparent"))
  
  
  gridExtra::grid.arrange(grobs=list(
      plot_Z1,
      plot_Z2,
      plot_X1,
      plot_X2,
      plot_R_12,
      
      plot_R_1A_2A,
      plot_R_1B_2B,
      plot_R_1A_2B,
      plot_R_1B_2A,
     
     
      plot_C_12,
      plot_C_1A_2A,
      plot_C_1B_2B,
      plot_C_1A_2B,
      plot_C_1B_2A
      
    ), ncol=6,
    layout_matrix=rbind(
      c(1, 1, 1,   3, 3, 3),
      c(2, 2, 2,   4, 4, 4),
      c(5, 5,    6, 6, 7, 7),
      c(5, 5,    8, 8, 9, 9),
      c(10,10,  11,11,12,12),
      c(10,10,  13,13,14,14)
    )
  )
}
