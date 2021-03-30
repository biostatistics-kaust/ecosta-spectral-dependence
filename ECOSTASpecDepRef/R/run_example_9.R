eig_base_colors <- function() c("#194a8d", "#8d4a19", "#df1c44", "#39a275", "#434468", "#901e06", "#19818d" )

plot_pca_vs_spca <- function(pca, spca, n_components, fs) {
  colors <- eig_base_colors()
  internal_plot <- function(y, pc_id, label_id, label) {
    slabel_id <- intToUtf8(utf8ToInt("A") + label_id - 1)
    #color1 <- colors[label_id %% (length(colors) - 1) + 1]
    #color2 <- colors[(label_id + 1) %% (length(colors) - 1) + 1]
    color1 <- colors[(label_id %/% 2) %% (length(colors) - 1) + 1]
    color2 <- colors[(label_id %/% 2) %% (length(colors) - 1) + 1]
    t <- (1:length(y)) / fs
    plots[[plot_id]] <<- qplot(t, y, geom="line", color=I(color1)) +
      labs(x="time [sec]", y="amplitude",
           title=sprintf("%s) %s - PC%d", slabel_id, label, pc_id)) + 
      csda_theme()
    
    slabel_id <- intToUtf8(utf8ToInt("A") + label_id)
    Sx <- estimate_smooth_spectrum(y, fcut=0.25) / sqrt(length(t))
    fx <- (1:length(Sx)) / (2 * length(Sx)) * fs
    plots[[plot_id + 1]] <<- qplot(fx, Sx, geom="line", color=I(color2)) +
      labs(x="frequency [Hz]", y="amplitude",
        title=sprintf("%s) Spectrum PC%d", slabel_id, pc_id)) + 
      csda_theme()
  }
  pca_pcs <- pca$PCs
  spca_pcs <- spca$PCs
  order <- matrix(NA, ncol=4, nrow=n_components)
  plots <- list()
  plot_id <- 1
  plot_ids <- c()
  for(k in 1:n_components) {
    internal_plot(y=pca_pcs[, k], pc_id=k, label_id=plot_id, label="PCA")
    plot_ids <- c(plot_ids, c(plot_id, plot_id, plot_id+1))
    plot_id <- plot_id + 2
    internal_plot(y=spca_pcs[, k], pc_id=k, label_id=plot_id, label="SPCA")
    plot_ids <- c(plot_ids, c(plot_id, plot_id, plot_id+1))
    plot_id <- plot_id + 2
  }
  #print(t(matrix(plot_ids, nrow=6)))
  q <- gridExtra::grid.arrange(grobs=plots,
    ncol=6,
    layout_matrix=t(matrix(plot_ids, nrow=6))
  )
  q
}


plot_vector_eig2d <- function(eigvecs, eigvals) {
  repsilon <- function(x) (0.10 * (max(x) - min(x))) 
  rmin <- function(x, dx=0.1) (min(x) - 0.3 * (max(x) - min(x))) 
  rmax <- function(x, dx=0.1) (max(x) + 0.3 * (max(x) - min(x))) 
  
  colors <- eig_base_colors()
  
  g <- ggplot(as.data.frame(eigvecs)) + 
    geom_vline(xintercept=0, colour="black", linetype="dotted") + 
    geom_hline(yintercept=0, colour="black", linetype="dotted") 
  for(i in 1:nrow(eigvecs)) {
    efactor <- sqrt(eigvals[i])
    #print(i); print(eigvecs[i,])
    #print(i); print(eigvals[i])
    text <- TeX(sprintf("$X_{%d}$", i))
    x <- eigvecs[i, 1] * efactor
    y <- eigvecs[i, 2] * efactor
    #print(sprintf("$X_{%d}$ %f %f", i, x, y))
    g <- g + geom_segment(
      x=0, y=0,
      xend=x,
      yend=y,
      color=colors[i],
      arrow=arrow(
        length=unit(0.03, "npc"),
        ends="last"
        #ends="both"
        #ends="first"
    )) + geom_text(
      x=(x + (if(x < 0) -1 else +1) * repsilon(eigvecs[, 1]*efactor)),
      y=(y + (if(y < 0) -1 else +1) * repsilon(eigvecs[, 2]*efactor) * 0.1),
      #y=y,
      color=colors[i],
      label=text,
      vjust="inward",
      hjust="inward")
  }
  g <- g + 
    csda_theme() +
    labs(
      x="PC1",
      y="PC2"
    ) +
    coord_cartesian(
      xlim=c(rmin(eigvecs[,1] * sqrt(eigvals)), rmax(eigvecs[,1] * sqrt(eigvals))),
      ylim=c(rmin(eigvecs[,2] * sqrt(eigvals)), rmax(eigvecs[,2] * sqrt(eigvals)))
    )
  g
}

plot_eigvec_pca_vs_spca <- function(pca, spca, fs, eeg_bands=c("all")) {
  slabel_id <- function(label_id) intToUtf8(utf8ToInt("A") + label_id - 1)
  #
  plots <- list()
  ncolumns <- dim(spca$C_factors[[1]])[2]
  #print(spca$spca$C_factors); print(ncolumns)
  vectPCA <- pca$eig_vectors[, 1:2]
  vectSPCA <- list()
  valsSPCA <- list()
  eeg_rhythms_vals <- eeg_rhythms()
  for(i in 1:length(eeg_bands)) {
    freq1 <- eeg_rhythms_vals[[eeg_bands[i]]][1]
    freq2 <- eeg_rhythms_vals[[eeg_bands[i]]][2]
    N <- nrow(spca$C_factors[[1]])
    stopifnot(freq1 / fs < 0.5)
    stopifnot(freq2 / fs < 0.5)
    k1 <- ceiling(N * freq1 / fs)
    k2 <- ceiling(N * freq2 / fs)
    vectY <- matrix(NA, nrow=ncolumns, ncol=2)
    valsY <- rep(0, ncolumns)
    for(k in 1:ncolumns) {
      vectY[k, 1:2]  <- colMeans(spca$C_factors[[k]][k1:k2,1:2])
      valsY[k]  <- mean(spca$lambda_factors[k1:k2])
    }
    vectSPCA[[i]] <- vectY
    valsSPCA[[i]] <- valsY
  }
  pca_eigvecs <- plot_vector_eig2d(vectPCA, pca$eig_values) +
    labs(title=sprintf("%s) PCA", slabel_id(1)))
  #plot(pca_eigvecs); readline("...")
  plots[[1]] <- pca_eigvecs
  #
  for(i in 1:length(eeg_bands)) {
    spca_eigvecs <- plot_vector_eig2d(vectSPCA[[i]], valsSPCA[[i]]) +
      labs(title=sprintf("%s) SPCA (%s band)", slabel_id(i+1), eeg_bands[i]))
    #plot(spca_eigvecs); readline("...")
    plots[[i+1]] <- spca_eigvecs
  }
  #print(eeg_bands); print(length(plots)); print(t(matrix(1:length(plots), nrow=2)))
  q <- gridExtra::grid.arrange(grobs=plots,
    ncol=2,
    layout_matrix=t(matrix(1:length(plots), nrow=2))
  )
  #plot(q)
  q
}

#' @export
run_example_9 <- function() {
  set.seed(1009)

  #x1 <- simulate_unstable_oscillator(N=N+mN, f=2, tau=4, fs=fs, sigma=1)
  #x2 <- simulate_unstable_oscillator(N=N+mN, f=15, tau=1, fs=fs, sigma=1)
  N <- 200
  x1 <- (1:N)
  x2 <- (1:N) ** 2
  x3 <- (1:N) ** 3
  #
  fs <- 128
  mixture_matrix <- rbind(
    c(1, 0, 1),
    c(0, 0, 1),
    c(1, 0, 1),
    c(0, 1, 0),
    c(1, 1, 0)
  )
  print(mixture_matrix)
  Z <- simulate_sparse_contemporaneous_mixture(
    N, 
    mixture_matrix, 
    fs, 
    lag_matrix=NULL, 
    f=c(2, 10, 35), 
    tau=c(1, 1, 1) * 3, 
    sigma=c(1, 1, 1)*2, 
    sigma_epsilon=c(1, 1, 1, 1, 1) 
  )
  X <- Z[[1]]
  Y <- Z[[2]]
  #
  Q1 <- PCA(t(Y))
  print(paste("explained_variance", Q1$explained_variance))
  print(paste("eig_values", Q1$eig_values))
  print(paste("eig_vectors", dim(Q1$eig_vectors)))
  Q2 <- SpectralPCA(
    t(Y),
    var_order=3,
    estimation_type="FFT"
  )
  q <- plot_pca_vs_spca(pca=Q1, spca=Q2, n_components=3, fs=fs)
  csda_save_image(q, "principal-components-pca-spca", height=15*1.2)
  
  q <- plot_eigvec_pca_vs_spca(pca=Q1, spca=Q2$spca, fs=fs, eeg_bands=c(
    "delta",
    "alpha",
    "gamma"
  ))
  csda_save_image(q, "eigvectors-pca-spca", height=12*1.5, width=18)
  
}



