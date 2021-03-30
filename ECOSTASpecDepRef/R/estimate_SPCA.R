VALID_WINDOWS <- c("bartlett", "blackman", "boxcar", "gausswin", "hamming", "hanning", "triang")

spectrum_model <- function(data, estimation_type="FFT", var_order, number_freq_points=200, variance_window_C_T=11, n_trials=10, win_spans=seq(2, 20, by=2), win="gausswin"){
    freqs <- (0:(number_freq_points-1)) / ( 2 *number_freq_points)
	if(estimation_type %in% c(NULL, "olsVARS", "LassoBigVAR", "RawLSE", "LASSLE")){
		model <- var_coefficients(data, var_order, estimation_type)
	}else if(estimation_type == "GSE"){
		if(!(win %in% VALID_WINDOWS)) print(paste("Invalid window", win))
		model <- general_shrinkage_estimator(data, var_lag=var_order, C_T=variance_window_C_T, n_trials=n_trials, win_spans=win_spans, win=match.fun(win))
	}else if(estimation_type == "FFT"){
		if(!(win %in% VALID_WINDOWS)) print(paste("Invalid window", win))
		model <- non_parametric_spectrum_FFT(data, win_span=win_spans[1], win=match.fun(win))
	}else if(estimation_type == "PURE"){
		if(!(win %in% VALID_WINDOWS)) print(paste("Invalid window", win))
		model <- non_parametric_spectrum_PURE(data, n_trials=n_trials, win_spans=win_spans, win=match.fun(win))
	}else{
        stop(paste("Estimation method:", estimation_type, "not defined!"))
    }
	model
}

#' @export
estimate_SPCA <- function(model, X, n_components=2){
	d <- apply(X, 2, fft)
	#
	n_freqs <- length(model$freqs)
  # print(sprintf("n_freqs=%d", n_freqs))
	n_positive_freqs <- floor(n_freqs/2) + 1
	spectrum_nrow <- ncol(X)#ceiling(sqrt((spectrum_coeff)))
	explained_variance <- matrix(rep(0, n_freqs * n_components), nrow=n_freqs)
	lambda_factors <- matrix(rep(0, n_freqs * n_components), nrow=n_freqs)
	unsorted_lambda_factors <- matrix(rep(0, n_freqs * n_components), nrow=n_freqs)
	C_factors <- list()
	for(index_factor in 1:n_components){
		C_factors[[index_factor]] <- matrix(rep(0, spectrum_nrow * n_freqs), nrow=n_freqs)
	}
  #print(dim(model$spectrum)); readline("[model$spectrum]")
	for(index in 1:n_freqs){
		negative_index <- if (index <= 1 ) -1 else (n_freqs - index + 2) # index should be > 1
		# Apply PCA on each frequency
		freq <- model$freqs[index]
		spectrum_coeff <- model$spectrum[index,]
		spectrum_coeff <- (matrix(spectrum_coeff, nrow=spectrum_nrow))
    #print(spectrum_coeff)
		spectrum_eig <- eigen(spectrum_coeff)
		unsorted_spectrum_eig <- spectrum_eig
    ### print(index); print(spectrum_eig);
		ordered_eigvals <- sort(spectrum_eig$values, decreasing=TRUE, index.return=TRUE)
		spectrum_eig$values <- ordered_eigvals$x
		spectrum_eig$vectors <- spectrum_eig$vectors[, ordered_eigvals$ix]
		explained_variance[index,] <- cumsum(spectrum_eig$values[1:n_components]) / sum(spectrum_eig$values)
		lambda_factors[index,] <- spectrum_eig$values
		unsorted_lambda_factors[index,] <- unsorted_spectrum_eig$values
		
    for(index_factor in 1:n_components){
			C_factors[[index_factor]][index, ] <- spectrum_eig$vectors[, index_factor]
			if(negative_index > 0){ # Get the equivalent in the negative domain
				C_factors[[index_factor]][negative_index, ] <- Conj(spectrum_eig$vectors[, index_factor])
			}
		}
		if(negative_index > 0){ # Get the equivalent in the negative domain
			explained_variance[negative_index,] <- explained_variance[index,]
			lambda_factors[negative_index,] <- lambda_factors[index,]
		}
    #print(explained_variance[index,])
    #print(lambda_factors[index,])
    #print(" .      .")
	}

	#Obtain factors:
	time_factors <- matrix(rep(0, n_freqs * n_components), nrow=n_freqs)
	for(index_factor in 1:n_components){
		time_factors[, index_factor] <- Re(fft(rowSums(C_factors[[index_factor]] * d), inverse=TRUE)) / n_freqs
	}
  print("factors:")
  print(dim(time_factors))
  print("explained_variance:")
  print(dim(explained_variance))
  print("lambda_factors:")
  print(dim(lambda_factors))
  print("C_factors:")
  print(length(C_factors))
  print(dim(C_factors[[1]]))
  print("===")
	list(
		freqs=model$freqs,
		spectrum=model$spectrum,
		factors=time_factors,
		explained_variance=explained_variance,
		unsorted_lambda_factors=unsorted_lambda_factors,
		lambda_factors=lambda_factors,
		C_factors=C_factors
	)
}


# X: rows x cols dimmension
raw_periodogram <- function(X){
	d <- apply(X, 2, fft)
	#I <- (d %*% Conj(t(d))) / nrow(d)
	# I <- matrix(apply(d, 1, function(y) abs(y %*% Conj(t(y)))), ncol=ncol(d) ^ 2) / nrow(d)
	I <- t( matrix(c(apply(d, 1, function(y) abs(y %*% Conj(t(y))))), nrow=ncol(d)^2) ) / nrow(d) * sqrt(2*pi)
	I
}

#  > G = matrix(1:40, ncol=4)
#  > G
#        [,1] [,2] [,3] [,4]
#   [1,]    1   11   21   31
#   [2,]    2   12   22   32
#   [3,]    3   13   23   33
#   [4,]    4   14   24   34
#   [5,]    5   15   25   35
#   [6,]    6   16   26   36
#   [7,]    7   17   27   37
#   [8,]    8   18   28   38
#   [9,]    9   19   29   39
#  [10,]   10   20   30   40
#  > stats::filter(G, signal::triang(3), circular=TRUE)
#  Time Series:
#  Start = 1
#  End = 10
#  Frequency = 1
#     [,1] [,2] [,3] [,4]
#   1    7   27   47   67
#   2    4   24   44   64
#   3    6   26   46   66
#   4    8   28   48   68
#   5   10   30   50   70
#   6   12   32   52   72
#   7   14   34   54   74
#   8   16   36   56   76
#   9   18   38   58   78
#  10   15   35   55   75
#  > stats::filter(G[,1], signal::triang(3), circular=TRUE)
#   [1]  7  4  6  8 10 12 14 16 18 15
smoothed_periodogram <- function(X, win_span=nrow(X)%/%10, win=signal::hanning){
	T <- nrow(X)
	win <- win(win_span)
	win_factor <- sum(abs(win))
	#print(win); print(paste("W=", win_factor))
	I_w <- raw_periodogram(X)
	f_tilde <- matrix(stats::filter(I_w, win, circular=TRUE), ncol=ncol(I_w)) / win_factor
	#f_tilde <- matrix(stats::filter(I_w, win, circular=TRUE), nrow=nrow(I_w))
	f_tilde
}

default_frequency_range <- function(n_freqs){
	n_positive_freqs <- floor(n_freqs/2) + 1
	freqs <-  c(1:n_positive_freqs-1, (-n_positive_freqs - (n_freqs %% 2) + 2):-1) / n_freqs
	freqs
}

non_parametric_spectrum_FFT <- function(X, win_span=nrow(X)%/%10, win=signal::hanning){
	###
  f_star <- smoothed_periodogram(X, win_span=win_span, win=win)
	###f_star <- raw_periodogram(X)
  ###
	freqs <- default_frequency_range(nrow(f_star))
	spectr <- list(spectrum=f_star, freqs=freqs)
	class(spectr) <- "modelNonParametric"
	spectr
}

#' @export
SpectralPCA <- function(data, var_order, estimation_type="FFT", variance_window_C_T=1, n_trials=10, win_spans=seq(2, 20, by=2), win="gausswin"){
	model <- spectrum_model(data, estimation_type=estimation_type, variance_window_C_T=variance_window_C_T, n_trials=n_trials, win_spans=win_spans, win=win)
	n_components <- ncol(data)
	pca <- estimate_SPCA(model, X=data, n_components=n_components)
  n_channels <- ncol(data)
	marginal_spectrum <- (1:n_channels) + (1:n_channels - 1) * n_channels
	list(
    freqs=pca$freqs[1:ceiling(nrow(data)/2)],
    spectra=pca$spectrum[1:ceiling(nrow(data)/2), ][, marginal_spectrum],
    PCs=pca$factors,
    spca=pca
  )
}
