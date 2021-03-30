var_coefficients <- function(data, var_order, var_type=NULL){
  T <- nrow(data)
  k <- ncol(data)
  if(is.null(colnames(data))) {
    colnames(data) <- paste0("X", 1:k)
  } 
  if(is.null(var_type) || var_type == "olsVARS"){
    model_results <- vars::VAR(
      data, 
      p=var_order, 
      type="const"
    )
    Phi <- Bcoef(model_results)[,1:(k * var_order)]
    intercepts <- Bcoef(model_results)[, (k * var_order + 1)]
    Sigma <- summary(model_results)$covres
  }else if(var_type == "LassoBigVAR"){
    model <- BigVAR::constructModel(
      as.matrix(data),
      p=var_order,
      struct="Basic",
      gran=c(150,10),
      RVAR=FALSE,
      h=1,
      cv="Rolling",
      MN=FALSE,
      verbose=FALSE,
      IC=TRUE)
    model_results <- BigVAR::cv.BigVAR(model)
    intercepts <- model_results@betaPred[, 1]
    Phi <- model_results@betaPred[, 2:(k * var_order + 1)]
    Sigma <- var(model_results@resids)
    #Bias?# Sigma <- Sigma * (T - k * var_order - 1)  / (T + var_order - 1)
  }else if(var_type == "RawLSE"){
    Z <- expand_lag(as.matrix(data), lags=var_order)
    fitted <- VAR(Z$X, Z$Y, LS.estimator=LSE)
    intercepts <- fitted$beta_intercept
    Phi <- fitted$beta
    Sigma <- fitted$sigma
  }else if(var_type == "LASSLE"){
    Z <- expand_lag(as.matrix(data), lags=var_order)
    fitted <- VAR(Z$X, Z$Y, LS.estimator=LASSLE)
    intercepts <- fitted$beta_intercept
    Phi <- fitted$beta
    Sigma <- fitted$sigma
  }else{
    stop(paste("VAR estimation method:", var_type, "not defined!"))
  }
  list(
    intercepts=intercepts,
    #Phi=t(Phi),
    #Sigma=t(Sigma)
    Phi=Phi,
    Sigma=Sigma
  )
}

apply_lag <- function(X, l, total_lags){
    Y <- X[-c(0:(total_lags-l), (-l):1 + 1 + nrow(X)), ]
    colnames(Y) <- paste0(colnames(X), ".", l)
    Y
}

expand_lag <- function(X, lags){
    t.X <- matrix(sapply(1:lags, function(l) apply_lag(X, l, lags)), ncol=lags*ncol(X))
    colnames(t.X) <- c(sapply(1:lags, function(l) paste0(colnames(X), ".L", l)))
    list(
        X=t.X,
        Y=apply_lag(X, 0, lags)
    )
}

#Z=expand_lag(XX, lags=2)
LSE <- function(X, Y, keep.resid=FALSE){
    M <- cbind(1, X)
    tM_M_inv <- solve(t(M) %*% M)
    H <- tM_M_inv %*% t(M)
    beta_hat <- H %*% Y
    z.err = (Y - M %*% H %*% Y)
    sse_err <- sapply(1:ncol(z.err), function(i) var(z.err[,i]))
    sse_err <- sse_err * (nrow(X)+1)/(nrow(X)-ncol(X)+1)
    beta_std <- sqrt(matrix(diag(tM_M_inv), ncol=1) %*% sse_err)
    colnames(beta_std) <- colnames(beta_hat)
    rownames(beta_std) <- rownames(beta_hat)
    t_val <- beta_hat / beta_std
    m1 <- nrow(beta_hat)
    #m0 <- (if(keep.intercept) 1 else 2)
    list(
        beta_intercept=beta_hat[1],
        std_intercept=beta_std[1],
        t_val_intercept=t_val[1],
        beta=t(beta_hat[2:m1,]),
        std=t(beta_std[2:m1,]),
        t_val=t(t_val[2:m1,]),
        resid=if(keep.resid) z.err else NULL
    )
}


VAR <- function(X, Y, LS.estimator=LSE) {
    copy_y_names <- function(V){
        rownames(V) <- colnames(Y)
        V
    }
    Omega <- t(sapply(1:ncol(Y), function(k) LS.estimator(X, Y[,k], keep.resid=TRUE)))
    sigma <- cov(do.call(cbind, Omega[, "resid"]))
    rownames(sigma) <- colnames(sigma) <- colnames(Y)
    k_resid = which.max(colnames(Omega) == "resid")
    Omega <- Omega[, -c(k_resid)]
    VAR.Omega <- lapply(colnames(Omega), function(p) copy_y_names(do.call(rbind, Omega[,p])))
    names(VAR.Omega) <- colnames(Omega)
    VAR.Omega$sigma <- sigma
    VAR.Omega
}
