LASSLE <- function(X, y, keep.resid=FALSE){
  cvfit <- glmnet::cv.glmnet(X, y, nfolds=10)
  #cvfit <- glmnet::cv.glmnet(X, y, nfolds=4)
	nonzero_cols <- as.logical(coef(cvfit, s="lambda.min") != 0)[-c(1)]
	nonzero_cols <- which(nonzero_cols)
	fit_lse <- LSE(matrix(X[,nonzero_cols], ncol=length(nonzero_cols)), y, keep.resid=keep.resid)
    beta_intercept <- fit_lse$beta_intercept
    std_intercept <- fit_lse$std_intercept
    t_val_intercept <- fit_lse$t_val_intercept
    #
    beta <- matrix(rep(0, each=ncol(X)), nrow=1)
    colnames(beta) <- colnames(X)
    beta[nonzero_cols] <- fit_lse$beta
    #
    std <- matrix(rep(0, each=ncol(X)), nrow=1)
    colnames(std) <- colnames(X)
    std[nonzero_cols] <- fit_lse$std
    #
    t_val <- matrix(rep(0, each=ncol(X)), nrow=1)
    colnames(t_val) <- colnames(X)
    t_val[nonzero_cols] <- fit_lse$t_val
    #
    list(
        beta_intercept=beta_intercept,
        std_intercept=std_intercept,
        t_val_intercept=t_val_intercept,
        beta=beta,
        std=std,
        t_val=t_val,
        resid=fit_lse$resid
    )
}
