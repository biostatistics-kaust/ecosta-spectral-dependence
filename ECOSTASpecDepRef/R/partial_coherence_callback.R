#' @export
coherence_callback <- function(max_lags, sign_lags="-") {
  function(X, fs, selected_band) {
    N <- ncol(X)
    print(sprintf("  => %s :: %s %s", N, nrow(X), ncol(X) ))
    metric <- matrix(1, nrow=N, ncol=N)
    for(i in 1:N){
      for(j in 1:N){
        if(i == j) next;
        ux <- extract_rhythms(X[,i], fs=fs, rhythms=c(selected_band), print_filters=FALSE)[[selected_band]]
        uy <- extract_rhythms(X[,j], fs=fs, rhythms=c(selected_band), print_filters=FALSE)[[selected_band]]
        metric[i,j] <- coherence_corr(ux, uy, max_lags=max_lags, sign_lags=sign_lags)$coherence
        metric[j,i] <- metric[i,j]
      }
    }
    metric
  }
}

#' @export
partial_coherence_callback <- function(max_lags, sign_lags="-") {
  function(X, fs, selected_band) {
    N <- ncol(X)
    #stopifnot(N == 3)
    metric <- matrix(1, nrow=N, ncol=N)
    for(i in 1:N){
      for(j in 1:N){
        if(i == j) next;
        other <- setdiff(1:N, c(i, j))[1]
        ux0 <- extract_rhythms(X[,i], fs=fs, rhythms=c(selected_band), print_filters=FALSE)[[selected_band]]
        uy0 <- extract_rhythms(X[,j], fs=fs, rhythms=c(selected_band), print_filters=FALSE)[[selected_band]]
        uz0 <- extract_rhythms(X[,other], fs=fs, rhythms=c(selected_band), print_filters=FALSE)[[selected_band]]
        ux <- resid(lm(ux0 ~ 1 + uz0))
        uy <- resid(lm(uy0 ~ 1 + uz0))
        metric[i,j] <- coherence_corr(ux, uy, max_lags=max_lags, sign_lags=sign_lags)$coherence
        metric[j,i] <- metric[i,j]
      }
    }
    metric
  }
}
