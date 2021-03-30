#https://github.com/ellchow/moRe/blob/master/mdls.R
#' @export
PCA <- function(X, center=T, scale=T, is.symm=FALSE){
  ## perform principal component analysis on x using spectral decomposition
  ## X: matrix-like on which to perform PCA
  ## center: center the columns
  ## scale: scale the columns

  ## head(spectral.pca.predict(spectral.pca.fit(as.matrix(USArrests)),USArrests))

  X <- scale(X,scale=scale, center=center)
  centers <- attr(X, "scaled:center")
  scales <- attr(X, "scaled:scale")
  ###stop.if(any(scales == 0), 'cannot rescale constant columns to unit variance (%s)', paste(which(scales == 0),collapse=','))

  if(!is.symm){
    XtX <- t(X) %*% X
    eig <- eigen(XtX, TRUE)
  }else
    eig <- eigen(X, TRUE)

  P <- t(eig$vectors)
  Xpca <- t(P %*% t(X))

  ordered_eigvals <- sort(eig$values, decreasing=TRUE, index.return=TRUE)
  spectrum_eig <- list()
  spectrum_eig$values <- ordered_eigvals$x
  
  spectrum_eig$vectors <- eig$vectors[, ordered_eigvals$ix]
  explained_variance <- cumsum(spectrum_eig$values) / sum(spectrum_eig$values)
  
  print("PCA::eig_values")
  print(eig$values)
  print("PCA::eig_vectors")
  print(eig$vectors)
  print("PCA::explained_variance")
  print(explained_variance)
  
  list(
    p=P, 
    PCs=Xpca,
    #centers=if(is.null(centers)) NA else centers, 
    #scales=if(is.null(scales)) NA else scales,
    eig_values=spectrum_eig$values,
    eig_vectors=spectrum_eig$vectors,
    explained_variance=explained_variance
  )
}
