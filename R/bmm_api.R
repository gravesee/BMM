#' Train Bernoulli Mixture Model
#' @param data A Matrix. Either dense or a sparse, column-compressed, pattern matrix
#' @param K The number of prototypes to train
#' @param max.iter Maximum number of iterations. The algorithm will quit if this is reached
#' before convergence
#' @param verbose Print convergence information while training
#' @param hbbmm Impose a beta prior on the prototype components. Calculated empircally.
#' @return A list with class attribute, \code{BMM}.
#' @export
setGeneric("BMM", function(data, K, max.iter=10L, verbose=TRUE, hbbmm=TRUE) standardGeneric("BMM"))

#' @export
setMethod(
  "BMM",
  c("matrix", "numeric", "numeric", "logical", "logical"),
  function(data, K, max.iter, verbose, hbbmm) {
    
    if (!identical(typeof(data), "integer")) mode(data) <- "integer"
    
    res <- .Call(
      C_bmm_dense_matrix,
      data,
      nrow(data),
      ncol(data),
      as.integer(K),
      as.integer(max.iter),
      as.integer(verbose),
      as.integer(hbbmm))
    
    res$prototypes <- t(res$prototypes)
    
    class(res) <- "BMM"
    
    res
    
  })


#' @export
setMethod(
  "BMM",
  c("ngCMatrix", "numeric", "numeric", "logical", "logical"),
  function(data, K, max.iter, verbose, hbbmm) {
    
    
    res <- .Call(
      C_bmm_sparse_matrix,
      data@p,
      data@i,
      nrow(data),
      ncol(data),
      as.integer(K),
      as.integer(max.iter),
      as.integer(verbose),
      as.integer(hbbmm))
    
    res$prototypes <- t(res$prototypes)
    
    class(res) <- "BMM"
    
    res
    
  })

setOldClass("BMM")

#' Predict Bernoulli Mixture Model on new data
#' @param object An object of class \code{"BMM"} created by the \code{\link{BMM}} method.
#' @param newdata A dense, integer matrix or an ngCMatrix
#' @return A list with two elemends: "z" a matrix of each prototype likelihood for each sample, and
#' "ll" the likelihood of the data currently hard-coded to 0.
#' @export
setMethod("predict", "BMM", function(object, newdata) {
  
  ## stop if columns don't match prototype dims
  
  cls <- class(newdata)
  
  switch(
    cls,
    "matrix" = {
      
      if (!identical(typeof(newdata), "integer")) mode(newdata) <- "integer"
      
      .Call(C_predict_dense_matrix,
            newdata,
            t(object$prototypes),
            object$pis,
            nrow(newdata),
            ncol(newdata),
            length(object$pis))
      
    },
    
    "ngCMatrix" = {
      
      .Call(C_predict_sparse_matrix,
            newdata@p,
            newdata@i,
            t(object$prototypes),
            object$pis,
            nrow(newdata),
            ncol(newdata),
            length(object$pis))
      
    },
    
    stop("No predict method matching newdata class: %s", cls)
  )
  
  
})