#' @export
setGeneric("BMM", function(data, K, max.iter=10L, verbose=1L, hbbmm=1L) standardGeneric("BMM"))

#' @export
setMethod(
  "BMM",
  c("matrix", "integer", "integer", "integer", "integer"),
  function(data, K, max.iter, verbose, hbbmm) {
    
    res <- .Call(C_bmm_dense_matrix, data, nrow(data), ncol(data),
          as.integer(K), as.integer(max.iter), as.integer(verbose), as.integer(hbbmm))
    
    res$prototypes <- t(res$prototypes)
    
    res
    
  })


#' @export
setMethod(
  "BMM",
  c("ngCMatrix", "integer", "integer", "integer", "integer"),
  function(data, K, max.iter, verbose, hbbmm) {
    
    
    res <- .Call(C_bmm_sparse_matrix, data@p, data@i, nrow(data), ncol(data),
                 as.integer(K), as.integer(max.iter), as.integer(verbose), as.integer(hbbmm))
    
    res$prototypes <- t(res$prototypes)
    
    res
    
  })


predict_dense_matrix <- function(data, protos, pis) {
  
  ## check everything makes sense here
  res <- .Call(C_predict_dense_matrix, data, t(protos), pis, nrow(data), ncol(data), length(pis))
  
  res
  
}


predict_sparse_matrix <- function(data, protos, pis) {
  
  ## check everything makes sense here
  res <- .Call(C_predict_sparse_matrix, data@p, data@i, t(protos), pis, nrow(data), ncol(data), length(pis))
  
  res
  
}

