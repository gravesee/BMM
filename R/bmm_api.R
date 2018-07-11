#' @export
setGeneric("BMM", function(data, K, max.iter=10L, verbose=1L) standardGeneric("BMM"))

#' @export
setMethod(
  "BMM",
  c("matrix", "integer", "integer", "integer"),
  function(data, K, max.iter, verbose) {
    
    res <- .Call(C_bmm_dense_matrix, data, nrow(data), ncol(data),
          as.integer(K), as.integer(max.iter), as.integer(verbose))
    
    res$prototypes <- t(res$prototypes)
    
    res
    
  })


#' @export
setMethod(
  "BMM",
  c("ngCMatrix", "integer", "integer", "integer"),
  function(data, K, max.iter, verbose) {
    
    
    res <- .Call(C_bmm_sparse_matrix, data@p, data@i, nrow(data), ncol(data),
                 as.integer(K), as.integer(max.iter), as.integer(verbose))
    
    res$prototypes <- t(res$prototypes)
    
    res
    
  })
