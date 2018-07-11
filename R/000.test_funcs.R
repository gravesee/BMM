
#' @export
test_dense_matrix <- function(m, n, d) {
    .Call(C_test_dense_matrix, m, n, d)
}




#' @export
test_sparse_matrix <- function(m, n, d) {
  .Call(C_test_sparse_matrix, m@p, m@i, n, d)
}