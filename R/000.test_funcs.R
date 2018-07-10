
#' @export
test_dense_matrix <- function(m, n, d) {
    .Call(C_test_dense_matrix, m, n, d)
}