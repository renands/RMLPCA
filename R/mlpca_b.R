#' Maximum likelihood principal component analysis for mode B error conditions
#'
#' @description Performs maximum likelihood principal components analysis for
#'   mode B error conditions (independent errors, homoscedastic within a column).
#'   Equivalent to perfoming PCA on data scaled by the error SD, but results are
#'   rescaled to the original space.
#'
#' @param X MxN matrix of measurements.
#' @param Xsd MxN matrix of measurements error standard deviations.
#' @param p Rank of the model's subspace, p must be than the minimum of M and N.
#'
#' @return The parameters returned are the results of SVD on the estimated
#'   subspace. The quantity Ssq represents the sum of squares of weighted
#'   residuals. All the results are nested in a list format.
#'
#' @details The returned parameters, U, S and V, are analogs to the
#'   truncated SVD solution, but have somewhat different properties since they
#'   represent the MLPCA solution. In particular, the solutions for different
#'   values of p are not necessarily nested (the rank 1 solution may not be in
#'   the space of the rank 2 solution) and the eigenvectors do not necessarily
#'   account for decreasing amounts of variance, since MLPCA is a subspace
#'   modeling technique and not a variance modeling technique.
#'
#' @references Wentzell, P. D.
#'   "Other topics in soft-modeling: maximum likelihood-based soft-modeling
#'   methods." (2009): 507-558.
#'
#' @export
#'
#' @examples
#'
#' library(RMLPCA)
#' data(data_clean)
#' data(data_error_b)
#' data(sds_b)
#'
#' # data that you will usually have on hands
#' data_noisy <- data_clean + data_error_b
#'
#' # run mlpca_b with rank p = 2
#' results <- RMLPCA::mlpca_b(
#'   X = data_noisy,
#'   Xsd = sds_b,
#'   p = 2
#' )
#'
#' # estimated clean dataset
#' data_cleaned_mlpca <- results$U %*% results$S %*% t(results$V)
mlpca_b <- function(X, Xsd, p) {
  m <- base::dim(x = X)[1]
  n <- base::dim(x = X)[2]

  if (p > min(m, n)) {
    stop("mlpca_b:err1 - Invalid rank for MLPCA decomposition")
  }

  ml <- base::dim(x = Xsd)[1]
  nl <- base::dim(x = Xsd)[2]

  if (m != ml | n != nl) {
    stop("mlpca_b:err2 - Dimensions of data and standard deviations do not match")
  }


  if (base::isFALSE(base::all(Xsd > 0))) {
    stop("mlpca_b:err3 - Standard deviations must be positive")
  }

  if (base::isTRUE(base::any(Xsd == 0))) {
    stop("mlpca_b:err4 - Zero value(s) for standard deviations")
  }

  SclMat <- Xsd

  Xsc <- X / SclMat

  DecomXsc <- RSpectra::svds(Xsc, p) # Decompose adjusted matrix
  U <- DecomXsc$u
  S <- base::diag(DecomXsc$d,
    nrow = base::length(DecomXsc$d),
    ncol = base::length(DecomXsc$d)
  )
  V <- DecomXsc$v

  XCalc <- (U %*% S %*% t(V)) * SclMat

  DecomXCalc <- RSpectra::svds(XCalc, p) # Decompose adjusted matrix
  U <- DecomXCalc$u
  S <- base::diag(DecomXCalc$d,
    nrow = base::length(DecomXCalc$d),
    ncol = base::length(DecomXCalc$d)
  )
  V <- DecomXCalc$v

  Ssq <- base::sum(base::sum(((X - XCalc) / SclMat)^2))

  result <- base::list(
    "U" = U,
    "S" = S,
    "V" = V,
    "Ssq" = Ssq
  )

  return(result)
}
