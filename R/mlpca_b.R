#' Maximum Likelihood Principal Component Analysis for Mode B Error
#' Conditions
#'
#' @param X MxN matrix of measurements
#' @param Xsd MxN matrix of measurements error standard deviations
#' @param p Rank of the model's subspace, p must be than the minimum of M and N
#'
#' @return The parameters returned are the results of SVD on the estimated
#' subspace. The quantity Ssq represents the sum of squares of weighted
#' residuals.
#'
#' @export
#'
#' @examples
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
