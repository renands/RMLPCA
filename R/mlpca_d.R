#' Maximum likelihood principal component analysis for mode D error conditions
#'
#' @description Performs maximum likelihood principal components analysis for
#'   mode D error conditions (commom row covariance matrices).
#'   Employs rotation and scaling of the original data.
#'
#' @param X IxJ matrix of measurements
#' @param Cov JxJ matrix of measurement error covariance, which is commom to all
#'   rows
#' @param p Rank of the model's subspace
#'
#' @return  The parameters returned are the results of SVD on the estimated
#'   subspace. The quantity Ssq represents the sum of squares of weighted
#'   residuals.
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
#'   library(RMLPCA)
#'   data(data_clean)
#'   data(data_error_d)
#'   # covariance matrix
#'   data(cov_d)
#'   data(data_cleaned_mlpca_d)
#'   # data that you will usually have on hands
#'   data_noisy <- data_clean + data_error_d
#'
#'   # run mlpca_c with rank p = 5
#'   results <- RMLPCA::mlpca_d(
#'     X = data_noisy,
#'     Cov = cov_d,
#'     p = 2
#'   )
#'
#'   # estimated clean dataset
#'   data_cleaned_mlpca <- results$U %*% results$S %*% t(results$V)
mlpca_d <- function(X, Cov, p) {
  m <- base::dim(x = X)[1]
  n <- base::dim(x = X)[2]

  Df <- (m - p) * (n - p)

  DecomCov <- base::svd(Cov)

  S1 <- base::diag(DecomCov$d,
    nrow = base::nrow(Cov),
    ncol = base::nrow(Cov)
  )

  U1 <- DecomCov$u

  V1 <- DecomCov$v

  CovRank <- Matrix::rankMatrix(Cov)
  if (CovRank < n) { # If the covariance matrix is singular, then
    # expand uniformly in deficient directions

    Scale <- base::c(
      base::sqrt(
        base::diag(
          base::matrix(
            S1[1:CovRank, 1:CovRank], CovRank, CovRank
          )
        )
      ),
      (base::matrix(1, 1, n - CovRank) *
        base::sqrt(S1[CovRank, CovRank]) * 0.01)
    )
  } else {
    Scale <- base::sqrt(base::diag(S1))
  }

  Z <- X %*% U1 %*% base::diag(1 / Scale)

  decomZ <- base::svd(Z)

  S2 <- base::diag(decomZ$d,
    nrow = base::nrow(Z),
    ncol = base::nrow(Z)
  )
  U2 <- decomZ$u

  V2 <- decomZ$v

  ZCalc <- U2[, 1:p] %*% base::matrix(S2[1:p, 1:p], p, p) %*% t(V2[, 1:p])

  Ssq <- 0

  for (i in 1:m) {
    Ssq <- Ssq + base::norm(base::matrix(ZCalc[i, 1:CovRank],
      nrow = 1,
      ncol = CovRank
    ) -
      base::matrix(Z[i, 1:CovRank],
        nrow = 1,
        ncol = CovRank
      ),
    type = "2"
    )^2
  }

  XCalc <- ZCalc %*% base::diag(Scale) %*% base::t(U1)


  DecomXCalc <- base::svd(XCalc)

  S <- base::diag(DecomXCalc$d,
    nrow = base::nrow(XCalc),
    ncol = base::nrow(XCalc)
  )

  U <- DecomXCalc$u

  V <- DecomXCalc$v

  if (p == 1) {
    S <- base::diag(S[1:p, 1:p],
      ncol = p,
      nrow = p
    )
  } else {
    S <- S[1:p, 1:p]
  }


  U <- base::matrix(U[, 1:p],
    nrow = m,
    ncol = p
  )

  V <- base::matrix(V[, 1:p],
    nrow = n,
    ncol = p
  )


  result <- base::list(
    "U" = U,
    "S" = S,
    "V" = V,
    "Ssq" = Ssq
  )

  return(result)
}
