options(encoding = "UTF-8")

#' @title Maximum Likelihood Principal Component Analysis for Mode D Error
#' Conditions
#'
#' @description Performs maximum likelihood principal components analysis for
#' mode D error conditions (commom row covariance matrices).
#' Employs rotation and scaling of the original data.
#'
#' @author Renan Santos Barbosa
#'
#' @param X IxJ matrix of measurements
#' @param Cov JxJ matrix of measurement error covariance, which is commom to all
#' rows
#' @param p Rank of the model's subspace
#'
#' @return  The parameters returned are the results of SVD on the estimated
#' subspace. The quantity Ssq represents the sum of squares of weighted
#' residuals.
#'
#' @details The returned parameters, U, S and V, are analogs to the
#' truncated SVD solution, but have somewhat different properties since they
#' represent the MLPCA solution. In particular, the solutions for different
#' values of p are not necessarily nested (the rank 1 solution may not be in the
#' space of the rank 2 solution) and the eigenvectors do not necessarily account
#' for decreasing amounts of variance, since MLPCA is a subspace modeling
#' technique and not a variance modeling technique.
#'
#'
#' @references Wentzell, Peter D., et al. "Maximum likelihood principal
#' component analysis." Journal of Chemometrics: A Journal of the Chemometrics
#' Society 11.4 (1997): 339-366.
#'
#' @exportClass list
#' @import Matrix
#' @export
#' @encoding UTF-8
#'
#' @example set.seed(123)
#'
#' X = matrix(c(1,2,3,4,
#' 5,6,7,8,
#' 3,2,1,9),
#' nrow =3,
#' ncol=4,
#' byrow = TRUE)
#'
#' Cov = matrix(c(1,0.2,0.3,0.4,
#'                0.6,0.34,0.123,0.112,
#'                0.145,0.3451,233,0.651,
#'                0.2341,0.2341,0.5341,0.3544),
#'              nrow = 4,
#'              ncol = 4,
#'              byrow = TRUE
#' )
#'
#'
#' p = 1
#'
#' res <- mlpca_d(X,Cov,p)
#'
#' xhat <- res$U%*%res$S%*%t(res$V)
#'

mlpca_d <- function(X, Cov, p) {
  m <- dim(x = X)[1]
  n <- dim(x = X)[2]

  Df <- (m - p) * (n - p)

  DecomCov <- svd(Cov)

  S1 <- diag(DecomCov$d,
    nrow = nrow(Cov),
    ncol = nrow(Cov)
  )

  U1 <- DecomCov$u

  V1 <- DecomCov$v

  CovRank <- Matrix::rankMatrix(Cov)
  if (CovRank < n) { # If the covariance matrix is singular, then
    # expand uniformly in deficient directions

    Scale <- c(
      sqrt(diag(matrix(S1[1:CovRank, 1:CovRank], CovRank, CovRank))),
      (matrix(1, 1, n - CovRank) *
        sqrt(S1[CovRank, CovRank]) * 0.01)
    )
  } else {
    Scale <- sqrt(diag(S1))
  }

  Z <- X %*% U1 %*% diag(1 / Scale)

  decomZ <- svd(Z)

  S2 <- diag(decomZ$d,
    nrow = nrow(Z),
    ncol = nrow(Z)
  )
  U2 <- decomZ$u

  V2 <- decomZ$v

  ZCalc <- U2[, 1:p] %*% matrix(S2[1:p, 1:p], p, p) %*% t(V2[, 1:p])

  Ssq <- 0

  for (i in 1:m) {
    Ssq <- Ssq + norm(matrix(ZCalc[i, 1:CovRank],
      nrow = 1,
      ncol = CovRank
    ) -
      matrix(Z[i, 1:CovRank],
        nrow = 1,
        ncol = CovRank
      ),
    type = "2"
    )^2
  }

  XCalc <- ZCalc %*% diag(Scale) %*% t(U1)


  DecomXCalc <- svd(XCalc)

  S <- diag(DecomXCalc$d,
    nrow = nrow(XCalc),
    ncol = nrow(XCalc)
  )

  U <- DecomXCalc$u

  V <- DecomXCalc$v

  S <- diag(S[1:p, 1:p],
    ncol = p,
    nrow = p
  )

  U <- matrix(U[, 1:p],
    nrow = m,
    ncol = p
  )

  V <- matrix(V[, 1:p],
    nrow = n,
    ncol = p
  )


  result <- list(
    "U" = U,
    "S" = S,
    "V" = V,
    "Ssq" = Ssq
  )

  return(result)
}
