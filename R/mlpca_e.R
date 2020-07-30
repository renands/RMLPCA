options(encoding = "UTF-8")

#' @title Maximum Likelihood Principal Component Analysis for Mode E Error
#' Conditions
#'
#' @description  Performs maximum likelihood principal components analysis for
#' mode E error conditions (correlated errors, with a different covariance
#' matrix for each row, but no error correlation between the rows). Employs an
#' ALS algorithm.
#'
#' @author Renan Santos Barbosa
#'
#' @param X IxJ matrix of measurements
#' @param Cov JXJXI matrices of measurement error covariance
#' @param p Rank of the model's subspace, p must be than the minimum of I and J
#'
#' @return  The parameters returned are the results of SVD on the estimated
#' subspace. The quantity Ssq represents the sum of squares of weighted
#' residuals. ErrFlag indicates the convergence condition,
#' with 0 indicating normal termination and 1 indicating the maximum number of
#' iterations have been exceeded.
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
#' @import RSpectra
#' @import Matrix
#'
#' @export
#' @encoding UTF-8
#'
#' @example
#' X <- matrix(c(1, 2, 3, 5, 6, 7), nrow = 2, ncol = 3, byrow = TRUE)
#' Cov <- array(c(c(1, 4, 7, 2, 5, 8, 3, 6, 9),
#' c(10, 13, 16, 11, 14, 17, 12, 15, 18)),
#'              dim = c(3, 3, 2)
#' )
#' p <- 1
#'
#' res <- mlpca_e(X = X, Cov = Cov, p = p)
#' res$U %*% res$S %*% t(res$V)
#'
mlpca_e <- function(X, Cov, p) {
  m <- dim(x = X)[1]

  n <- dim(x = X)[2]

  if (p > min(m, n)) {
    stop("mlpca_e:err1 - Invalid rank for MLPCA decomposition")
  }

  ml <- dim(x = Cov)[1]
  nl <- dim(x = Cov)[2]
  q <- dim(x = Cov)[3]


  if (n != ml | n != nl | m != q) {
    stop("mlpca_e:err2 - Invalid dimensions of covariance matrix")
  }

  # Initialization -------------------------------------------------------------


  ConvLim <- 1e-10 # Convergence Limit
  MaxIter <- 20000 # Maximum no. of iterations

  # Calculate the inverse of the full covariance matrix blockwise --------------

  mn <- m * n

  Q <- matrix(
    0,
    mn,
    mn
  )

  for (i in 1:m) {
    Indx <- (i - 1) * n
    Tmp <- Cov[, , i]
    Q[(Indx + 1):(Indx + n), (Indx + 1):(Indx + n)] <- pracma::pinv(Tmp)
  }

  # Now find the commutation matrix for the covariance matrix and
  # apply to Q

  ix <- matrix(1:mn)
  iy <- array(t(array(ix, dim = c(m, n))), c(mn, 1))

  K <- Matrix::sparseMatrix(
    i = ix,
    j = iy,
    x = 1,
    dims = c(
      mn,
      mn
    )
  )
  K <- as.matrix(K)

  Q <- (t(K) %*% Q %*% K)

  # Generate initial estimates assuming an average covariance matrix

  Covavg <- rowMeans(Cov, dims = 2)

  res_mlpca_d <- mlpca_d(
    X = X,
    Cov = Covavg,
    p
  )

  U <- res_mlpca_d$U
  S <- res_mlpca_d$S
  V <- res_mlpca_d$V

  # Main loop to do alternating regression for MLPCA solution

  Count <- 0
  Sold <- 0
  ErrFlag <- -1

  XX <- X

  while (ErrFlag < 0) {
    Count <- Count + 1

    # Vectorize X and create big U.  A Kronecker product is probably
    # prettier for generating Ubig, but probably not as fast.

    VecX <- matrix(XX,
      nrow = mn,
      ncol = 1
    )

    Ubig <- matrix(
      0,
      mn,
      n * p
    )

    for (i in 1:n) {
      Indx1 <- (i - 1) * m
      Indx2 <- (i - 1) * p
      Ubig[(Indx1 + 1):(Indx1 + m), (Indx2 + 1):(Indx2 + p)] <- U
    }

    UbigCalc <- t(Ubig) %*% Q %*% Ubig

    FCalc <- pracma::pinv(UbigCalc)
    VecMlx <- Ubig %*% (FCalc %*% t(Ubig) %*% (Q %*% VecX))

    Dx <- VecX - VecMlx
    Sobj <- t(Dx) %*% Q %*% Dx
    Mlx <- matrix(VecMlx,
      nrow = m,
      ncol = n
    )

    # Check for convergence or excessive iterations

    if (Count %% 2 == 1) { # check on odd iterations only
      ConvCalc <- abs((Sold - Sobj)) / abs(Sobj) # Convergence Criterion
      if (ConvCalc < ConvLim) {
        ErrFlag <- 0
      }
      if (Count > MaxIter) { # Maximum iterations

        ErrFlag <- 1
        warning("mlpca_e:err3 - Maximum iterations exceeded")
      }
    }

    if (ErrFlag < 0) {
      Sold <- Sobj
      DecomMlx <- svd(Mlx,
        nu = nrow(Mlx),
        nv = ncol(Mlx)
      )
      S <- diag(DecomMlx$d,
        nrow = nrow(Mlx),
        ncol = ncol(Mlx)
      )
      U <- DecomMlx$u
      V <- DecomMlx$v
      XX <- t(XX)
      Q <- K %*% Q %*% t(K)
      K <- t(K)
      m <- nrow(XX)
      n <- ncol(XX)
      U <- matrix(V[, 1:p], ncol = p)
    }
  }

  # All done, clean up and go home

  DecomFinal <- svd(Mlx)
  S <- diag(DecomFinal$d,
    nrow = nrow(Mlx),
    ncol = ncol(Mlx)
  )
  U <- DecomFinal$u
  V <- DecomFinal$v

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

  Ssq <- Sobj

  result <- list(
    "U" = U,
    "S" = S,
    "V" = V,
    "Ssq" = Ssq
  )

  return(result)
}
