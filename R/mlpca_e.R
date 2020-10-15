#' Maximum likelihood orincipal component analysis for mode E error conditions
#'
#' @description  Performs maximum likelihood principal components analysis for
#'   mode E error conditions (correlated errors, with a different covariance
#'   matrix for each row, but no error correlation between the rows). Employs an
#'   ALS algorithm.
#'
#' @author Renan Santos Barbosa
#'
#' @param X IxJ matrix of measurements
#' @param Cov JXJXI matrices of measurement error covariance
#' @param p Rank of the model's subspace, p must be than the minimum of I and J
#' @param MaxIter Maximum no. of iterations
#'
#' @return  The parameters returned are the results of SVD on the estimated
#'   subspace. The quantity Ssq represents the sum of squares of weighted
#'   residuals. ErrFlag indicates the convergence condition,
#'   with 0 indicating normal termination and 1 indicating the maximum number of
#'   iterations have been exceeded.
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
#' @examples
#'
#' library(RMLPCA)
#' data(data_clean_e)
#' data(data_error_e)
#' # covariance matrix
#' data(cov_e)
#' data(data_cleaned_mlpca_e)
#' # data that you will usually have on hands
#' data_noisy <- data_clean_e + data_error_e
#'
#' # run mlpca_e with rank p = 1
#' results <- RMLPCA::mlpca_e(
#'   X = data_noisy,
#'   Cov = cov_e,
#'   p = 1
#' )
#'
#' # estimated clean dataset
#' data_cleaned_mlpca <- results$U %*% results$S %*% t(results$V)
mlpca_e <- function(X, Cov, p, MaxIter = 20000) {
  m <- base::dim(x = X)[1]

  n <- base::dim(x = X)[2]

  if (p > base::min(m, n)) {
    stop("mlpca_e:err1 - Invalid rank for MLPCA decomposition")
  }

  ml <- base::dim(x = Cov)[1]
  nl <- base::dim(x = Cov)[2]
  q <- base::dim(x = Cov)[3]


  if (n != ml | n != nl | m != q) {
    stop("mlpca_e:err2 - Invalid dimensions of covariance matrix")
  }

  # Initialization -------------------------------------------------------------


  ConvLim <- 1e-10 # Convergence Limit

  # Calculate the inverse of the full covariance matrix blockwise --------------

  mn <- m * n

  Q <- base::matrix(
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

  ix <- base::matrix(1:mn)
  iy <- base::array(
    base::t(base::array(ix,
      dim = base::c(m, n)
    )),
    base::c(mn, 1)
  )

  K <- Matrix::sparseMatrix(
    i = ix,
    j = iy,
    x = 1,
    dims = c(
      mn,
      mn
    )
  )
  K <- base::as.matrix(K)

  Q <- (base::t(K) %*% Q %*% K)

  # Generate initial estimates assuming an average covariance matrix

  Covavg <- base::rowMeans(Cov, dims = 2)

  res_mlpca_d <- RMLPCA::mlpca_d(
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

    VecX <- base::matrix(XX,
      nrow = mn,
      ncol = 1
    )

    Ubig <- base::matrix(
      0,
      mn,
      n * p
    )

    for (i in 1:n) {
      Indx1 <- (i - 1) * m
      Indx2 <- (i - 1) * p
      Ubig[(Indx1 + 1):(Indx1 + m), (Indx2 + 1):(Indx2 + p)] <- U
    }

    UbigCalc <- base::t(Ubig) %*% Q %*% Ubig

    FCalc <- pracma::pinv(UbigCalc)
    VecMlx <- Ubig %*% (FCalc %*% base::t(Ubig) %*% (Q %*% VecX))

    Dx <- VecX - VecMlx
    Sobj <- base::t(Dx) %*% Q %*% Dx
    Mlx <- base::matrix(VecMlx,
      nrow = m,
      ncol = n
    )

    # Check for convergence or excessive iterations

    if (Count %% 2 == 1) { # check on odd iterations only
      ConvCalc <- base::abs((Sold - Sobj)) / base::abs(Sobj) # Convergence Criterion
      if (ConvCalc < ConvLim) {
        ErrFlag <- 0
      }
      if (Count > MaxIter) { # Maximum iterations

        ErrFlag <- 1
        stop("mlpca_e:err3 - Maximum iterations exceeded")
      }
    }

    if (ErrFlag < 0) {
      Sold <- Sobj
      DecomMlx <- base::svd(Mlx,
        nu = base::nrow(Mlx),
        nv = base::ncol(Mlx)
      )
      S <- base::diag(DecomMlx$d,
        nrow = base::nrow(Mlx),
        ncol = base::ncol(Mlx)
      )
      U <- DecomMlx$u
      V <- DecomMlx$v
      XX <- base::t(XX)
      Q <- K %*% Q %*% base::t(K)
      K <- base::t(K)
      m <- base::nrow(XX)
      n <- base::ncol(XX)
      U <- base::matrix(V[, 1:p], ncol = p)
    }
  }

  # All done, clean up and go home

  DecomFinal <- base::svd(Mlx)
  S <- base::diag(DecomFinal$d,
    nrow = base::nrow(Mlx),
    ncol = base::ncol(Mlx)
  )
  U <- DecomFinal$u
  V <- DecomFinal$v

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

  Ssq <- Sobj

  result <- base::list(
    "U" = U,
    "S" = S,
    "V" = V,
    "Ssq" = Ssq
  )

  return(result)
}
