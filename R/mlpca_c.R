options(encoding = "UTF-8")

#' @title Maximum Likelihood Principal Component Analysis for Mode C Error
#' Conditions
#'
#' @description  R implementation of Maximum Likelihood Principal Component
#' Analisys proposed in Wentzell, Peter D., et al. "Maximum likelihood principal
#' component analysis." Journal of Chemometrics: A Journal of the Chemometrics
#' Society 11.4 (1997): 339-366.
#' @author Renan Santos Barbosa
#'
#' @param X MxN matrix of measurements
#' @param Xsd MxN matrix of measurements error standard deviations
#' @param p Rank of the model's subspace, p must be than the minimum of M and N
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
#' @export
#' @encoding UTF-8
#'
#' @example set.seed(123)
#'
#' X <- matrix(rnorm(300,
#'                   mean = 10,
#'                   sd = 5),
#'             nrow = 3,
#'             ncol = 100,
#'             byrow = FALSE)
#'
#' Xsd <- matrix(rgamma(300,
#'                      shape = 1),
#'               nrow = 3,
#'               ncol = 100,
#'               byrow = FALSE)
#'
#' p <- 2
#'
#' results <- RMLPCA::mlpca_c(X,Xsd,p)

mlpca_c <- function(X,Xsd,p){

  m <- dim(x = X)[1]
  n <- dim(x = X)[2]

  if (p > min(m,n)) {
    stop("mlpca_c:err1 - Invalid rank for MLPCA decomposition")
  }

  ml <- dim(x = Xsd)[1]
  nl <- dim(x = Xsd)[2]

  if (m != ml | n != nl){
    stop("mlpca_c:err2 - Dimensions of data and standard deviations do not matchn")
  }


  if(isFALSE(all(Xsd > 0 ))){

    stop("mlpca_c:err3 - Standard deviations must be positive")

  }

  if(isTRUE(any(Xsd == 0))){

    stop("mlpca_c:err4 - Zero value(s) for standard deviations")

  }

  # Initialization -------------------------------------------------------------


  ConvLim <- 1e-10 # Convergence Limit
  MaxIter <- 2000 # Maximum no. of iterations
  VarMUlt <- 1000 # Multiplier for missing data
  VarX <- Xsd^2 # Convert sd's to variances
  IndX <- which(is.na(VarX)) # Find missing values
  VarMax <- max(VarX,na.rm = TRUE) # Maximum variance
  VarX[IndX] <- VarMax*VarMUlt # Give missing values large variance

  # Generate Initial estimates assuming homocedastic errors --------------------

  DecomX <- RSpectra::svds(X,p) # Decompose adjusted matrix
  U <- DecomX$u
  S <- diag(DecomX$d)
  V <- DecomX$v

  Count <- 0 # Loop counter
  Sold <- 0 # Holds last value of objective function
  ErrFlag <- -1 # Loop flag

  while(ErrFlag < 0){

    Count <- Count + 1 # Loop counter

    # Evaluate objective function ----------------------------------------------

    Sobj <- 0 # Initialize sum
    MLX <- matrix(data = 0,
                  nrow = dim(X)[1],
                  ncol = dim(X)[2])

    for(i in 1:n){
      Q = diag(1/VarX[,i]) # Inverse of error covariance matrix
      FInter <- solve(t(U)%*%Q%*%U) # Intermediate calculation
      MLX[,i] <- U %*% (FInter %*% (t(U) %*% (Q %*% X[,i]))) # Max.Lik Estimates
      Dx <- matrix(data = X[,i] - MLX[,i]) # Residual Vector
      Sobj <- Sobj + t(Dx) %*% Q %*% Dx # update objective function
    }


    # Check for convergence or excessive iterations ----------------------------

    if(Count%%2 == 1){ # check on odd iterations only
      ConvCalc <- abs(Sold - Sobj)/Sobj # Convergence Criterion
      if( ConvCalc < ConvLim){

        ErrFlag <- 0

      }
      if( Count > MaxIter){ # Maximum iterations

        ErrFlag <- 1
        warning("mlpca_c:err5 - Maximum iterations exceeded")

      }

    }

    # Now flip matrices for alternating regression -----------------------------

    if( ErrFlag < 0 ){ # Only do this part if not done

      Sold <- Sobj # Save most recent objective function
      DecomMLX <- RSpectra::svds(MLX,p) # Decompose Model values
      U <- DecomMLX$u
      S <- diag(DecomMLX$d)
      V <- DecomMLX$v

      X <- t(X) # Flip matrix
      VarX <- t(VarX) # And the variances
      n <- ncol(X) # Adjust no. of columns
      U <- V # V becomes U in for transpose

    }

    # All done -----------------------------------------------------------------

  }

  DecomFinal <- RSpectra::svds(MLX,p)
  U <- DecomFinal$u
  S <- diag(DecomFinal$d)
  V <- DecomFinal$v
  Ssq <- Sobj

  result <- list("U" = U,
                 "S" = S,
                 "V" = V,
                 "Ssq" = Sobj,
                 "ErrFlag" = ErrFlag)

  return(result)

}
