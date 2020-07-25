options(encoding = "UTF-8")

#' @title Maximum Likelihood Principal Component Analysis for Mode B Error
#' Conditions
#'
#' @description Performs maximum likelihood principal components analysis for
#' mode B error conditions (independent errors, homoscedastic within a column).
#' Equivalent to perfoming PCA on data scaled by the error SD, but results are
#' rescaled to the original space.
#'
#' @author Renan Santos Barbosa
#'
#' @param X MxN matrix of measurements
#' @param Xsd MxN matrix of measurements error standard deviations
#' @param p Rank of the model's subspace, p must be than the minimum of M and N
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
#' @import RSpectra
#' @export
#' @encoding UTF-8
#'
#' @example set.seed(123)
#'
#' X = matrix(c(1,2,3,4,4,4,5,4,23,43,55,4),
#' nrow = 3,
#' ncol = 4,
#' byrow = TRUE)
#' Xsd = matrix(c(0.1,0.1,0.1,0.1,
#'               0.1,0.1,0.1,0.1,
#'               0.1,0.1,0.1,0.1),
#'             nrow = 3,
#'             ncol = 4)
#'
#' p = 1
#'
#' mlpca_b(X,Xsd,p)

mlpca_b <- function(X,Xsd,p){

  m <- dim(x = X)[1]
  n <- dim(x = X)[2]

  if (p > min(m,n)) {
    stop("mlpca_b:err1 - Invalid rank for MLPCA decomposition")
  }

  ml <- dim(x = Xsd)[1]
  nl <- dim(x = Xsd)[2]

  if (m != ml | n != nl){
    stop("mlpca_b:err2 - Dimensions of data and standard deviations do not matchn")
  }


  if(isFALSE(all(Xsd > 0 ))){

    stop("mlpca_b:err3 - Standard deviations must be positive")

  }

  if(isTRUE(any(Xsd == 0))){

    stop("mlpca_b:err4 - Zero value(s) for standard deviations")

  }

  SclMat <- Xsd

  Xsc <- X/SclMat

  DecomXsc <- RSpectra::svds(Xsc,p) # Decompose adjusted matrix
  U <- DecomXsc$u
  S <- diag(DecomXsc$d,
            nrow = length(DecomXsc$d),
            ncol = length(DecomXsc$d))
  V <- DecomXsc$v

  XCalc <- (U%*%S%*%t(V)) * SclMat

  DecomXCalc <- RSpectra::svds(XCalc,p) # Decompose adjusted matrix
  U <- DecomXCalc$u
  S <- diag(DecomXCalc$d,
            nrow = length(DecomXCalc$d),
            ncol = length(DecomXCalc$d))
  V <- DecomXCalc$v

  Ssq <- sum(sum(((X - XCalc)/SclMat)^2))

  result <- list("U" = U,
                 "S" = S,
                 "V" = V,
                 "Ssq" = Ssq
                 )

  return(result)

}
