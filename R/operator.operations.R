## operator.operations.R
##
##   Copyright (C) 2018, 2019, David Bolin
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Operations with the Pr and Pl operators
#'
#' @description
#' Functions for multiplying and solving with the \eqn{P_r} and \eqn{P_l} operators as well as
#' the latent precision matrix \eqn{Q = P_l C^{-1}P_l} and covariance matrix
#' \eqn{\Sigma = P_r Q^{-1} P_r^T}.
#' These operations are done without first assembling \eqn{P_r}, \eqn{P_l} in order to avoid
#' numerical problems caused by ill-conditioned matrices.
#'
#' @param obj rSPDE object
#' @param v vector to apply the operation to
#' @param transpose set to TRUE if the operation should be performed with the transposed object
#'
#' @return A vector with the values of the operation
#'
#' @details \code{Pl.mult}, \code{Pr.mult}, and \code{Q.mult} multiplies the vector with the respective object.
#' Changing \code{mult} to \code{solve} in the function names multiplies the vector with the inverse of the
#' object. \code{Qsqrt.mult} and \code{Qsqrt.solve} performs the operations with the square-root type object
#' \eqn{Q_r = C^{-1/2}P_l} defined so that \eqn{Q = Q_r^T Q_r}.
#'
#' @name operator.operations
NULL

#' @rdname operator.operations
#' @export
Pr.mult <- function(obj, v, transpose = FALSE) {
  if (class(obj) != "rSPDEobj") {
    stop("obj is not of class rSPDE.obj")
  }
  if (transpose) {
    v <- obj$Pr.factors$Phi %*% v
  }
  I <- Matrix::Diagonal(dim(obj$C)[1])
  if (!is.null(obj$Pr.factors$roots)) {
    for (i in 1:length(obj$Pr.factors$roots)) {
      if (transpose) {
        v <- t(t(v) %*% obj$Pr.factors$factor[[i]])
      } else {
        v <- obj$Pr.factors$factor[[i]] %*% v
      }
    }
  }
  if (!transpose) {
    v <- obj$Pr.factors$Phi %*% v
  }
  return(v)
}


#' @rdname operator.operations
#' @export
Pr.solve <- function(obj, v, transpose = FALSE) {
  if (class(obj) != "rSPDEobj") {
    stop("obj is not of class rSPDE.obj")
  }
  I <- Matrix::Diagonal(dim(obj$C)[1])
  if (!transpose) {
    v <- solve(obj$Pr.factors$Phi, v)
  }
  if (!is.null(obj$Pr.factors$roots)) {
    for (i in 1:length(obj$Pr.factors$roots)) {
      if (transpose) {
        v <- solve(t(obj$Pr.factors$factor[[i]]), v)
      } else {
        v <- solve(obj$Pr.factors$factor[[i]], v)
      }
    }
  }
  if (transpose) {
    v <- solve(obj$Pr.factors$Phi, v)
  }
  return(v)
}

#' @rdname operator.operations
#' @export
Pl.mult <- function(obj, v, transpose = FALSE) {
  if (class(obj) != "rSPDEobj") {
    stop("obj is not of class rSPDE.obj")
  }
  I <- Matrix::Diagonal(dim(obj$C)[1])
  if (transpose) {
    M <- t(obj$CiL)
    v <- obj$C %*% v
    if (obj$Pl.factors$k > 0) {
      for (i in 1:obj$Pl.factors$k) {
        v <- t(t(v) %*% obj$CiL)
      }
    }
  }
  if (!is.null(obj$Pl.factors$roots)) {
    for (i in 1:length(obj$Pl.factors$roots)) {
      if (transpose) {
        v <- t(t(v) %*% obj$Pl.factors$factor[[i]])
      } else {
        v <- obj$Pl.factors$factor[[i]] %*% v
      }
    }
  }
  if (!transpose) {
    if (obj$Pl.factors$k > 0) {
      for (i in 1:obj$Pl.factors$k) {
        v <- obj$CiL %*% v
      }
    }
    v <- obj$C %*% v
  }
  v <- obj$Pl.factors$scaling * v
  return(v)
}

#' @rdname operator.operations
#' @export
Pl.solve <- function(obj, v, transpose = FALSE) {
  if (class(obj) != "rSPDEobj") {
    stop("obj is not of class rSPDE.obj")
  }
  I <- Matrix::Diagonal(dim(obj$C)[1])
  if (!transpose) {
    v <- obj$Ci %*% v
    if (obj$Pl.factors$k > 0) {
      for (i in 1:obj$Pl.factors$k) {
        v <- solve(obj$CiL, v)
      }
    }
  }
  if (!is.null(obj$Pl.factors$roots)) {
    for (i in 1:length(obj$Pl.factors$roots)) {
      if (transpose) {
        v <- solve(t(obj$Pl.factors$factor[[i]]), v)
      } else {
        v <- solve(obj$Pl.factors$factor[[i]], v)
      }
      
    }
  }
  if (transpose) {
    if (obj$Pl.factors$k > 0) {
      for (i in 1:obj$Pl.factors$k) {
        v <- solve(t(obj$CiL), v)
      }
    }
    v <- obj$Ci %*% v
  }
  v <-  v / obj$Pl.factors$scaling
  return(v)
}

#' @rdname operator.operations
#' @export
Q.mult <- function(obj, v) {
  if (class(obj) != "rSPDEobj") {
    stop("obj is not of class rSPDE.obj")
  }
  v <- Pl.mult(obj, v)
  v <- obj$Ci %*% v
  v <- Pl.mult(obj, v, transpose = TRUE)
  return(v)
}

#' @rdname operator.operations
#' @export
Q.solve <- function(obj, v) {
  if (class(obj) != "rSPDEobj") {
    stop("obj is not of class rSPDE.obj")
  }
  v <- Pl.solve(obj, v, transpose = TRUE)
  v <- obj$C %*% v
  v <- Pl.solve(obj, v)
  return(v)
}

#' @rdname operator.operations
#' @export
Qsqrt.mult <- function(obj, v, transpose = FALSE) {
  if (class(obj) != "rSPDEobj") {
    stop("obj is not of class rSPDE.obj")
  }
  if (transpose) {
    v <- sqrt(obj$Ci) %*% v
    v <- Pl.mult(obj, v, transpose = TRUE)
  } else {
    v <- Pl.mult(obj, v)
    v <- sqrt(obj$Ci) %*% v
  }
  return(v)
}

#' @rdname operator.operations
#' @export
Qsqrt.solve <- function(obj, v, transpose = FALSE) {
  if (class(obj) != "rSPDEobj") {
    stop("obj is not of class rSPDE.obj")
  }
  if (transpose) {
    v <- Pl.solve(obj, v, transpose = TRUE)
    v <- sqrt(obj$C) %*% v
  } else {
    v <- sqrt(obj$C) %*% v
    v <- Pl.solve(obj, v)
  }
  return(v)
}

#' @rdname operator.operations
#' @export
Sigma.mult <- function(obj, v) {
  if (class(obj) != "rSPDEobj") {
    stop("obj is not of class rSPDE.obj")
  }
  v <- Pr.mult(obj, v, transpose = TRUE)
  v <- Q.solve(obj, v)
  v <- Pr.mult(obj, v)
  return(v)
}

#' @rdname operator.operations
#' @export
Sigma.solve <- function(obj, v) {
  if (class(obj) != "rSPDEobj") {
    stop("obj is not of class rSPDE.obj")
  }
  v <- Pr.solve(obj, v)
  v <- Q.mult(obj, v)
  v <- Pr.solve(obj, v, transpose = TRUE)
  return(v)
}
