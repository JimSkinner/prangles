#' Compute Principal Angles between subspaces
#'
#' This function takes two matrices, A and B, each defining a subspace via a
#' set of basis vectors in columns. The Principal Angles between these
#' subspaces are then computed and returned in radians. These matrices do not
#' need to have the same number of rows, but must have the same number of
#' columns; i.e., the subspaces defined must live in the same vector space, but
#' do not need to share dimensionality. Vectors may be provided, and these are
#' converted into single column matrices.
#'
#' @param A (d x n1) matrix containin basis for subspace A.
#' @param B (d x n2) matrix containin basis for subspace B.
#'
#' @return A vector of min(n1, n2) principal angles in descending order (units
#'   in radians).
#'
#' @export
#' @examples
#' A <- matrix(c(1, 0, 0,
#'               0, 1, 0), nrow=3, ncol=2)
#' B <- matrix(c(1, 0, 0), nrow=3, ncol=1)
#' C <- matrix(c(0, 0, 1), nrow=3, ncol=1)
#' D <- matrix(c(B, C), nrow=3, ncol=2)
#'
#' # Subspace B intersects with A, so gives the single PA of 0
#' anglesAB <- prangles(A, B)
#' stopifnot(anglesAB == 0)
#'
#' # Subspace C is orthogonal to A, so has a PA of pi/2
#' anglesAC <- prangles(A, C)
#' stopifnot(all.equal(anglesAC, pi/2))
#'
#' # A and D have PAs {pi/2, 0}, so they intersect along 1 dimension
#' anglesAD <- prangles(A, D)
#' stopifnot(all.equal(anglesAD, c(pi/2, 0)))
prangles <- function(A, B) {
  if (is.vector(A)) A = matrix(A, ncol=1)
  if (is.vector(B)) B = matrix(B, ncol=1)

  Qa = qr.Q(qr(A))
  Qb = qr.Q(qr(B))

  C = svd(crossprod(Qa, Qb))$d
  C = vapply(C, function(x) min(x, 1), numeric(1))
  angles = sort(acos(C), decreasing=TRUE)
  return(angles)
}
