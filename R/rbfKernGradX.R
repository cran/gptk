rbfKernGradX <-
function (kern, X, X2) {

  gX = array(0, c(dim(as.array(X2))[1], dim(as.array(X2))[2], dim(as.array(X))[1]))
  for (i in 1:dim(X)[1]) {
    gX[, , i] = rbfKernGradXpoint(kern, X[i, ], X2)
  }

  return (gX)
}

