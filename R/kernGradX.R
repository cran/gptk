kernGradX <-
function (kern, x, x2) {
  funcName <- paste(kern$type, "KernGradX", sep="")
  func <- get(funcName, mode="function")
  k <- func(kern, x, x2)
  return (k)
}

