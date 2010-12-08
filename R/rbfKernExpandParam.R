rbfKernExpandParam <-
function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$inverseWidth <- params[1]
  kern$variance <- params[2]

  return (kern)
}

