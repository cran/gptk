.complexLog <-
function (x) {
  if ( is.double(x) & x>0 ) {
      y <- log(x)
  } else {
      if ( is.double(x) & x<0 )
          warning("Log of negative real number, using complex log!")
      y <- log(x+0i)
  }
  return ( y )
}
.dist2 <-
function (x, x2) {
  xdim <- dim(as.matrix(x))
  x2dim <- dim(as.matrix(x2))

  xMat <- array(apply(as.matrix(x*x),1,sum), c(xdim[1], x2dim[1]))
  x2Mat <- t(array(apply(as.matrix(x2*x2),1,sum), c(x2dim[1], xdim[1])))

  if ( xdim[2] != x2dim[2] )
    stop("Data dimensions are not matched.")

  n2 <-   xMat+x2Mat-2*tcrossprod(x, x2)

  return (n2)
}
.distfit <-
function(data, dist = "normal") {
  if (dist == "gamma") {
    cdf <- qgamma 
  }

  else if (dist == "normal") {
    cdf <- qnorm
  }

  else {
    stop("Unknown distribution.")
  }

  t <- optim(c(1, 1), fn=.distfit_obj, gr=NULL, data, cdf)

  return (t)
}
.distfit_obj <-
function(theta, y, cdf) {
  p <- c(.05, .25, .50, .75, .95)
  x <- cdf(p, theta[1], theta[2])
  r <- .5 * sum((x - y)^2)

  return (r)
}
.fn_line <-
function (linemin, fun, para0, direction, ...) {
  ## y = fn (x)
  func <- function(x, ...) fun(x, ...)
  
  para <- para0 + linemin * direction
  
  ans <- func(para, ...)
  
  return (ans)
}
.gradLnDiffErfs <-
function(x1, x2, fact1, fact2) {
  m <- pmin(as.matrix(x1)^2, as.matrix(x2)^2)
  dlnPart <- 2/sqrt(pi) * (exp(-x1^2 + m) * fact1 - exp(-x2^2 + m) * fact2)

  g <- list(dlnPart=dlnPart, m=m)
  return (g)

}
.jitChol <-
function ( M, Num=10, silent=FALSE ) {
  jitter <- 0
  jitter1 <- abs(mean(diag(M)))*1e-6
  eyeM <- diag( 1, nrow=length(M[,1]), ncol=length(M[1,]) )

  for ( i in 1:Num ) {
    ## clear the last error message
    try(stop(""),TRUE)

    Ch <- try( chol( M + jitter*eyeM ), silent=TRUE )

    nPos <- grep("not positive definite",  geterrmessage())

    if ( length(nPos) != 0 ) {
      jitter1 <- jitter1*10
      jitter <- jitter1

      if (! silent) {
        warnmsg <- paste("Matrix is not positive definite, adding",
                         signif(jitter,digits=4), "jitter!")
        warning(warnmsg)
      }
    }
    else break
  }

  return (list(chol=Ch, jitter=jitter))
}
.jitCholInv <-
function ( M, Num=10, silent=FALSE ) {
  jitter <- 0
  jitter1 <- abs(mean(diag(M)))*1e-6
  eyeM <- diag( 1, nrow=length(M[,1]), ncol=length(M[1,]) )

  for ( i in 1:Num ) {

    ## clear the last error message
    try(stop(""),TRUE)

    Ch <- try( chol( M + jitter*eyeM ), silent=TRUE )

    nPos <- grep("not positive definite",  geterrmessage())

    if ( length(nPos) != 0 ) {
      jitter1 <- jitter1*10
      jitter <- jitter1

      if (! silent) {
        warnmsg <- paste("Matrix is not positive definite, adding",
                         signif(jitter,digits=4), "jitter!")
        warning(warnmsg)
      }
    }
    else break
  }

  invCh <- try (solve( Ch, eyeM ), silent=TRUE)

  if ( class(invCh) == "try-error" ) {
    return (NaN)
  }
  else {
    invM <- invCh %*% t(invCh)

    if ( jitter == 0 ) {
      ans <- list(invM=invM, jitter=jitter, chol=Ch)
    }
    else ans <- list(invM=invM, jitM=M+jitter*eyeM , jitter=jitter, chol=Ch)

    return (ans)
  }
}
.kernFactors <-
function (kern, factorType) {
  factors <- list()

  if ( length(kern$transforms) > 0 ) {
    funcName <- paste(kern$type, "KernExtractParam", sep="")
    func <- get(funcName, mode="function")
    params <- func(kern)

    for (i in seq(along=kern$transforms)) {
      factors[[i]] <- list()
      factors[[i]]$index <- kern$transforms[[i]]$index
      funcName <- optimiDefaultConstraint(kern$transforms[[i]]$type)
      func <- get(funcName$func, mode="function")
      if (funcName$hasArgs)
        factors[[i]]$val <- func(params[factors[[i]]$index], factorType, kern$transformArgs[[i]])
      else
        factors[[i]]$val <- func(params[factors[[i]]$index], factorType)
    }
  }
  return (factors)
}
.kernTestCombinationFunction <-
function (kern1, kern2) {
  if (kern1$type == "selproj" && kern2$type == "selproj")
    funcName <- paste(kern1$comp[[1]]$type, "X", kern2$comp[[1]]$type, "KernCompute", sep="")
  else
    funcName <- paste(kern1$type, "X", kern2$type, "KernCompute", sep="")

  if ( !exists(funcName, mode="function") ) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}
.multiKernCacheBlock <-
function(kern, fhandle, i, j, x1, x2=NULL, arg1, arg2=NULL) {

  global_cache <- get("cache", envir = kern$cache)
  if ((length(global_cache) >= i) && (length(global_cache[[i]]) >= j))
    cache <- global_cache[[i]][[j]]
  else
    cache <- list()
  key <- c(x1, x2)

  for (k in seq(along=cache)) {
    if (length(key) == length(cache[[k]]$key) && all(key == cache[[k]]$key)) {
      #cat("multiKernCacheBlock: cache hit ", i, j, x1, x2, "\n")
      return (cache[[k]]$value)
    }
  }

  #cat("multiKernCacheBlock: cache miss ", i, j, x1, x2, "\n")
  # No match if we get all the way here
  if (is.null(arg2)) {
    if (is.null(x2))
      K <- fhandle(arg1, x1)
    else
      K <- fhandle(arg1, x1, x2)
  }
  else {
    if (is.null(x2))
      K <- fhandle(arg1, arg2, x1)
    else
      K <- fhandle(arg1, arg2, x1, x2)
  }
  cache <- append(cache, list(list(key=key, value=K)))
  if (length(global_cache) < i)
    global_cache[[i]] <- list()
  global_cache[[i]][[j]] <- cache
  assign("cache", global_cache, envir=kern$cache)
  return(K)
}
.multiKernComputeBlock <-
function (kern, i, j, x1, x2=NULL) {
  if ( i==j ) {
    funcName <- paste(kern$comp[[i]]$type, "KernCompute", sep="")
    transpose <- 0
    arg1 <- kern$comp[[i]]

    func <- get(funcName, mode="function")
    if (kern$fixBlocks[[i]] && kern$fixBlocks[[j]]) {
      K <- .multiKernCacheBlock(kern, func, i, j, arg1=arg1, x1=x1, x2=x2)
    }
    else {
      if (is.null(x2))
        K <- func(arg1, x1)
      else
        K <- func(arg1, x1, x2)
    }
  } else {

    if ( j<i ) {
      funcName <- paste(kern$block[[i]]$cross[j], "KernCompute", sep="")
      transpose <- kern$block[[i]]$transpose[j]
    } else {
      funcName <- paste(kern$block[[j]]$cross[i], "KernCompute", sep="")
      transpose <- !kern$block[[j]]$transpose[i]
    }

    if ( transpose ) {
      arg1 <- kern$comp[[j]]
      arg2 <- kern$comp[[i]]
    } else {
      arg1 <- kern$comp[[i]]
      arg2 <- kern$comp[[j]]      
    }
    
    func <- get(funcName, mode="function")
    if (kern$fixBlocks[[i]] && kern$fixBlocks[[j]]) {
      K <- .multiKernCacheBlock(kern, func, i, j, arg1=arg1, arg2=arg2, x1=x1, x2=x2)
    }
    else {
      if (is.null(x2))
        K <- func(arg1, arg2, x1)
      else
        K <- func(arg1, arg2, x1, x2)
    }
  }
  return (K)
}
.multiKernGradientBlock <-
function (kern, x, x2, covGrad, i, j) {
  if ( nargs()<6 ) {
    j <- i
    i <- covGrad
    covGrad <- x2
    x2 <- array()
  }

  if (kern$fixBlocks[[i]] && kern$fixBlocks[[j]]) {
    if (i==j)
      return (0)
    else
      return (list(g1=0, g2=0))
  }

  if ( i==j ) {
    funcName <- paste(kern$comp[[i]]$type, "KernGradient", sep="")
    transpose <- 0
    arg1 <- kern$comp[[i]]
    factors <- .kernFactors(kern$comp[[i]], "gradfact")

    func <- get(funcName, mode="function")

    if ( is.na(x2) ) {
      g <- func(arg1, x, covGrad)
    } else {
      g <- func(arg1, x, x2, covGrad)
    }
    for (i in seq(along=factors))
      g[factors[[i]]$index] <- g[factors[[i]]$index]*factors[[i]]$val
    
  } else {
    if ( j<i ) {
      funcName <- paste(kern$block[[i]]$cross[j], "KernGradient", sep="")
      transpose <- kern$block[[i]]$transpose[j]
    } else {
      funcName <- paste(kern$block[[j]]$cross[i], "KernGradient", sep="")
      transpose <- kern$block[[j]]$transpose[i]
    }

    if ( transpose ) {
      arg1 <- kern$comp[[j]]
      factors1 <- .kernFactors(kern$comp[[j]], "gradfact")
      arg2 <- kern$comp[[i]]
      factors2 <- .kernFactors(kern$comp[[i]], "gradfact")
    } else {
      arg1 <- kern$comp[[i]]
      factors1 <- .kernFactors(kern$comp[[i]], "gradfact")      
      arg2 <- kern$comp[[j]]
      factors2 <- .kernFactors(kern$comp[[j]], "gradfact")
    }

    func <- get(funcName, mode="function")
    if ( any(is.na(x2)) ) {
      gList <- func(arg1, arg2, x, covGrad)
    } else {
      gList <- func(arg1, arg2, x, x2, covGrad)
    }

    g1 <- gList$g1
    g2 <- gList$g2
    
    for (i in seq(along=factors1))
      g1[factors1[[i]]$index] <- g1[factors1[[i]]$index]*factors1[[i]]$val

    for (i in seq(along=factors2))
      g2[factors2[[i]]$index] <- g2[factors2[[i]]$index]*factors2[[i]]$val

    if ( transpose ) {
      g <- g2
      g2 <- g1
      g1 <- g
    }
    g <- list(g1=g1, g2=g2)   
    
  }
  return (g)
}
.Random.seed <-
c(403L, 352L, 1059855081L, -1953086024L, 1330504294L, -780792063L, 
-385390189L, -1326874046L, -461502988L, -1870389257L, -1615125523L, 
568134788L, -664849750L, 1822132245L, 824914095L, 697933494L, 
1546063072L, 1743581219L, 1703661377L, 600388144L, -376647602L, 
1223790777L, 1696359627L, 354888058L, 1118876796L, -1369849665L, 
-1359813227L, -1600449108L, 1260187522L, 969741661L, -1405473705L, 
-983004354L, 1234359160L, -653049765L, 1855621049L, -1942055416L, 
77669558L, -890584079L, 1414696483L, -951505294L, -1036118844L, 
-129080697L, 1619925597L, -1826579372L, 129726074L, -1307343419L, 
982951711L, 678806182L, -1797528528L, -1706010925L, 1771107569L, 
-421507040L, 261785278L, -542174519L, 1743054971L, -1405417910L, 
2006598316L, 1275221423L, 1081652581L, -1346862308L, 1005037426L, 
68778989L, -1268200409L, -1432980402L, -451782488L, 396900011L, 
-872209463L, -570572840L, -1043835258L, -154661023L, -1293630733L, 
-370131614L, 987454420L, -826663209L, -131947315L, 531833956L, 
914008650L, -316854347L, 1624986383L, 484896790L, 1437563584L, 
-1103291261L, 169584801L, 913637008L, -1856193170L, -1673000167L, 
-1413558357L, 788540826L, 1490737500L, -914853793L, 980697141L, 
-28634164L, -901146846L, 138619261L, 1759003255L, -1326761058L, 
786557528L, 1067893691L, 544605529L, -1357488408L, 1879006998L, 
-1454968687L, -439464637L, -308635566L, 122825444L, 1565093159L, 
-1381091971L, 994090100L, -305232294L, 1165017765L, 935659711L, 
1801501510L, -1799353904L, -250995725L, 1720730385L, -1164455104L, 
48167710L, -1116912663L, -1062680293L, 264561450L, -1617615796L, 
667036687L, -1262319803L, -1272434436L, 570982226L, -1544359731L, 
546046087L, 1290593134L, -235039416L, 18653515L, -2093214551L, 
869197304L, 404639014L, 517451329L, 1969563603L, -1252229758L, 
-1080159820L, 1788198199L, -1330120659L, -196942524L, 721281642L, 
-1193931563L, 1896860271L, 956792694L, -2040381280L, -1123985053L, 
1224776193L, 1679886704L, -1391642354L, 1236863865L, -54263029L, 
587172922L, 2015472444L, -2002845569L, 1018635221L, -1040953876L, 
90005314L, -615096163L, -661858025L, -1246515074L, -1193670600L, 
-1823334629L, -315560839L, -680659768L, -2016078602L, -2032794319L, 
1081292131L, -1177823438L, 36438148L, -134259513L, -253669475L, 
1711242004L, 1769277882L, 1942396933L, -1030302241L, -840914842L, 
1898485744L, -757210861L, -1986052431L, -467469984L, 1494044414L, 
109521545L, 94782907L, 1071300490L, 358655596L, 611784687L, -1524387803L, 
-288145188L, -205183310L, 1403144621L, 482950503L, 1954376334L, 
998339432L, -1016157333L, -150956791L, -22385768L, 107735750L, 
-1559179231L, -733834061L, 1056671266L, 2053519380L, -85345129L, 
-1773594995L, 819754916L, 1249143178L, -1428072971L, -1715681201L, 
-741482666L, 624408064L, -48560317L, -697770271L, 703282256L, 
-230718290L, -156000167L, -737886101L, -1081182758L, 450551452L, 
-1820921697L, 327846645L, -1775298676L, -261605790L, 2083621437L, 
697939127L, -98571426L, 1336069784L, -1686327813L, -1444221031L, 
174082088L, -201609002L, 471029986L, 546538664L, 2048315084L, 
1334827516L, -307862606L, 91814896L, -1414506044L, 915598680L, 
-927970950L, -1947158096L, -2088697412L, -116507052L, 1754917314L, 
-230782176L, 728931260L, -2123259824L, 227664866L, -1491302552L, 
-1052762292L, -260047108L, 905248498L, -125895424L, 1767847028L, 
1862311432L, -1074963014L, 863782352L, 1514037228L, -856720492L, 
-980725550L, 2123910240L, 521672012L, 528395040L, 1273805442L, 
724574696L, 80782316L, 1975551548L, 1372601970L, 581817200L, 
-1791881948L, 178768440L, 156496282L, 2103297488L, 1209693628L, 
536618644L, 369214082L, -1166072384L, -740962500L, 181203600L, 
1823213346L, 32183816L, 1546933004L, 2133494684L, 223618066L, 
57521024L, 591129364L, -238967512L, 153563194L, 1444058064L, 
1065945708L, -2053575628L, -138487662L, 900596064L, -803267252L, 
-1599009312L, -1615531102L, 1959222184L, -1165857652L, 1454130108L, 
-202331022L, 629311472L, 358020164L, 576538712L, -29581254L, 
-134606736L, 1014728252L, 180160276L, 904157890L, -2060909856L, 
364571772L, 1738075344L, 495713186L, -1431245272L, -1204451060L, 
805918652L, -1307453134L, -33865536L, 181891636L, 1115572936L, 
1183230650L, 700298384L, -1972922644L, -302547372L, 1272678034L, 
369805088L, 2113238988L, -1749499680L, 1236926146L, 1894822952L, 
778967980L, 1455766204L, -936637966L, 885839408L, 1525495716L, 
-1484295752L, -994475174L, -1339296624L, -1981535364L, -1838815468L, 
342936898L, -1410175232L, -392703044L, 182439760L, 728762658L, 
1762419272L, 846794764L, 1276957788L, -137907758L, -1864188416L, 
-1677024940L, -77770520L, -924360454L, 335464144L, 1250218284L, 
922247412L, 2060094866L, -111519392L, -193898484L, -1117469984L, 
-893998494L, -1001373016L, 485288524L, -1209700484L, 1025316146L, 
2044978800L, -476941500L, -897816616L, -962525318L, 401782320L, 
666134460L, 1146926676L, 1594380354L, -1189580640L, -844974404L, 
1210673616L, 905462626L, -224746776L, -350502708L, -731422212L, 
-541462030L, -1076098432L, 115514868L, 1937028744L, 32726842L, 
535024464L, 1928064108L, 1439302036L, -634997678L, 1542455776L, 
-272900916L, 1279345824L, -1752894078L, 266389736L, -1966115092L, 
686511420L, 1393626226L, 1532881008L, 1489862692L, 1414253624L, 
-2037993702L, -2064945968L, 603766588L, -1881151596L, -1831447678L, 
70470080L, 1158970940L, 497534352L, -1774633438L, 1000514824L, 
1164914316L, -2054119652L, 113295506L, -1781089920L, 1535403540L, 
-1240228056L, -105309638L, 1531395536L, -1910691476L, -1004928076L, 
670193682L, -178653728L, 1929299916L, 1450565856L, -1676109790L, 
800743720L, 836738572L, 1075947324L, -2089892878L, -2117985552L, 
2060209732L, 87278168L, -820031686L, -272636432L, 1029920444L, 
2027993364L, -170900798L, 709129312L, 1417675644L, 1037834192L, 
-1059364574L, -23856472L, -1911972980L, 968276156L, -1766723790L, 
1800514112L, 1873476788L, 349250120L, 268363194L, 1830060560L, 
-1854606612L, -1683250860L, 1084014866L, 921357728L, 493920972L, 
1814109792L, -949454910L, -1413257432L, -1953009748L, -753283441L, 
1123092953L, -191529362L, 1717345788L, 135243629L, 827756823L, 
1027128544L, 848097982L, 572309643L, -1587077363L, 1865199418L, 
-820811912L, -1650262751L, -1697604141L, -960520972L, -2058381950L, 
-1865969641L, -169249055L, 100134950L, -1716607708L, -1764438571L, 
-1654443649L, -1192644584L, 850415414L, 881875395L, 73979909L, 
-2039620702L, -324368496L, 2102977337L, 2619371L, -242575140L, 
1200687370L, -1750140929L, 580975177L, -2112872482L, 487026572L, 
-1334010019L, 1620166503L, 654507152L, 1919783598L, -918171813L, 
93331261L, -1695441494L, 916845128L, -829173135L, -111270941L, 
2117717732L, 602560018L, -76881433L, 710925745L, -52139018L, 
643753492L, -1391523547L, 126841391L, -718624792L, -1274965818L, 
-1011235469L, 6704597L, -295652942L, -1856073856L, 1998033769L, 
238385115L, 1582776108L, -1738952390L, 1718357999L, 327424633L, 
2123271438L, 324543708L, 2142157645L, 21615287L, -1825529792L, 
270180638L, 1954565099L, -1126090899L, -1985675814L, 1675898776L, 
-521372159L, -1120825421L, -1264894252L, 1375963106L, -107271561L, 
-523414591L, 1711526214L, 1472987268L, 515229237L, -839643617L, 
-1040336520L, 50903126L, -155298589L, 395281061L, 38027074L, 
-456598224L, -730407591L, 2095779531L, -1291360452L, -641929366L, 
-1879089761L, 754283561L, -868391234L, 1761707948L, 437357821L, 
400881735L, -1545222352L, -2075811762L, 68201595L, -1415544867L, 
1626533386L, -664620248L, -166068847L, -346897533L, 611362948L, 
9914034L, -2038408697L, 466659921L, -1837162026L, 1647103156L, 
677427397L, -804402673L, 1504813704L, -195274842L, 160440147L, 
2024169653L, -816234990L, 2053495776L, -199168439L, -2130383749L, 
1163758924L, -1635174182L, 104745167L, -902556519L, -1324427090L, 
-46037572L, -1263657427L, 326990295L, 877826080L, -602180226L, 
546040651L, 2111418701L, -974155910L, 1359798712L, -71728287L, 
-1029893357L, 121236276L, 1442015554L, -1764020905L, -1664661343L, 
577773926L, -397946012L, 124385941L, -938263617L, -129790120L, 
1912522102L, -1673844861L, 1058536901L, -1818721566L, 1980148560L, 
964343417L, 774781355L, 313167644L, 269703242L, 2024754239L, 
621530249L, 34318494L, 2045486412L, -1870691043L, 354893991L, 
-1167708208L, 2099996526L, 704417307L, 41275771L)
