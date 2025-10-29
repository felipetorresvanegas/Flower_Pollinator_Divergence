constructGradientRRR = function (hM, focalVariable, non.focalVariables = list(), ngrid = 20, 
          coordinates = list()) 
{
  if (!is.null(coordinates)) {
    if (!all(names(coordinates) %in% hM$rLNames)) 
      stop("all names in coordinates should be names of random levels")
  }
  if (is.null(hM$XData)) 
    stop("needs model defined using 'XData' and 'XFormula'")
  if (is.numeric(non.focalVariables) && length(non.focalVariables) == 
      1) {
    defType <- non.focalVariables
    if (!(defType %in% 1:2)) {
      warning("only type 1 and 2 are allowed as a single number: setting to 2")
      defType <- 2
    }
    non.focalVariables <- NULL
  }
  else {
    defType <- 2
  }
  non.focalNames = names(non.focalVariables)
  Mode <- function(x, na.rm = FALSE) {
    if (na.rm) {
      x = x[!is.na(x)]
    }
    ux <- unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
  }
  vars = hM$covNames[-1]
  nvars = length(vars)
  factors = rep(FALSE, nvars)
  focal = NA
  non.focals = NULL
  types = NULL
  vals = list()
  for (i in seq_len(nvars)) {
    if (vars[i] == focalVariable) {
      focal = i
    }
    else {
      non.focals = c(non.focals, i)
      found = FALSE
      for (j in seq_len(length(non.focalVariables))) {
        if (vars[i] == non.focalNames[[j]]) {
          found = TRUE
          type = as.numeric(non.focalVariables[[j]][[1]])
          types = c(types, type)
          if (type == 3) {
            vals[[length(vals) + 1]] = non.focalVariables[[j]][[2]]
          }
          else {
            vals[[length(vals) + 1]] = NA
          }
        }
      }
      if (!found) {
        types = c(types, defType)
        vals[[length(vals) + 1]] = NA
      }
    }
    switch(class(hM$XData)[1L], matrix = {
      if (is.factor(hM$XData[, vars[i]])) {
        factors[i] = TRUE
      }
    }, data.frame = {
      if (is.factor(hM$XData[, vars[i]])) {
        factors[i] = TRUE
      }
    }, list = {
      if (is.factor(hM$XData[[1]][, vars[i]])) {
        factors[i] = TRUE
      }
    })
  }
  switch(class(hM$XData)[1L], matrix = {
    nz = 1
  }, data.frame = {
    nz = 1
  }, list = {
    nz = hM$ns
    XDataNewList = list()
  })
  for (case in 1:nz) {
    switch(class(hM$XData)[1L], matrix = {
      XData = hM$XData
    }, data.frame = {
      XData = hM$XData
    }, list = {
      XData = hM$XData[[case]]
    })
    f.focal = factors[focal]
    v.focal = XData[, vars[focal]]
    if (f.focal) {
      xx = factor(levels(v.focal), levels = levels(v.focal))
      ngrid = length(xx)
    }
    else {
      mi = min(v.focal)
      ma = max(v.focal)
      xx = seq(mi, ma, length.out = ngrid)
    }
    XDataNew = data.frame(xx, stringsAsFactors = TRUE)
    colnames(XDataNew) = vars[focal]
    for (i in seq_len(length(non.focals))) {
      non.focal = non.focals[i]
      type = types[i]
      val = vals[[i]]
      f.non.focal = factors[non.focal]
      v.non.focal = XData[, vars[non.focal]]
      if (f.non.focal) {
        if (type == 1) {
          XDataNew[, vars[non.focal]] = Mode(v.non.focal)
        }
        if (type == 2) {
          mymnm = multinom(v.non.focal ~ v.focal)
          yy = predict(mymnm, newdata = data.frame(v.focal = xx, 
                                                   stringsAsFactors = TRUE))
          XDataNew[, vars[non.focal]] = yy
        }
        if (type == 3) {
          XDataNew[, vars[non.focal]] = rep(val, ngrid)
        }
      }
      if (!f.non.focal) {
        v.non.focal = XData[, vars[non.focal]]
        if (type == 1) {
          XDataNew[, vars[non.focal]] = mean(v.non.focal)
        }
        if (type == 2) {
          mylm = lm(v.non.focal ~ v.focal)
          yy = predict(mylm, newdata = data.frame(v.focal = xx, 
                                                  stringsAsFactors = TRUE))
          XDataNew[, vars[non.focal]] = yy
        }
        if (type == 3) {
          XDataNew[, vars[non.focal]] = rep(val, ngrid)
        }
      }
    }
    if (inherits(hM$XData, "list")) {
      XDataNewList[[case]] = XDataNew
    }
  }
  if (inherits(hM$XData, "list")) {
    XDataNew = XDataNewList
  }
  dfPiNew = matrix("new_unit", ngrid, hM$nr)
  colnames(dfPiNew) = hM$rLNames
  dfPiNew = as.data.frame(dfPiNew, stringsAsFactors = TRUE)
  rLNew = vector("list", hM$nr)
  names(rLNew) = hM$rLNames
  for (r in seq_len(hM$nr)) {
    rLname <- hM$rLNames[r]
    coord <- coordinates[[rLname]]
    if (!is.null(coord) && !(is.numeric(coord) || coord %in% 
                             c("c", "i"))) 
      stop("'coordinates' must be 'c', 'i' or numeric in random level ", 
           sQuote(rLname))
    rL1 = hM$rL[[r]]
    xydata = rL1$s
    if (!is.null(xydata)) {
      if (is(xydata, "Spatial")) {
        if (!is.null(coord) && coord == "i") 
          centre <- rep(Inf, NCOL(coordinates(xydata)))
        else if (is.numeric(coord)) 
          centre <- coord
        else centre <- as.data.frame(t(colMeans(coordinates(xydata))))
        rownames(centre) <- "new_unit"
        coordinates(centre) <- colnames(centre)
        proj4string(centre) <- proj4string(xydata)
        xydata <- rbind(xydata, centre)
      }
      else {
        if (!is.null(coord) && coord == "i") 
          centre <- rep(Inf, NCOL(xydata))
        else if (is.numeric(coord)) 
          centre <- coord
        else centre <- colMeans(xydata)
        xydata = rbind(xydata, new_unit = centre)
      }
    }
    rL1$s = xydata
    distMat = rL1$distMat
    if (!is.null(distMat)) {
      units1 = c(rownames(distMat), "new_unit")
      if (is.numeric(coord)) {
        stop("numeric coordinates are not enabled for 'distMat'")
      }
      else if (!is.null(coord) && coord == "i") {
        newdist <- rep(Inf, NCOL(distMat))
      }
      else {
        newdist <- distMat^2/2
        newdist <- sweep(newdist, 2L, colMeans(newdist), 
                         check.margin = FALSE)
        newdist <- sweep(newdist, 1L, rowMeans(newdist), 
                         check.margin = FALSE)
        newdist <- sqrt(diag(-newdist))
      }
      distMat1 = cbind(distMat, newdist)
      distMat1 = rbind(distMat1, c(newdist, 0))
      rownames(distMat1) = units1
      colnames(distMat1) = units1
      rL1$distMat = distMat1
    }
    rL1$pi = c(rL1$pi, "new_unit")
    rL1$N = rL1$N + 1
    rLNew[[r]] = rL1
  }
  Gradient = list(XDataNew = XDataNew, studyDesignNew = dfPiNew, 
                  rLNew = rLNew)
  return(Gradient)
}
