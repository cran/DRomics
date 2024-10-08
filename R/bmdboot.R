## Perform non parametric bootstrap on a selection of items in order to give
## a confidence interval for the BMD values
bmdboot <- function(r, items = r$res$id, niter = 1000,
                    conf.level = 0.95,
                    tol = 0.5, progressbar = TRUE,
                    parallel = c("no", "snow", "multicore"), ncpus) {
  
  # Checks
  if (!inherits(r, "bmdcalc")) {
    stop("Use only with 'bmdcalc' objects, created with the function bmdcalc.")
  }
  
  # bootmethod <- match.arg(bootmethod, c("nonparam", "param"))
  bootmethod <- "nonparam"
  
  if (niter < 1000) {
    warning(strwrap(prefix = "\n", initial = "\n",
                    "A small number of iterations (less than 1000) may not be sufficient
      to ensure a good quality of bootstrap confidence intervals."))
  }
  
  parallel <- match.arg(parallel, c("no", "snow", "multicore"))
  if (parallel == "multicore" && .Platform$OS.type == "windows") {
    parallel <- "snow"
    warning(strwrap(prefix = "\n", initial = "\n",
                    "As the multicore option is not supported on Windows it was replaced by snow."))
  }
  
  if ((parallel == "snow" || parallel == "multicore") && missing(ncpus)) {
    stop("You have to specify the number of available processors to parallelize the bootstrap.")
  }
  
  if (parallel != "no") {
    progressbar <- FALSE
  }
  
  if (progressbar) {
    cat(strwrap("The bootstrap may be long if the number of items and the number of bootstrap iterations is high."), fill = TRUE)
  }
  
  i.items <- match(items, r$res$id)
  nitems <- length(items)
  
  if (!is.numeric(tol) || (tol > 1) || (tol < 0)) {
    stop("Wrong argument 'tol'. If not omitted it must be a number between 0 and 1 (the proportion
    of failure of model fit among bootstrap samples).")
  }
  
  if (!is.numeric(conf.level) || (conf.level >= 1) || (conf.level <= 0)) {
    stop("Wrong argument 'conf.level'. If not omitted it must be a number between 0 and 1 
    (the confidence level of the bootstrap confidence intervals).")
  }
  
  prob.lower <- (1 - conf.level) / 2
  prob.upper <- conf.level + prob.lower
  
  z <- r$z
  x <- r$x
  xdiv100 <- r$x / 100
  
  dose <- r$omicdata$dose
  dosemin <- min(dose[dose != 0])
  dosemax <- max(dose)
  minBMD <- r$minBMD
  ratio2switchinlog <- r$ratio2switchinlog
  
  # progress bar
  if (progressbar) {
    pb <- utils::txtProgressBar(min = 0, max = nitems, style = 3)
  }
  
  ##### Bootstrap for one item ####################
  bootoneitem <- function(i) {
    resitem <- r$res[i.items[i], ]
    # parameter estimates as a named list for starting values
    estimpar <- unlist(resitem[c("b", "c", "d", "e", "f")])
    lestimpar <- as.list(estimpar[!is.na(estimpar)])
    
    # residual standard error
    SDresi <- unlist(resitem["SDres"])
    
    modeli <- as.character(unlist(resitem["model"]))
    nbpari <- unlist(resitem["nbpar"])
    
    # dataset
    datai <- r$omicdata$data[resitem$irow, ]
    dset <- data.frame(signal = datai, dose = dose)
    # removing lines with NA values for the signal
    dset <- dset[stats::complete.cases(dset$signal), ]
    ndata <- nrow(dset)
    
    ############## Model expo ###########
    if (modeli == "exponential") {
      b1 <- lestimpar$b
      d1 <- lestimpar$d
      e1 <- lestimpar$e
      fitted1 <- fExpo(x = dset$dose, d = d1, b = b1, e = e1)
      resid1 <- dset$signal - fitted1
      
      dsetboot <- dset
      fboot <- function(i) {
        
        if (bootmethod == "param") {
          dsetboot[, 1] <- fitted1 + stats::rnorm(ndata, mean = 0, sd = SDresi)
        } else {
          dsetboot[, 1] <- fitted1 + sample(scale(resid1, scale = FALSE), replace = TRUE)
        }
        
        # fit
        if (e1 < 0) {
          nlsboot <- suppressWarnings(try(stats::nls(formula = formExp3p, data = dsetboot, start = lestimpar,
                                                     lower = c(-Inf, -Inf, -Inf),
                                                     upper = c(Inf, Inf, 0), algorithm = "port"),
                                          silent = TRUE))
        } else {
          nlsboot <- suppressWarnings(try(stats::nls(formula = formExp3p, data = dsetboot, start = lestimpar,
                                                     lower = c(-Inf, -Inf, 0), algorithm = "port"),
                                          silent = TRUE))
        }
        
        if (inherits(nlsboot, "nls")) {
          SDresboot <- sqrt(sum(stats::residuals(nlsboot)^2) / (ndata - nbpari))
          bboot <- stats::coef(nlsboot)["b"]
          dboot <- stats::coef(nlsboot)["d"]
          eboot <- stats::coef(nlsboot)["e"]
          y0boot <- dboot
          # ydosemaxboot <- fExpo(x = dosemax, b = bboot, d = dboot, e = eboot) # local not used variable (linter)
          ypboot <- y0boot * (1 + xdiv100 * sign(eboot * bboot))
          BMDpboot <- pmax(invExpo(ypboot, b = bboot, d = dboot, e = eboot), minBMD)
          ysdboot <- y0boot + z * SDresboot * sign(eboot * bboot)
          BMDsdboot <- pmax(invExpo(ysdboot, b = bboot, d = dboot, e = eboot), minBMD)
          return(list(BMDp = BMDpboot, BMDsd = BMDsdboot))
        }
      } # end fboot
      ############ END model expo ###################
      
      ############## Model Hill ###########
    } else if (modeli == "Hill") {
      b1 <- lestimpar$b
      c1 <- lestimpar$c
      d1 <- lestimpar$d
      e1 <- lestimpar$e
      fitted1 <- fHill(x = dset$dose, b = b1, c = c1, d = d1,  e = e1)
      resid1 <- dset$signal - fitted1
      
      dsetboot <- dset
      fboot <- function(i) {
        if (bootmethod == "param") {
          dsetboot[, 1] <- fitted1 + stats::rnorm(ndata, mean = 0, sd = SDresi)
        } else {
          dsetboot[, 1] <- fitted1 + sample(scale(resid1, scale = FALSE), replace = TRUE)
        }
        
        # fit
        nlsboot <- suppressWarnings(try(stats::nls(formula = formHill, data = dsetboot, start = lestimpar,
                                                   lower = c(0, -Inf, -Inf, 0), algorithm = "port"),
                                        silent = TRUE))
        if (inherits(nlsboot, "nls")) {
          SDresboot <- sqrt(sum(stats::residuals(nlsboot)^2) / (ndata - nbpari))
          bboot <- stats::coef(nlsboot)["b"]
          cboot <- stats::coef(nlsboot)["c"]
          dboot <- stats::coef(nlsboot)["d"]
          eboot <- stats::coef(nlsboot)["e"]
          y0boot <- dboot
          # ydosemaxboot <- fHill(x = dosemax, b = bboot, c = cboot, d = dboot, e = eboot) # local not used variable (linter)
          ypboot <- y0boot * (1 + xdiv100 * sign(cboot * dboot))
          BMDpboot <- pmax(invHill(ypboot, b = bboot, c = cboot, d = dboot, e = eboot), minBMD)
          ysdboot <- y0boot + z * SDresboot * sign(cboot * dboot)
          BMDsdboot <- pmax(invHill(ysdboot, b = bboot, c = cboot, d = dboot, e = eboot), minBMD)
          
          return(list(BMDp = BMDpboot, BMDsd = BMDsdboot))
        }
      } # end fboot
      ############ END model Hill ###################
      
      ############## Linear model ###########
    } else if (modeli == "linear") {
      b1 <- lestimpar$b
      d1 <- lestimpar$d
      fitted1 <- flin(x = dset$dose, b = b1, d = d1)
      resid1 <- dset$signal - fitted1
      
      dsetboot <- dset
      fboot <- function(i) {
        if (bootmethod == "param") {
          dsetboot[, 1] <- fitted1 + stats::rnorm(ndata, mean = 0, sd = SDresi)
        } else {
          dsetboot[, 1] <- fitted1 + sample(scale(resid1, scale = FALSE), replace = TRUE)
        }
        
        # fit
        linboot <- stats::lm(signal ~ dose, data = dsetboot)
        SDresboot <- sqrt(sum(stats::residuals(linboot)^2) / (ndata - nbpari))
        bboot <- stats::coef(linboot)[2]
        dboot <- stats::coef(linboot)[1]
        y0boot <- dboot
        # ydosemaxboot <- flin(x = dosemax, b = bboot, d = dboot) # local not used variable (linter)
        ypboot <- y0boot * (1 + xdiv100 * sign(bboot))
        BMDpboot <- pmax(invlin(ypboot, b = bboot, d = dboot), minBMD)
        ysdboot <- y0boot + z * SDresboot * sign(bboot)
        BMDsdboot <- pmax(invlin(ysdboot, b = bboot, d = dboot), minBMD)
        
        return(list(BMDp = BMDpboot, BMDsd = BMDsdboot))
      } # end fboot
      
      ############ END linear model ###################
      
      ##### Model Gauss-probit #######
    } else if (modeli == "Gauss-probit") {
      b1 <- lestimpar$b
      c1 <- lestimpar$c
      d1 <- lestimpar$d
      e1 <- lestimpar$e
      f1 <- lestimpar$f
      fitted1 <- fGauss5p(x = dset$dose, c = c1, d = d1, b = b1, e = e1, f = f1)
      resid1 <- dset$signal - fitted1
      
      dsetboot <- dset
      fboot <- function(i) {
        if (bootmethod == "param") {
          dsetboot[, 1] <- fitted1 + stats::rnorm(ndata, mean = 0, sd = SDresi)
        } else {
          dsetboot[, 1] <- fitted1 + sample(scale(resid1, scale = FALSE), replace = TRUE)
        }
        
        # fit
        if (nbpari == 5) {
          nlsboot <- suppressWarnings(try(stats::nls(formula = formGauss5p, data = dsetboot, start = lestimpar,
                                                     lower = c(0, -Inf, -Inf, 0, -Inf), algorithm = "port"),
                                          silent = TRUE))
        } else {
          if (f1 == 0) {
            lestimpar.f0 <- list(b = lestimpar$b, c = lestimpar$c, d = lestimpar$d, e = lestimpar$e)
            nlsboot <- suppressWarnings(try(stats::nls(formula = formprobit, data = dsetboot, start = lestimpar.f0,
                                                       lower = c(0, -Inf, -Inf, 0), algorithm = "port"),
                                            silent = TRUE))
          } else {
            lestimpar.4p <- list(b = lestimpar$b, d = lestimpar$d, e = lestimpar$e, f = lestimpar$f)
            nlsboot <- suppressWarnings(try(stats::nls(formula = formGauss4p, data = dsetboot, start = lestimpar.4p,
                                                       lower = c(0, -Inf, 0, -Inf), algorithm = "port"),
                                            silent = TRUE))
          }
        }
        
        if (inherits(nlsboot, "nls")) {
          SDresboot <- sqrt(sum(stats::residuals(nlsboot)^2) / (ndata - nbpari))
          bboot <- stats::coef(nlsboot)["b"]
          if (nbpari == 5) {
            cboot <- stats::coef(nlsboot)["c"]
          } else {
            cboot <- stats::coef(nlsboot)["d"]
          }
          dboot <- stats::coef(nlsboot)["d"]
          eboot <- stats::coef(nlsboot)["e"]
          fboot <- stats::coef(nlsboot)["f"]
          y0boot <- fGauss5p(x = 0, b = bboot, c = cboot, d = dboot, e = eboot, f = fboot)
          ydosemaxboot <- fGauss5p(x = dosemax, b = bboot, c = cboot, d = dboot, e = eboot, f = fboot)
          
          if (f1 != 0) {
            xextrboot <- eboot + (cboot - dboot) * bboot / (fboot * sqrt(2 * pi))
            yextrboot <- fGauss5p(x = xextrboot, b = bboot, c = cboot, d = dboot, e = eboot, f = fboot)
            
            deltapboot <- abs(y0boot) * xdiv100
            deltasdboot <- z * SDresboot
            
            resBMDp <- suppressWarnings(try(
              calcBMD(y0 = y0boot, delta = deltapboot, xext = xextrboot, yext = yextrboot,
                      dosemin = dosemin, dosemax = dosemax, ydosemax = ydosemaxboot,
                      func = fGauss5pBMR, func_xinlog = fGauss5pBMR_xinlog,
                      b = bboot, c = cboot, d = dboot, e = eboot, g = fboot,
                      minBMD = minBMD, ratio2switchinlog = ratio2switchinlog),
              silent = TRUE))
            
            resBMDsd <- suppressWarnings(try(
              calcBMD(y0 = y0boot, delta = deltasdboot, xext = xextrboot, yext = yextrboot,
                      dosemin = dosemin, dosemax = dosemax, ydosemax = ydosemaxboot,
                      func = fGauss5pBMR, func_xinlog = fGauss5pBMR_xinlog,
                      b = bboot, c = cboot, d = dboot, e = eboot, g = fboot,
                      minBMD = minBMD, ratio2switchinlog = ratio2switchinlog),
              silent = TRUE))
            # return a value only if no problem with nls AND uniroot
            if (!inherits(resBMDsd, "try-error") &&  !inherits(resBMDp, "try-error")) {
              return(list(BMDp = resBMDp$BMD, BMDsd = resBMDsd$BMD))
            }
          } else { # if f1 == 0
            
            ypboot <- y0boot * (1 + xdiv100 * sign(cboot * dboot))
            BMDpboot <- pmax(invprobit(ypboot, b = bboot, c = cboot, d = dboot, e = eboot), minBMD)
            ysdboot <- y0boot + z * SDresboot * sign(cboot * dboot)
            BMDsdboot <- pmax(invprobit(ysdboot, b = bboot, c = cboot, d = dboot, e = eboot), minBMD)
            return(list(BMDp = BMDpboot, BMDsd = BMDsdboot))
          }
        } # end if (inherits(nlsboot, "nls"))
      } # end fboot
      ############ END model Gauss probit ###################
      
      ############ Model log Gauss-probit #################
    } else if (modeli == "log-Gauss-probit") {
      b1 <- lestimpar$b
      c1 <- lestimpar$c
      d1 <- lestimpar$d
      e1 <- lestimpar$e
      f1 <- lestimpar$f
      fitted1 <- fLGauss5p(x = dset$dose, c = c1, d = d1, b = b1, e = e1, f = f1)
      resid1 <- dset$signal - fitted1
      
      dsetboot <- dset
      fboot <- function(i) {
        if (bootmethod == "param") {
          dsetboot[, 1] <- fitted1 + stats::rnorm(ndata, mean = 0, sd = SDresi)
        } else {
          dsetboot[, 1] <- fitted1 + sample(scale(resid1, scale = FALSE), replace = TRUE)
        }
        
        # fit
        if (nbpari == 5) {
          nlsboot <- suppressWarnings(try(stats::nls(formula = formLGauss5p, data = dsetboot, start = lestimpar,
                                                     lower =  c(0, -Inf, -Inf, 0, -Inf), algorithm = "port"),
                                          silent = TRUE))
        } else {
          if (f1 == 0) {
            lestimpar.f0 <- list(b = lestimpar$b, c = lestimpar$c, d = lestimpar$d, e = lestimpar$e)
            nlsboot <- suppressWarnings(try(stats::nls(formula = formLprobit, data = dsetboot, start = lestimpar.f0,
                                                       lower = c(0, -Inf, -Inf, 0), algorithm = "port"),
                                            silent = TRUE))
          } else {
            lestimpar.4p <- list(b = lestimpar$b, d = lestimpar$d, e = lestimpar$e, f = lestimpar$f)
            nlsboot <- suppressWarnings(try(stats::nls(formula = formLGauss4p, data = dsetboot, start = lestimpar.4p,
                                                       lower = c(0, -Inf, 0, -Inf), algorithm = "port"),
                                            silent = TRUE))
          }
        }
        
        if (inherits(nlsboot, "nls")) {
          SDresboot <- sqrt(sum(stats::residuals(nlsboot)^2) / (ndata - nbpari))
          bboot <- stats::coef(nlsboot)["b"]
          if (nbpari == 5 || f1 == 0) {
            cboot <- stats::coef(nlsboot)["c"]
          } else {
            cboot <- stats::coef(nlsboot)["d"]
          }
          dboot <- stats::coef(nlsboot)["d"]
          eboot <- stats::coef(nlsboot)["e"]
          if (f1 == 0) {
            fboot <- 0
          } else {
            fboot <- stats::coef(nlsboot)["f"]
          }
          y0boot <- dboot
          ydosemaxboot <- fLGauss5p(x = dosemax, b = bboot, c = cboot, d = dboot, e = eboot, f = fboot)
          
          if (f1 != 0) {
            xextrboot <-  exp(log(eboot) + (cboot - dboot) * bboot / (fboot * sqrt(2 * pi)))
            yextrboot <-  fLGauss5p(x = xextrboot, b = bboot, c = cboot, d = dboot, e = eboot, f = fboot)
            
            deltapboot <- abs(y0boot) * xdiv100
            deltasdboot <- z * SDresboot
            
            resBMDp <- suppressWarnings(try(
              calcBMD(y0 = y0boot, delta = deltapboot, xext = xextrboot, yext = yextrboot,
                      dosemin = dosemin, dosemax = dosemax, ydosemax = ydosemaxboot,
                      func = fLGauss5pBMR, func_xinlog = fLGauss5pBMR_xinlog,
                      b = bboot, c = cboot, d = dboot, e = eboot, g = fboot,
                      minBMD = minBMD, ratio2switchinlog = ratio2switchinlog),
              silent = TRUE))
            
            resBMDsd <- suppressWarnings(try(
              calcBMD(y0 = y0boot, delta = deltasdboot, xext = xextrboot, yext = yextrboot,
                      dosemin = dosemin, dosemax = dosemax, ydosemax = ydosemaxboot,
                      func = fLGauss5pBMR, func_xinlog = fLGauss5pBMR_xinlog,
                      b = bboot, c = cboot, d = dboot, e = eboot, g = fboot,
                      minBMD = minBMD, ratio2switchinlog = ratio2switchinlog),
              silent = TRUE))
            
            # return a value only if no problem with nls AND uniroot
            if (!inherits(resBMDsd, "try-error") && !inherits(resBMDp, "try-error")) {
              return(list(BMDp = resBMDp$BMD, BMDsd = resBMDsd$BMD))
            }
            
          } else { # if f1 == 0
            ypboot <- y0boot * (1 + xdiv100 * sign(cboot * dboot))
            BMDpboot <- pmax(invLprobit(ypboot, b = bboot, c = cboot, d = dboot, e = eboot), minBMD)
            ysdboot <- y0boot + z * SDresboot * sign(cboot * dboot)
            BMDsdboot <- pmax(invLprobit(ysdboot, b = bboot, c = cboot, d = dboot, e = eboot), minBMD)
            return(list(BMDp = BMDpboot, BMDsd = BMDsdboot))
          } # end if f1 == 0
        } # end if (inherits(nlsboot, "nls"))
      } # end fboot
    }
    ############ END model log Gauss probit ###################
    
    
    ########### Bootstrap iterations on item i ################
    l1 <- lapply(1:niter, fboot)
    nboot.successful <- niter - sum(sapply(l1, is.null))
    if (nboot.successful < niter * tol) {
      # warning(strwrap(paste0("Procedure aborted: the fit only converged for ", nboot.successful,
      #               " iterations during bootstrapping for item ", items[i], ".")), fill = TRUE)
      return(c(NA, NA, NA, NA, nboot.successful))
    } else {
      BMDpbooti <- sapply(l1[!sapply(l1, is.null)], function(z) z$BMDp)
      BMDsdbooti <- sapply(l1[!sapply(l1, is.null)], function(z) z$BMDsd)
      
      BMDpbooti[is.na(BMDpbooti) | BMDpbooti > dosemax] <- Inf
      BMDsdbooti[is.na(BMDsdbooti) | BMDsdbooti > dosemax] <- Inf
      BMDp.CI <- stats::quantile(BMDpbooti, probs = c(prob.lower, prob.upper))
      BMDplower <- BMDp.CI[1]
      BMDpupper <- BMDp.CI[2]
      
      BMDsd.CI <- stats::quantile(BMDsdbooti, probs = c(prob.lower, prob.upper))
      BMDsdlower <- BMDsd.CI[1]
      BMDsdupper <- BMDsd.CI[2]
      
      if (progressbar) {
        utils::setTxtProgressBar(pb, i)
      }
      return(c(BMDsdlower, BMDsdupper, BMDplower, BMDpupper, nboot.successful))
    }
    
  }
  ##### END bootstrap for one item ####################
  
  # Loop on items
  # parallel or sequential computation
  if (parallel != "no") {
    if (parallel == "snow") {
      type <- "PSOCK"
    } else if (parallel == "multicore") {
      type <- "FORK"
    }
    clus <- parallel::makeCluster(ncpus, type = type)
    res <- parallel::parSapply(clus, 1:nitems, bootoneitem)
    parallel::stopCluster(clus)
  } else {
    res <- sapply(1:nitems, bootoneitem)
  }
  
  # close progress bar
  if (progressbar) {
    close(pb)
  }
  
  dres <- as.data.frame(t(res))
  colnames(dres) <- c("BMD.zSD.lower", "BMD.zSD.upper", "BMD.xfold.lower", "BMD.xfold.upper", "nboot.successful")
  
  dres <- cbind(r$res[i.items, ], dres)
  reslist <- list(res = dres,
                  z = z, x = x, tol = tol, niter = niter)
  return(structure(reslist, class = "bmdboot"))
}
############################# END of bmdboot


print.bmdboot <- function(x, ...) {
  if (!inherits(x, "bmdboot")) {
    stop("Use only with 'bmdboot' objects.")
  }
  
  ntot <- nrow(x$res)
  nNA.BMDboot <- sum(x$res$nboot.successful < x$tol * x$niter)
  if (nNA.BMDboot == 0) {
    cat(strwrap(paste0("Bootstrap confidence interval computation was successful on ", ntot, " items among", ntot, ".")), fill = TRUE)
  } else {
    cat(strwrap(paste0("Bootstrap confidence interval computation failed on ", nNA.BMDboot, " items among ", ntot,
                       " due to lack of convergence of the model fit for a fraction of the 
                       bootstrapped samples greater than ", x$tol, ".")), fill = TRUE)
  }
  
  nInf.BMD.zSD.upper <- sum(is.infinite(x$res$BMD.zSD.upper))
  nInf.BMD.xfold.upper <- sum(is.infinite(x$res$BMD.xfold.upper))
  cat(strwrap(paste0("For ", nInf.BMD.zSD.upper, " BMD.zSD values and ", nInf.BMD.xfold.upper,
                     " BMD.xfold values among ", ntot,
                     " at least one bound of the 95 percent confidence interval could not be 
                     computed due to some bootstrapped BMD values not reachable due to model asymptotes 
                     or reached outside the range of tested doses (bounds coded Inf)).")), fill = TRUE)
}


plot.bmdboot <- function(x, BMDtype = c("zSD", "xfold"), remove.infinite = TRUE,
                         by = c("none", "trend", "model", "typology"),
                         CI.col = "blue", BMD_log_transfo = TRUE, ...) {
  
  if (!inherits(x, "bmdboot")) {
    stop("Use only with 'bmdboot' objects.")
  }
  
  BMDtype <- match.arg(BMDtype, c("zSD", "xfold"))
  by <- match.arg(by, c("none", "trend", "model", "typology"))
  
  res <- x$res
  
  if (BMDtype == "zSD") {
    dwithNA <- data.frame(BMD = res$BMD.zSD, BMD.lower = res$BMD.zSD.lower,
                          BMD.upper = res$BMD.zSD.upper)
  } else {
    dwithNA <- data.frame(BMD = res$BMD.xfold, BMD.lower = res$BMD.xfold.lower,
                          BMD.upper = res$BMD.xfold.upper)
  }
  nrow(dwithNA)
  
  if (by == "trend") {
    dwithNA$by <- res$trend
  } else if (by == "model") {
    dwithNA$by <- res$model
  } else if (by == "typology") {
    dwithNA$by <- res$typology
  }
  
  # Remove NA values if needed
  d <- dwithNA[!is.na(dwithNA$BMD) & !is.na(dwithNA$BMD.lower) & !is.na(dwithNA$BMD.upper), ]
  nrow(d)
  
  # remove BMD with infinite lower bounds
  d <- d[is.finite(d$BMD.lower), ]
  nrow(d)
  
  # remove BMD with infinite upper bounds
  if (remove.infinite) {
    # remove BMD with infinite upper values
    d <- d[is.finite(d$BMD.upper), ]
    define.xlim <- FALSE
    nrow(d)
  } else {
    ind.infinite <- !is.finite(d$BMD.upper)
    if (any(ind.infinite)) {
      allBMDval <- c(d$BMD, d$BMD.upper, d$BMD.lower)
      BMDmax <- max(allBMDval[is.finite(allBMDval)])
      BMDlimmax <- BMDmax * 1.5
      d$BMD.upper[ind.infinite] <- BMDlimmax # I did not manage to fix it at a higher value
      # (not plotted by geom_errorbarh)
      define.xlim <- TRUE
    } else {
      define.xlim <- FALSE
    }
  }
  
  nplotted <- nrow(d)
  nremoved <- nrow(dwithNA) - nplotted
  
  if (nremoved > 0) {
    if (remove.infinite) {
      warning(strwrap(prefix = "\n", initial = "\n", paste0(nremoved,
                                                            " BMD values for which lower and upper bounds were coded NA or with lower or upper infinite bounds were removed before plotting.")))
    } else {
      warning(strwrap(prefix = "\n", initial = "\n", paste0(nremoved,
                                                            " BMD values for which lower and upper bounds were coded NA or with lower and upper infinite bounds were removed before plotting.")))
    }
  }
  
  if (by != "none") {
    uniqueby <- unique(d$by)
    n.uniqueby <- length(uniqueby)
    d$ECDF <- rep(0, nplotted) # initialization
    for (i in 1:n.uniqueby) {
      indi <- which(d$by == uniqueby[i])
      ntoti <- length(indi)
      d$ECDF[indi] <- (rank(d$BMD[indi], ties.method = "first") - 0.5) / ntoti
    }
    
    if (!define.xlim) {
      g <- ggplot(data = d, mapping = aes(x = .data$BMD, y = .data$ECDF)) +
        facet_wrap(~ by) +
        geom_errorbarh(aes(xmin = .data$BMD.lower, xmax = .data$BMD.upper), col = CI.col,
                       alpha = 0.5, height = 0) + geom_point()
    } else {
      g <- ggplot(data = d, mapping = aes(x = .data$BMD, y = .data$ECDF)) +
        facet_wrap(~ by) +
        geom_errorbarh(aes(xmin = .data$BMD.lower, xmax = .data$BMD.upper), col = CI.col,
                       alpha = 0.5, height = 0) + geom_point() + xlim(0, BMDlimmax)
      
    }
  } else { # global plot of BMDs
    d$ECDF <- (rank(d$BMD, ties.method = "first") - 0.5) / nplotted
    if (!define.xlim) {
      g <- ggplot(data = d, mapping = aes(x = .data$BMD, y = .data$ECDF)) +
        geom_errorbarh(aes(xmin = .data$BMD.lower, xmax = .data$BMD.upper), col = CI.col,
                       alpha = 0.5,  height = 0) + geom_point()
    } else {
      g <- ggplot(data = d, mapping = aes(x = .data$BMD, y = .data$ECDF)) +
        geom_errorbarh(aes(xmin = .data$BMD.lower, xmax = .data$BMD.upper), col = CI.col,
                       alpha = 0.5,  height = 0) + geom_point() + xlim(0, BMDlimmax)
    }
  }
  
  if (BMD_log_transfo) {
    g <- g + scale_x_log10() + xlab("BMD (in log scale)")
  }
  
  return(g)
}
