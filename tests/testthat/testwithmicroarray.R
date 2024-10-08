test_that("Test DRomics on microarray data", {
  skip_on_cran()
  
  # importation and check of data and normalization if needed
  # options to put in shiny : norm.method (4 methods)
  ## sample of the transcripto data set
  datafilename <- system.file("extdata", "transcripto_very_small_sample.txt", package="DRomics")
  (o <- microarraydata(datafilename, check = TRUE, norm.method = "cyclicloess"))
  plot(o)
  (o.2 <- microarraydata(datafilename, check = TRUE, norm.method = "none"))
  (o.3 <- microarraydata(datafilename, check = TRUE, norm.method = "quantile"))
  (o.4 <- microarraydata(datafilename, check = TRUE, norm.method = "scale"))
  
  plot(o.2)
  plot(o.3)
  plot(o.4)
  
  # item selection using the quadratic method
  # options to put in shiny : select.method (3 methods), FDR (numerical positive value < 1)
  (s_quad <- itemselect(o, select.method = "quadratic", FDR = 0.001))
  
  (s_lin <- itemselect(o, select.method = "linear", FDR = 0.001))
  (s_ANOVA <- itemselect(o, select.method = "ANOVA", FDR = 0.001))
  
  # fit of dose response models and choice of the best fit for each item
  # no options in shiny
  (f <- drcfit(s_quad, progressbar = TRUE))
  f$fitres
  f$unfitres
  plot(f)
  # Alternative plots
  # with a chosen number of first items
  plot(f, items = 12) 
  # with chosen items in a specified order
  plot(f, items = c("15", "27.1", "7.1"))
  # residual plots
  plot(f, items = 12, plot.type = "fitted_residuals") 
  plot(f, items = 12, plot.type = "dose_residuals") 
  # plot with dose in log
  plot(f, items = 12, plot.type = "dose_fitted", dose_log_transfo = FALSE) 
  plot(f, items = 12, plot.type = "dose_residuals", dose_log_transfo = FALSE) 
  plot(f, items = 12, plot.type = "fitted_residuals", dose_log_transfo = FALSE) 
  
  # evaluate the impact of preventsfitsoutofrange
  datafilename <- system.file("extdata", "transcripto_sample.txt", package="DRomics")
  (o1 <- microarraydata(datafilename, check = TRUE, norm.method = "cyclicloess"))
  (s_quad1 <- itemselect(o1, select.method = "quadratic", FDR = 0.001))
  
  (f1 <- drcfit(s_quad1, 
                preventsfitsoutofrange = FALSE,
                enablesfequal0inGP  = FALSE,
                enablesfequal0inLGP  = FALSE,
                progressbar = TRUE))
  (f1bis <- drcfit(s_quad1, 
                   preventsfitsoutofrange = TRUE,
                   enablesfequal0inGP  = FALSE,
                   enablesfequal0inLGP  = FALSE,
                   progressbar = TRUE))
  (f1ter <- drcfit(s_quad1, 
                   preventsfitsoutofrange = TRUE,
                   enablesfequal0inGP  = TRUE,
                   enablesfequal0inLGP  = TRUE,
                   progressbar = TRUE))
  
  (idremovedinf1bis <- f1$fitres$id[!is.element(f1$fitres$id, f1bis$fitres$id)])
  
  (idchanged <- f1bis$fitres$id[which(f1bis$fitres$model != f1ter$fitres$model | 
                                        f1bis$fitres$f != f1ter$fitres$f)])
  # no impact on this dataset  
  
  
  # calculation of benchmark doses
  # options in shiny : z (numerical positive value), x (numerical positive value : percentage)
  (r <- bmdcalc(f, z = 1, x = 10))
  (r.2 <- bmdcalc(f, z = 2, x = 50))
  
  # plot of BMD
  # options in shiny : BMDtype (2 possibilities), plottype (3 possibilities), by (3 possibilities)
  # hist.bins (integer for hist only)
  plot(r, BMDtype = "zSD", plottype = "ecdf", by = "none") 
  plot(r, BMDtype = "xfold", plottype = "ecdf", by = "none") 
  
  plot(r, plottype = "hist", by = "none") 
  plot(r, plottype = "hist", by = "none", hist.bins = 10) 
  plot(r, plottype = "density", by = "none") 
  
  plot(r, plottype = "hist", by = "trend", hist.bins = 10) 
  plot(r, plottype = "hist", by = "model", hist.bins = 10) 
  plot(r, plottype = "hist", by = "typology", hist.bins = 10) 
  
  # Calculation of confidence intervals on BMDs by Bootstrap
  # niter <- 1000
  niter <- 10
  b <- bmdboot(r, niter = niter) # niter should be fixed at least at 1000 to get a reasonable precision
  plot(b) 
  
})