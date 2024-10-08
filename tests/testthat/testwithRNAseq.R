test_that("Test DRomics on RNAseq data", {
  skip_on_cran()
  
  # importation and check of RNAseq data and normalization
  # with respect to library size and transformation 
  # options to put in shiny : transfo.method (2 methods, rlog or vst)
  datafilename <- system.file("extdata", "RNAseq_sample.txt", package="DRomics")
  # small data set 'less than 1000 items (999)
  (o.vst <- RNAseqdata(datafilename, check = TRUE, transfo.method = "vst"))
  plot(o.vst)
  
  plot(o.vst, range = 1.5) # boxplot visualizing outliers
  (o.vst.notblind <- RNAseqdata(datafilename, check = TRUE, transfo.method = "vst",
                                transfo.blind = FALSE))
  plot(o.vst.notblind)
  
  (o.rlog <- RNAseqdata(datafilename, check = TRUE, transfo.method = "rlog"))
  plot(o.rlog)
  
  (o.rlog.notblind <- RNAseqdata(datafilename, check = TRUE, transfo.method = "rlog",
                                 transfo.blind = FALSE))
  plot(o.rlog.notblind)
  
  data(Zhou_kidney_pce)
  
  # variance stabilizing tranformation
  (o1 <- RNAseqdata(Zhou_kidney_pce, check = TRUE, transfo.method = "vst"))
  plot(o1)
  
  # regularized logarithm
  (o2 <- RNAseqdata(Zhou_kidney_pce, check = TRUE, transfo.method = "rlog"))
  plot(o2)
  
  # variance stabilizing tranformation (not blind to the experimental design)
  (o3 <- RNAseqdata(Zhou_kidney_pce, check = TRUE, transfo.method = "vst",
                    transfo.blind = FALSE))
  plot(o3)
  
  # regularized logarithm (not blind to the experimental design)
  (o4 <- RNAseqdata(Zhou_kidney_pce, check = TRUE, transfo.method = "rlog",
                    transfo.blind = FALSE))
  plot(o4)
  
  # item selection using the quadratic method
  # options to put in shiny : select.method (3 methods), FDR (numerical positive value < 1)
  o <- o.rlog
  (s_quad <- itemselect(o, select.method = "quadratic", FDR = 0.01))
  (s_lin <- itemselect(o, select.method = "linear", FDR = 0.01))
  (s_ANOVA <- itemselect(o, select.method = "ANOVA", FDR = 0.01))
  
  (f <- drcfit(s_quad, progressbar = TRUE))
  f$fitres
  plot(f)
  r <- bmdcalc(f)
  
  # various plot of fitted curves (without data)
  curvesplot(r$res, xmax = max(r$omicdata$dose), 
             facetby = "model", colorby = "model")
  
  curvesplot(r$res, xmax = max(r$omicdata$dose), 
             facetby = "typology")
  # plot of selection of curves
  curvesplot(r$res[r$res$trend == "bell", ], xmax = max(r$omicdata$dose), 
             facetby = "id")
  
  # evaluate the impact of preventsfitsoutofrange, enablesfequal0inGP, enablesfequal0inlGP
  data(Zhou_kidney_pce)
  (o1 <- RNAseqdata(Zhou_kidney_pce, check = TRUE, transfo.method = "rlog"))
  s_quad1 <- itemselect(o1, select.method = "quadratic", FDR = 0.1)
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
  targetplot(items = idremovedinf1bis, f1, dose_log_transfo = FALSE) 
  
  # (idchanged <- f1bis$fitres$id[which(f1bis$fitres$model != f1ter$fitres$model | 
  #                                       f1bis$fitres$f != f1ter$fitres$f)])
  # targetplot(items = idchanged, f1bis, dose_log_transfo = TRUE)
  # targetplot(items = idchanged, f1ter, dose_log_transfo = TRUE)
  
  # f1bis$fitres[f1bis$fitres$id %in% idchanged, ]
  # f1ter$fitres[f1ter$fitres$id %in% idchanged, ]
  
  
  # calculation of benchmark doses
  # options in shiny : z (numerical positive value), x (numerical positive value : percentage)
  (r <- bmdcalc(f, z = 1, x = 10))
  (r.2 <- bmdcalc(f, z = 2, x = 50))
  
  # plot of BMD
  # options in shiny : BMDtype (2 possibilities), plottype (3 possibilities), by (3 possibilities)
  # hist.bins (integer for hist only)
  plot(r, BMDtype = "zSD", plottype = "ecdf", by = "none") 
  
  plot(r, BMDtype = "xfold", plottype = "ecdf", by = "none") 
  
  plot(r, plottype = "hist", by = "none", hist.bins = 10) 
  plot(r, plottype = "density", by = "none") 
  
  plot(r, plottype = "hist", by = "trend", hist.bins = 10) 
  
  # Calculation of confidence intervals on BMDs by Bootstrap
  # niter <- 1000
  niter <- 10
  b <- bmdboot(r, niter = niter) # niter should be fixed at least at 1000 to get a reasonable precision
  plot(b)
  
  data(Zhou_kidney_pce)
  
  # exploration of an extreme case (BMD at 0)
  d <- Zhou_kidney_pce
  (o <- RNAseqdata(d))
  plot(o)
  (s <- itemselect(o, select.method = "quadratic", FDR = 0.01))
  (f <- drcfit(s, progressbar = TRUE))
  head(f$fitres)
  
  r <- bmdcalc(f, z = 1)
  plot(r) 
  
  bmdplotwithgradient(r$res, BMDtype = "zSD") 
  bmdplotwithgradient(r$res, BMDtype = "zSD", BMD_log_transfo = FALSE) 
  
  # no more 0 BMD values using argument minBMD 
  # res0 <- r$res[r$res$BMD.zSD == 0, ]
  # curvesplot(res0, xmin =0.0000000001, xmax = max(f$omicdata$dose), 
  #            colorby = "model", dose_log_transfo = TRUE)
  # plot(f, items = r$res[r$res$BMD.zSD == 0, ]$id)
  # plot(f, items = r$res[r$res$BMD.zSD == 0, ]$id, dose_log_trans = TRUE)
})
