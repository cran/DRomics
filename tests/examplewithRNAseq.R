library(DRomics)
# importation and check of RNAseq data and normalization
# with respect to library size and transformation 
# options to put in shiny : transfo.method (2 methods, rlog or vst)
datafilename <- system.file("extdata", "RNAseq_sample.txt", package="DRomics")
(o <- RNAseqdata(datafilename, check = TRUE, transfo.method = "rlog"))
plot(o)

# item selection using the quadratic method
# options to put in shiny : select.method (3 methods), FDR (numerical positive value < 1)
(s_quad <- itemselect(o, select.method = "quadratic", FDR = 0.001))

(s_lin <- itemselect(o, select.method = "linear", FDR = 0.001))
(s_ANOVA <- itemselect(o, select.method = "ANOVA", FDR = 0.001))

# no options in shiny
(f <- drcfit(s_quad, progressbar = TRUE))
f$fitres
plot(f)

# various plot of fitted curves (without data)
curvesplot(f$fitres, xmax = max(f$omicdata$dose), 
           facetby = "model", colorby = "model")

curvesplot(f$fitres, xmax = max(f$omicdata$dose), 
           facetby = "typology")

# plot of selection of curves
curvesplot(f$fitres[f$fitres$trend == "bell", ], xmax = max(f$omicdata$dose), 
           facetby = "id")
curvesplot(f$fitres[f$fitres$trend == "U", ], xmax = max(f$omicdata$dose), 
           facetby = "id")


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
b <- bmdboot(r, niter = 100) # niter should be fixed at least at 1000 to get a reasonable precision
plot(b)