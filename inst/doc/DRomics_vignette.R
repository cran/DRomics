## ----setup, echo=FALSE, message=FALSE, warning=FALSE--------------------------
require(DRomics)
require(ggplot2)
set.seed(1234)
options(digits = 3)
knitr::opts_chunk$set(echo = TRUE,
                      eval = TRUE, 
                      message=FALSE, 
                      warning=FALSE,
                      cache=FALSE,
                      fig.width = 7, 
                      fig.height = 5)

## ----eval = FALSE-------------------------------------------------------------
#  # Installation of required shiny packages
#  install.packages(c("shiny", "shinyBS", "shinycssloaders", "shinyjs", "shinyWidgets", "sortable"))
#  
#  # Launch of the first shiny application DRomics-shiny
#  shiny::runApp(system.file("DRomics-shiny", package = "DRomics"))
#  
#  # Launch of the second shiny application DRomicsInterpreter-shiny
#  shiny::runApp(system.file("DRomicsInterpreter-shiny", package = "DRomics"))

## -----------------------------------------------------------------------------
# Import the text file just to see what will be automatically imported
datafilename <- system.file("extdata", "RNAseq_sample.txt", package = "DRomics")
# datafilename <- "yourchosenname.txt" # for a local file
# Have a look of what information is coded in this file
d <- read.table(file = datafilename, header = FALSE)
nrow(d)
head(d)

## -----------------------------------------------------------------------------
# Load and look at the dataset directly coded as an R object
data(Zhou_kidney_pce)
nrow(Zhou_kidney_pce)
head(Zhou_kidney_pce)

## -----------------------------------------------------------------------------
# Load and look at the data as initially coded
data(zebraf)
str(zebraf)

# Formatting of data for use in DRomics
# 
data4DRomics <- formatdata4DRomics(signalmatrix = zebraf$counts, 
                           dose = zebraf$dose)
# Look at the dataset coded as an R object
nrow(data4DRomics)
head(data4DRomics)

## -----------------------------------------------------------------------------
RNAseqfilename <- system.file("extdata", "RNAseq_sample.txt", package = "DRomics")
# RNAseqfilename <- "yourchosenname.txt" # for a local file

## -----------------------------------------------------------------------------
(o.RNAseq <- RNAseqdata(RNAseqfilename, transfo.method = "vst"))
plot(o.RNAseq, cex.main = 0.8, col = "green")

## -----------------------------------------------------------------------------
microarrayfilename <- system.file("extdata", "transcripto_sample.txt", package = "DRomics")
# microarrayfilename <- "yourchosenname.txt" # for a local file

## -----------------------------------------------------------------------------
(o.microarray <- microarraydata(microarrayfilename, norm.method = "quantile"))
plot(o.microarray, cex.main = 0.8, col = "green")

## -----------------------------------------------------------------------------
metabolofilename <- system.file("extdata", "metabolo_sample.txt", package = "DRomics")
# metabolofilename <- "yourchosenname.txt" # for a local file

## -----------------------------------------------------------------------------
(o.metabolo <- continuousomicdata(metabolofilename))
plot(o.metabolo, col = "green")

## -----------------------------------------------------------------------------
anchoringfilename <- system.file("extdata", "apical_anchoring.txt", package = "DRomics")
# anchoringfilename <- "yourchosenname.txt" # for a local file

## ----fig.width = 7, fig.height = 3--------------------------------------------
(o.anchoring <- continuousanchoringdata(anchoringfilename, backgrounddose = 0.1))
plot(o.anchoring) + theme_bw()

## -----------------------------------------------------------------------------
datafilename <- system.file("extdata", "insitu_RNAseq_sample.txt", package="DRomics")
# Importation of data specifying that observed doses below the background dose
# fixed here to 2e-2 will be considered as null dose to have a control 
(o.insitu <- RNAseqdata(datafilename, backgrounddose = 2e-2, transfo.method = "vst"))
plot(o.insitu)

## -----------------------------------------------------------------------------
# Load of data
data(zebraf)
str(zebraf)

# Look at the design of this dataset
xtabs(~ zebraf$dose + zebraf$batch)

## -----------------------------------------------------------------------------
# Formating of data using the formatdata4DRomics() function
data4DRomics <- formatdata4DRomics(signalmatrix = zebraf$counts, 
                           dose = zebraf$dose)

# Importation of data just to use DRomics functions
# As only raw data will be given to ComBat_seq after
(o <- RNAseqdata(data4DRomics))

# PCA plot to visualize the batch effect
PCAdataplot(o, batch = zebraf$batch) + theme_bw()

## ----results = "hide", message = FALSE----------------------------------------
# Batch effect correction using ComBat_seq{sva}
require(sva)
BECcounts <- ComBat_seq(as.matrix(o$raw.counts), 
                        batch = as.factor(zebraf$batch), 
                        group = as.factor(o$dose)) 

## -----------------------------------------------------------------------------
# Formating of data after batch effect correction
BECdata4DRomics <- formatdata4DRomics(signalmatrix = BECcounts, 
                                   dose = o$dose)
o.BEC <- RNAseqdata(BECdata4DRomics)

# PCA plot after batch effect correction
PCAdataplot(o.BEC, batch = zebraf$batch) + theme_bw()

## -----------------------------------------------------------------------------
(s_quad <- itemselect(o.microarray, select.method = "quadratic", FDR = 0.01))

## ----results = "hide"---------------------------------------------------------
require(VennDiagram)
s_lin <- itemselect(o.microarray, select.method = "linear", FDR = 0.01)
index_quad <- s_quad$selectindex
index_lin <- s_lin$selectindex
plot(c(0,0), c(1,1), type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
draw.pairwise.venn(area1 = length(index_quad), area2 = length(index_lin), 
                   cross.area = length(which(index_quad %in% index_lin)), 
                   category = c("quadratic trend test", "linear trend test"),
                   cat.col=c("cyan3", "darkorange1"), col=c("black", "black"), 
                   fill = c("cyan3", "darkorange1"), lty = "blank", cat.pos = c(1,11))


## -----------------------------------------------------------------------------
(f <- drcfit(s_quad, progressbar = FALSE))

## -----------------------------------------------------------------------------
head(f$fitres)

## -----------------------------------------------------------------------------
plot(f) 

## -----------------------------------------------------------------------------
targetitems <- c("88.1", "1", "3", "15")
targetplot(targetitems, f = f)

## ----fig.width = 7, fig.height = 5--------------------------------------------
plot(f, plot.type = "dose_residuals")

## ----echo = FALSE, results = "hide", fig.height=8, fig.width = 8--------------
par(mfrow = c(4,4), mar = c(0,0,0,0), xaxt = "n", yaxt = "n")
x <- seq(0,10, length.out = 50)
# linear
plot(x, DRomics:::flin(x, b = 1, d = 1), type = "l", lwd = 2, col = "red")
legend("topleft", legend = "linear, b > 0", bty = "n")
plot(x, DRomics:::flin(x, b = -1, d = 1), type = "l", lwd = 2, col = "red")
legend("bottomleft", legend = "linear, b < 0", bty = "n")

# expo
plot(x, DRomics:::fExpo(x, b = 1, d = 1, e = 3), type = "l", lwd = 2, col = "red")
legend("topleft", legend = "exponential, e > 0 and b > 0", bty = "n")
plot(x, DRomics:::fExpo(x, b = -1, d = 1, e = 3), type = "l", lwd = 2, col = "red")
legend("bottomleft", legend = "exponential, e > 0 and b < 0", bty = "n")
plot(x, DRomics:::fExpo(x, b = 1, d = 1, e = -3), type = "l", lwd = 2, col = "red")
legend("topright", legend = "exponential, e < 0 and b > 0", bty = "n")
plot(x, DRomics:::fExpo(x, b = -1, d = 1, e = -3), type = "l", lwd = 2, col = "red")
legend("bottomright", legend = "exponential, e < 0 and b < 0", bty = "n")

# Hill
plot(x, DRomics:::fHill(x, b = 10, c = 3, d = 1, e = 3), type = "l", lwd = 2, col = "red")
legend("bottomright", legend = "Hill, c > d", bty = "n")
plot(x, DRomics:::fHill(x, b = 10, c = 1, d = 3, e = 3), type = "l", lwd = 2, col = "red")
legend("topright", legend = "Hill, c < d", bty = "n")

# Gauss-probit
plot(x, DRomics:::fGauss5p(x, b = 2, c = 3, d = 1, e = 3, f = 2), type = "l", lwd = 2, col = "red")
legend("bottomright", legend = "Gauss-probit, c > d, f > 0", bty = "n")
plot(x, DRomics:::fGauss5p(x, b = 2, c = 1, d = 3, e = 3, f = 2), type = "l", lwd = 2, col = "red")
legend("topright", legend = "Gauss-probit, c < d, f > 0", bty = "n")
plot(x, DRomics:::fGauss5p(x, b = 2, c = 3, d = 1, e = 3, f = -2), type = "l", lwd = 2, col = "red")
legend("bottomright", legend = "Gauss-probit, c > d, f < 0", bty = "n")
plot(x, DRomics:::fGauss5p(x, b = 2, c = 1, d = 3, e = 3, f = -2), type = "l", lwd = 2, col = "red")
legend("topright", legend = "Gauss-probit, c < d, f < 0", bty = "n")

# LGauss-probit
x <- seq(0,100, length.out = 50)
plot(x, DRomics:::fLGauss5p(x, b = 0.5, c = 3, d = 1, e = 20, f = 4), type = "l", lwd = 2, col = "red")
legend("bottomright", legend = "log-Gauss-probit, c > d, f > 0", bty = "n")
plot(x, DRomics:::fLGauss5p(x, b = 0.5, c = 1, d = 3, e = 20, f = 4), type = "l", lwd = 2, col = "red")
legend("topright", legend = "log-Gauss-probit, c < d, f > 0", bty = "n")
plot(x, DRomics:::fLGauss5p(x, b = 0.5, c = 3, d = 1, e = 20, f = -4), type = "l", lwd = 2, col = "red")
legend("bottomright", legend = "log-Gauss-probit, c > d, f < 0", bty = "n")
plot(x, DRomics:::fLGauss5p(x, b = 0.5, c = 1, d = 3, e = 20, f = -4), type = "l", lwd = 2, col = "red")
legend("topright", legend = "log-Gauss-probit, c < d, f < 0", bty = "n")


## ----echo = FALSE,  fig.height = 4, fig.width = 7, results = "hide", out.width="80%"----
par(mar = c(0.1, 0.1, 0.1, 0.1))
datafilename <- system.file("extdata", "apical_anchoring.txt", package = "DRomics")
o_ls <- continuousanchoringdata(datafilename, check = TRUE, backgrounddose = 0.1)
s_ls <- itemselect(o_ls)
f_ls <- drcfit(s_ls)
growth <- f_ls$fitres[1,]
#plot(f)
plot(o_ls$dose, o_ls$data[1,], xlab = "dose", ylab = "signal", 
     pch = 16, xlim = c(0, 30), ylim = c(-20, 80))
xfin <- seq(0, 80, length.out = 100)
#plot(x, x+100, ylim = c(0, 7), xlab = "dose", ylab = "signal")
valb <- growth$b; valc <- growth$c; vald <- growth$d 
vale <- growth$e; valf <- growth$f
ytheo <- DRomics:::fGauss5p(xfin, valb, valc, vald, vale, valf)
lines(xfin, ytheo, col = "red", lwd = 2)

# Ajout de lois normales en vertical
doseu <- sort(unique(o_ls$dose))
ytheou <- DRomics:::fGauss5p(doseu, valb, valc, vald, vale, valf)
sy <- growth$SDres
npts <- 50 # nb de points par normale
coefsurx <- 12
tracenormale <- function(indice)
{
  x <- doseu[indice]
  my <- ytheou[indice]
  yplot <- seq(my - 2*sy, my+2*sy, length.out = npts)
  xplot <- dnorm(yplot, mean = my, sd = sy)
  lines(coefsurx*xplot+x, yplot, col = "blue", lwd = 2)
  segments(x, my - 2*sy, x, my + 2*sy, lty = 3, col = "blue")
}
sapply(1:7, tracenormale)


## -----------------------------------------------------------------------------
(r <- bmdcalc(f, z = 1, x = 10))

## -----------------------------------------------------------------------------
head(r$res)

## -----------------------------------------------------------------------------
plot(r, BMDtype = "zSD", plottype = "ecdf") + theme_bw() 

## -----------------------------------------------------------------------------
bmdplotwithgradient(r$res, BMDtype = "zSD",
                    facetby = "trend", 
                    shapeby = "model",
                    line.size = 1.2,
                    scaling = TRUE) 

## -----------------------------------------------------------------------------
(b <- bmdboot(r, niter = 50, progressbar = FALSE))

## -----------------------------------------------------------------------------
head(b$res)

## ----fig.height = 3-----------------------------------------------------------
# Plot of BMDs with no filtering
subres <- bmdfilter(b$res, BMDfilter = "none")
bmdplot(subres, BMDtype = "xfold", point.size = 2, point.alpha = 0.4, 
        add.CI = TRUE, line.size = 0.4) + theme_bw()

# Plot of items with defined BMD point estimate
subres <- bmdfilter(b$res, BMDtype = "xfold", BMDfilter = "definedBMD")
bmdplot(subres, BMDtype = "xfold", point.size = 2, point.alpha = 0.4, 
        add.CI = TRUE, line.size = 0.4) + theme_bw()

# Plot of items with defined BMD point estimate and CI bounds
subres <- bmdfilter(b$res, BMDtype = "xfold", BMDfilter = "definedCI")
bmdplot(subres, BMDtype = "xfold", point.size = 2, point.alpha = 0.4, 
        add.CI = TRUE, line.size = 0.4) + theme_bw()

# Plot of items with finite BMD point estimate and CI bounds
subres <- bmdfilter(b$res, BMDtype = "xfold", BMDfilter = "finiteCI") 
bmdplot(subres, BMDtype = "xfold", point.size = 2, point.alpha = 0.4, 
        add.CI = TRUE, line.size = 0.4) + theme_bw()

## -----------------------------------------------------------------------------
# If you do not want to add the confidence intervals just replace b
# the output of bmdboot() by r the output of bmdcalc()
plot(f, BMDoutput = b) 

## -----------------------------------------------------------------------------
tested.doses <- unique(f$omicdata$dose)
g <- curvesplot(r$res, addBMD = TRUE, xmax = max(tested.doses), colorby = "trend",
           line.size = 0.8, line.alpha = 0.3, point.size = 2, point.alpha = 0.6) +
  geom_vline(xintercept = tested.doses, linetype = 2) + theme_bw()
print(g)

## ----eval = FALSE-------------------------------------------------------------
#  if (require(plotly))
#  {
#    ggplotly(g)
#  }
#  

## -----------------------------------------------------------------------------
str(b$res)

## -----------------------------------------------------------------------------
# code to import the file for this example stored in our package
resfilename <- system.file("extdata", "triclosanSVmetabres.txt", package = "DRomics")
res <- read.table(resfilename, header = TRUE, stringsAsFactors = TRUE)

# to see the first lines of this data frame
head(res)

## -----------------------------------------------------------------------------
# code to import the file for this example in our package
annotfilename <- system.file("extdata", "triclosanSVmetabannot.txt", package = "DRomics")
# annotfilename <- "yourchosenname.txt" # for a local file
annot <- read.table(annotfilename, header = TRUE, stringsAsFactors = TRUE)

# to see the first lines of this data frame
head(annot)

## -----------------------------------------------------------------------------
# Merging
extendedres <- merge(x = res, y = annot, by.x = "id", by.y = "metab.code")

# to see the first lines of the merged data frame
head(extendedres)

## -----------------------------------------------------------------------------
bmdplot(extendedres, BMDtype = "zSD", add.CI = TRUE, 
                    facetby = "path_class", 
                    colorby = "trend") + theme_bw()

## ----eval = FALSE-------------------------------------------------------------
#  ecdfplotwithCI(variable = extendedres$BMD.zSD,
#                 CI.lower = extendedres$BMD.zSD.lower,
#                 CI.upper = extendedres$BMD.zSD.upper,
#                 by = extendedres$path_class,
#                 CI.col = extendedres$trend) + labs(col = "trend")

## -----------------------------------------------------------------------------
bmdplotwithgradient(extendedres, BMDtype = "zSD",
                    scaling = TRUE, 
                    facetby = "path_class", 
                    shapeby = "trend") 

## -----------------------------------------------------------------------------
extendedres_lipid <- extendedres[extendedres$path_class == "Lipid metabolism",] 

bmdplotwithgradient(extendedres_lipid, BMDtype = "zSD",
                    scaling = TRUE,
                    facetby = "path_class", 
                    add.label = TRUE,
                    xmin = 0, xmax = 6,
                    label.size = 3,
                    line.size = 2) 


## -----------------------------------------------------------------------------
sensitivityplot(extendedres, BMDtype = "zSD",
                group = "path_class",
                BMDsummary = "first.quartile") + theme_bw()

## -----------------------------------------------------------------------------
sensitivityplot(extendedres, BMDtype = "zSD",
                group = "path_class",
                BMDsummary = "median.and.IQR") + theme_bw()

## -----------------------------------------------------------------------------
psens <- sensitivityplot(extendedres, BMDtype = "zSD",
                         group = "path_class",
                         BMDsummary = "first.quartile")

psens + 
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_text(aes(label = paste("   ", psens$data$groupby, "   ")), 
            size = 3, hjust = "inward")

## -----------------------------------------------------------------------------
trendplot(extendedres, group = "path_class") + theme_bw()

## -----------------------------------------------------------------------------
# Plot of all the scaled dose-reponse curves split by path class
curvesplot(extendedres, facetby = "path_class", scaling = TRUE, npoints = 100, 
           colorby = "trend", xmin = 0, xmax = 6.5, 
           addBMD = TRUE) + theme_bw()

## -----------------------------------------------------------------------------
# Plot of the unscaled dose-reponses for one chosen group, split by metabolite
LMres <- extendedres[extendedres$path_class == "Lipid metabolism", ]
curvesplot(LMres, facetby = "id", npoints = 100, 
           point.size = 1.5, line.size = 1,
           colorby = "trend",
           xmin = 0, xmax = 6.5) + theme_bw() 

## -----------------------------------------------------------------------------
# 1. Import the data frame with DRomics results to be used
contigresfilename <- system.file("extdata", "triclosanSVcontigres.txt", package = "DRomics")
contigres <- read.table(contigresfilename, header = TRUE, stringsAsFactors = TRUE)

# 2. Import the data frame with biological annotation (or any other descriptor/category 
# you want to use, here KEGG pathway classes) 
contigannotfilename <- system.file("extdata", "triclosanSVcontigannot.txt", package = "DRomics")
# contigannotfilename <- "yourchosenname.txt" # for a local file
contigannot <- read.table(contigannotfilename, header = TRUE, stringsAsFactors = TRUE)

# 3. Merging of both previous data frames   
contigextendedres <- merge(x = contigres, y = contigannot, by.x = "id", by.y = "contig")
# to see the first lines of the data frame
head(contigextendedres)

## -----------------------------------------------------------------------------
metabextendedres <- extendedres

## -----------------------------------------------------------------------------
extendedres <- rbind(metabextendedres, contigextendedres)
extendedres$explevel <- factor(c(rep("metabolites", nrow(metabextendedres)),
                              rep("contigs", nrow(contigextendedres))))
# to see the first lines of the data frame
head(extendedres)

## -----------------------------------------------------------------------------
(t.pathways <- table(extendedres$path_class, extendedres$explevel)) 
original.par <- par()
par(las = 2, mar = c(4,13,1,1))
barplot(t(t.pathways), beside = TRUE, horiz = TRUE, 
        cex.names = 0.7, legend.text = TRUE, 
        main = "Frequencies of pathways")
par(original.par)

## -----------------------------------------------------------------------------
unique.items <- unique(extendedres$id)
ggplot(extendedres[match(unique.items, extendedres$id), ], aes(x = BMD.zSD, color = explevel)) +
      stat_ecdf(geom = "step") + ylab("ECDF") + theme_bw()

## -----------------------------------------------------------------------------
# BMD ECDF plot split by molecular level, after removing items redundancy
bmdplot(extendedres[match(unique.items, extendedres$id), ], BMDtype = "zSD", 
                    facetby = "explevel", point.alpha = 0.4) + theme_bw()

# BMD ECDF plot colored by molecular level and split by path class
bmdplot(extendedres, BMDtype = "zSD", 
                    facetby = "path_class", 
                    colorby = "explevel",
                    point.alpha = 0.4) + 
              labs(col = "molecular level") + theme_bw()

## -----------------------------------------------------------------------------
# Preliminary optional alphabetic ordering of path_class groups
extendedres$path_class <- factor(extendedres$path_class, 
                levels = sort(levels(extendedres$path_class), decreasing = TRUE))

# Trend plot
trendplot(extendedres, group = "path_class", facetby = "explevel") +
  theme_bw()

## -----------------------------------------------------------------------------
sensitivityplot(extendedres, BMDtype = "zSD",
                group = "path_class", colorby = "explevel",
                BMDsummary = "first.quartile") + theme_bw()

## -----------------------------------------------------------------------------
selectedres <- selectgroups(extendedres, 
                         group = "path_class",
                         explev = "explevel",
                         BMDmax = 0.75,
                         BMDtype = "zSD", 
                         BMDsummary = "first.quartile",
                         nitems = 3,
                         keepallexplev = TRUE)

# BMDplot on this selection
bmdplot(selectedres, BMDtype = "zSD", add.CI = TRUE,
                    facetby = "path_class", facetby2 = "explevel",
                    colorby = "trend") + theme_bw()

## -----------------------------------------------------------------------------
# Manual selection of groups on which to focus
chosen_path_class <- c("Nucleotide metabolism", 
                       "Membrane transport", 
                       "Lipid metabolism", 
                       "Energy metabolism")
selectedres2 <- extendedres[extendedres$path_class %in% chosen_path_class, ]
bmdplotwithgradient(selectedres2, BMDtype = "zSD", scaling = TRUE,
               facetby = "path_class", facetby2 = "explevel")

## -----------------------------------------------------------------------------
# Plot of the unscaled dose-response curves for the "lipid metabolism" path class
# using transparency to get an idea of density of curves with the shame shape
LMres <- extendedres[extendedres$path_class == "Lipid metabolism", ]
curvesplot(LMres, facetby = "explevel", free.y.scales = TRUE, npoints = 100, 
           addBMD = TRUE, line.alpha = 0.4, line.size = 1, colorby = "trend",
           xmin = 0, xmax = 6.5) + labs(col = "DR trend") + theme_bw()

# Plot of the scaled dose-response curves for previously chosen path classes
curvesplot(selectedres2, scaling = TRUE,
           facetby = "path_class", facetby2 = "explevel",
           addBMD = TRUE, npoints = 100, line.size = 1, line.alpha = 0.4,
           colorby = "trend",
           xmin = 0, xmax = 6.5) + labs(col = "DR trend") + theme_bw()

