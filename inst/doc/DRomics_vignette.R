## ----setup, echo=FALSE, message=FALSE, warning=FALSE--------------------------
# output:
#  html_vignette
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

## -----------------------------------------------------------------------------
# Load and look at the first line of the R object
data(Zhou_kidney_pce)
nrow(Zhou_kidney_pce)
head(Zhou_kidney_pce)

# Import the text file just to see what will be automatically imported
datafilename <- system.file("extdata", "RNAseq_sample.txt", package = "DRomics")
# for your local file datafilename choose the name: "yourchosenname.txt"
d <- read.table(file = datafilename, header = FALSE)
nrow(d)
head(d)

## -----------------------------------------------------------------------------
RNAseqfilename <- system.file("extdata", "RNAseq_sample.txt", package = "DRomics")

## -----------------------------------------------------------------------------
(o.RNAseq <- RNAseqdata(RNAseqfilename, transfo.method = "vst"))
plot(o.RNAseq, cex.main = 0.8, col = "green", range = 1e6)

## -----------------------------------------------------------------------------
microarrayfilename <- system.file("extdata", "transcripto_sample.txt", package = "DRomics")

## -----------------------------------------------------------------------------
(o.microarray <- microarraydata(microarrayfilename, norm.method = "quantile"))
plot(o.microarray, cex.main = 0.8, col = "green", range = 1e6)


## -----------------------------------------------------------------------------
metabolofilename <- system.file("extdata", "metabolo_sample.txt", package = "DRomics")

## -----------------------------------------------------------------------------
(o.metabolo <- continuousomicdata(metabolofilename))
plot(o.metabolo, col = "green", range = 1e6)

## -----------------------------------------------------------------------------
anchoringfilename <- system.file("extdata", "apical_anchoring.txt", package = "DRomics")

## -----------------------------------------------------------------------------
(o.anchoring <- continuousanchoringdata(anchoringfilename))
plot(o.anchoring)

## -----------------------------------------------------------------------------
(s_quad <- itemselect(o.microarray, select.method = "quadratic", FDR = 0.01))

## -----------------------------------------------------------------------------
(f <- drcfit(s_quad, progressbar = FALSE))

## -----------------------------------------------------------------------------
head(f$fitres, 10)

## -----------------------------------------------------------------------------
plot(f) 

## -----------------------------------------------------------------------------
targetitems <- c("88.1", "1", "3", "15")
targetplot(targetitems, f = f)

## ---- fig.width = 7, fig.height = 5-------------------------------------------
plot(f, plot.type = "dose_residuals")

## -----------------------------------------------------------------------------
(r <- bmdcalc(f, z = 1, x = 10))

## -----------------------------------------------------------------------------
head(r$res, 10)

## -----------------------------------------------------------------------------
plot(r, BMDtype = "zSD", plottype = "ecdf") 

## -----------------------------------------------------------------------------
plot(r, BMDtype = "zSD", plottype = "density", by = "trend") 

## -----------------------------------------------------------------------------
bmdplotwithgradient(r$res, BMDtype = "zSD",
                    facetby = "trend", 
                    shapeby = "model",
                    line.size = 1.2) + labs(shape = "model")

## -----------------------------------------------------------------------------
(b <- bmdboot(r, niter = 50, progressbar = FALSE))

## -----------------------------------------------------------------------------
head(b$res, 10)

## -----------------------------------------------------------------------------
plot(b, BMDtype = "zSD", by = "trend") 

## -----------------------------------------------------------------------------
# code to import the file for this example in our package
resfilename <- system.file("extdata", "triclosanSVmetabres.txt", package = "DRomics")
res <- read.table(resfilename, header = TRUE, stringsAsFactors = TRUE)

# to see the structure of this file
str(res)

## -----------------------------------------------------------------------------
# code to import the file for this example in our package
annotfilename <- system.file("extdata", "triclosanSVmetabannot.txt", package = "DRomics")
annot <- read.table(annotfilename, header = TRUE, stringsAsFactors = TRUE)

# to see the structure of this file
str(annot)

## -----------------------------------------------------------------------------
annotres <- merge(x = res, y = annot, by.x = "id", by.y = "metab.code")
head(annotres)

## -----------------------------------------------------------------------------
bmdplot(annotres, BMDtype = "zSD", add.CI = TRUE, 
                    facetby = "path_class", 
                    colorby = "trend") + labs(col = "trend")

## ---- eval = FALSE------------------------------------------------------------
#  ecdfplotwithCI(variable = annotres$BMD.zSD,
#                 CI.lower = annotres$BMD.zSD.lower,
#                 CI.upper = annotres$BMD.zSD.upper,
#                 by = annotres$path_class,
#                 CI.col = annotres$trend) + labs(col = "trend")

## -----------------------------------------------------------------------------
bmdplotwithgradient(annotres, BMDtype = "zSD",
                    facetby = "path_class", 
                    shapeby = "trend") + labs(shape = "trend")

## -----------------------------------------------------------------------------
bmdplotwithgradient(annotres, BMDtype = "zSD",
                    facetby = "path_class", 
                    add.label = TRUE) +
  theme(strip.text.x = element_text(size = 6))


## -----------------------------------------------------------------------------
annotres_lipid <- annotres[annotres$path_class == "Lipid metabolism",] 

bmdplotwithgradient(annotres_lipid, BMDtype = "zSD",
                    facetby = "path_class", 
                    add.label = TRUE,
                    limits4colgradient = c(-0.8, 0.8),
                    xmin = 0, xmax = 6.5,
                    label.size = 3) 


## -----------------------------------------------------------------------------
sensitivityplot(annotres, BMDtype = "zSD",
                    group = "path_class",
                 BMDsummary = "first.quartile") 

## -----------------------------------------------------------------------------
trendplot(annotres, group = "path_class") 

## -----------------------------------------------------------------------------
LMres <- annotres[annotres$path_class == "Lipid metabolism", ]
curvesplot(LMres, facetby = "id", npoints = 100, line.size = 1,
           colorby = "trend",
           xmin = 0, xmax = 8) + labs(col = "trend")

## -----------------------------------------------------------------------------
# 1. Import the dataframe with DRomics results to be used
contigresfilename <- system.file("extdata", "triclosanSVcontigres.txt", package = "DRomics")
contigres <- read.table(contigresfilename, header = TRUE, stringsAsFactors = TRUE)
str(contigres)

# 2. Import the dataframe with functional annotation (or any other descriptor/category 
# you want to use, here KEGG pathway classes) 
contigannotfilename <- system.file("extdata", "triclosanSVcontigannot.txt", package = "DRomics")
contigannot <- read.table(contigannotfilename, header = TRUE, stringsAsFactors = TRUE)
str(contigannot)

# 3. Merging of both previous dataframes   
contigextendedres <- merge(x = contigres, y = contigannot, by.x = "id", by.y = "contig")
# to see the structure of this dataframe
str(contigextendedres)

## -----------------------------------------------------------------------------
metabextendedres <- annotres

## -----------------------------------------------------------------------------
extendedres <- rbind(metabextendedres, contigextendedres)
extendedres$level <- factor(c(rep("metabolites", nrow(metabextendedres)),
                              rep("contigs", nrow(contigextendedres))))
str(extendedres)

## -----------------------------------------------------------------------------
(t.pathways <- table(extendedres$path_class, extendedres$level)) 
original.par <- par()
par(las = 2, mar = c(4,13,1,1))
barplot(t(t.pathways), beside = TRUE, horiz = TRUE, 
        cex.names = 0.7, legend.text = TRUE, 
        main = "Frequencies of pathways")
par(original.par)

## -----------------------------------------------------------------------------
(t.prop.pathways <- prop.table(t.pathways, margin = 2)) 
original.par <- par()
par(las = 2, mar = c(4,13,1,1))
barplot(t(t.prop.pathways), beside = TRUE, horiz = TRUE, 
        cex.names = 0.7, legend.text = TRUE, 
        main = "Proportion of pathways")
par(original.par)

## -----------------------------------------------------------------------------
if (require(ggplot2))
{
   ggplot(extendedres, aes(x = BMD.zSD, color = level)) +
      stat_ecdf(geom = "step") + ylab("ECDF")
   
}

## -----------------------------------------------------------------------------
if (require(ggplot2))
{
   ggplot(extendedres, aes(x = BMD.zSD)) + facet_wrap(~ level) +
      stat_ecdf(geom = "step") + ylab("ECDF")
   
}

## -----------------------------------------------------------------------------
bmdplot(extendedres, BMDtype = "zSD", 
                    facetby = "path_class", 
                    colorby = "level") + labs(col = "molecular level")

## -----------------------------------------------------------------------------
chosen_path_class <- c("Membrane transport", "Lipid metabolism", "Energy metabolism")
selectedres <- extendedres[extendedres$path_class %in% chosen_path_class, ]
bmdplot(selectedres, BMDtype = "zSD", add.CI = TRUE, 
                    facetby = "path_class", 
                    facetby2 = "level") + labs(col = "trend")

## -----------------------------------------------------------------------------
chosen_path_class <- c("Membrane transport", "Lipid metabolism", "Energy metabolism")
selectedres <- extendedres[extendedres$path_class %in% chosen_path_class, ]
bmdplotwithgradient(selectedres, BMDtype = "zSD",
               facetby = "path_class", facetby2 = "level")

## -----------------------------------------------------------------------------
sensitivityplot(extendedres, BMDtype = "zSD",
                group = "path_class", colorby = "level",
                BMDsummary = "first.quartile") 

## -----------------------------------------------------------------------------
# Preliminary optional alphabetic ordering of path_class groups
extendedres$path_class <- factor(extendedres$path_class, 
                levels = sort(levels(extendedres$path_class), decreasing = TRUE))
trendplot(extendedres, group = "path_class", facetby = "level") 

## -----------------------------------------------------------------------------
# Plot of the dose-response curves for a specific metabolic pathway
# in this example the "lipid metabolism" pathclass
LMres <- extendedres[extendedres$path_class == "Lipid metabolism", ]
curvesplot(LMres, facetby = "level", free.y.scales = TRUE, npoints = 100, line.size = 1,
           colorby = "trend",
           xmin = 0, xmax = 8) + labs(col = "DR_trend")

## -----------------------------------------------------------------------------
# Plot of the dose-response curves for a specific metabolic pathway
# in this example the "lipid metabolism" pathclass
curvesplot(selectedres, facetby = "path_class", facetby2 = "level",
           free.y.scales = TRUE, npoints = 100, line.size = 1,
           colorby = "trend",
           xmin = 0, xmax = 8) + labs(col = "DR_trend")

