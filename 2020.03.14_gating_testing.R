## -----
## <<Required_Packages>> 
## check for Bioconductor and install if not available 
ifelse(!requireNamespace("BiocManager", quietly = TRUE),
       install.packages("BiocManager",
                        dependencies = TRUE,
                        repos = "http://cran.wustl.edu/",
                        quiet = TRUE),
       paste0("Bioconductor available"))
require("BiocManager")

## requireNamespace checks whether a package is available and loads if it is
## the return value is logical and the function throws an error if not available 
## if(!requireNamespace("DNAcopy")) paste0("package not available")
## check that the output of requireNamespace is truly logical:
## requireNamespace("dygraphs") == requireNamespace("lattice")     ## TRUE
## requireNamespace("dygraphs") == requireNamespace("fakepackage") ## FALSE
## ifelse(!requireNamespace("fakepackage"),
##        paste0("no such package"),
##        paste0("there is a package"))


## -----
## load packages or install if not available
## have to split these out by bioconductor vs. non-bioconductor
## non-bioconductor
package_installer <- function(x){
    if(!requireNamespace(x, quietly = TRUE))
        install.packages(x, dependencies = TRUE,
                         repos = "http://cran.wustl.edu/",
                         quiet = TRUE, INSTALL_opts = '--no-lock')}
packages <- c("colorspace", "lattice", "ggvis", "dygraphs")
sapply(X = packages, FUN = package_installer)
sapply(X = packages, FUN = require, character.only = TRUE)


## -----
## bioconductor
bioc_package_installer <- function(x){if(!requireNamespace(x))
                                          BiocManager::install(x, INSTALL_opts = '--no-lock')}
bioc_packages <-  c("flowCore", "flowViz", "flowUtils", "flowStats", "flowFP", "geneplotter", "ggcyto")
sapply(X = bioc_packages, FUN = bioc_package_installer)
sapply(X = bioc_packages, FUN = require, character.only = TRUE)


## -----
## required for merging flowsets into a single flowframe 
source(file = "https://raw.githubusercontent.com/mac230/flow_scripts/master/set2frame.R")


##-----
## <<Reading_FCS_Files>>
## user-specified options - these will change for each analysis depending on strains/reporters
##############
## USER INPUT:
##############
## no trailing '/' at the end!
base.dir       <- "~/data/flow/2020.03.14_new_gate_testing"
setwd(base.dir)
needed.dirs <- c("/fcs", "/results", "/tables")
dir.maker <- function(x){if(!dir.exists(paths = paste0(base.dir, x)))
                             dir.create(path = paste0(base.dir, x))}
sapply(X = needed.dirs, FUN = dir.maker)
work.dir       <- paste0(base.dir, "/fcs")
results.dir    <- paste0(base.dir, "/results")
tables.dir     <- paste0(base.dir, "/tables")


##-----
## [x]
## now set regex for getting flowsets of the different strains
## generally, should name fcs files as follows:
## strain    - by, rm, rpn4, rpn10
## reporter  - PSV, TFT, untagged
## replicate - 001, 002, etc... per strain
##############
## USER INPUT:
##############
no.reporter   <- ".*untagged.*fcs"
by.odc   <- "BY.*ODC.*.fcs"
rm.arg.tft   <- "RM.*Arg.*.fcs"
rpn4.thr.tft <- "rpn10.*Thr.*.fcs"


##############
## USER INPUT:
##############
## for later use in plots
strain.names <- c("no reporter", "BY ODC TFT", "RM Arg TFT", "rpn10 Thr TFT")


##-----
## [x]
## bind all regex to a list and use the list to read files
## the result here is a list of ungated flowSets
## each flowset has 'n' tubes (flowframes), where n is the number of replicates
## access a single flowFrame/tube w/ .e.g. "all.set[[1]][[1]]", which would be strain 1, tube 1
##############
## USER INPUT:
##############
setwd(work.dir)
all.strains <- list(no.reporter,
                    by.odc,
                    rm.arg.tft, 
                    rpn4.thr.tft)


all.set     <- lapply(all.strains, function(x){read.flowSet(files = NULL, path = ".", pattern = x, alter.names = T, min.limit = 1)})
sapply(X = 1:length(all.set), FUN = function(x){all.set[[x]]@phenoData@data$name})


##################
## END USER INPUT:
##################


##-----
## [x]
## write strain/replicate groupings to a table for inspection
setwd(tables.dir)
cat("File, Group, Strain", "\n", file = "strain_replicate_groupings.txt", append = F)
strain.group    <- as.list(seq(from = 1, to = length(all.set), by = 1))
replicates.out  <- unlist(lapply(1:length(all.set),
                                 function(x)
                                 {paste0(all.set[[x]]@phenoData@data$name, ", ",
                                         strain.group[[x]], ", ", strain.names[[x]])}))
replicate.table <- function(x){cat(c(x, "\n"), file = "strain_replicate_groupings.txt", append = T, sep = ", ")}
sapply(X = replicates.out, FUN = replicate.table)


##-----
## <<Transformation_and_Gating>>
## use the transform function to get the TFT/PSV parameters we want
## start by converting 0's in fluors to 1's via truncate transform
trunc.trans   <- truncateTransform("Convert 0's to 1's.", a = 1)
trunc.fluors  <- function(x){
    transform(x,
              `eGFP.A` = trunc.trans(`eGFP.A`),
              `mCherry.A` = trunc.trans(`mCherry.A`))}
all.set <- lapply(all.set, fsApply, trunc.fluors)

PSV.TFT.transform <- function(x){
    transform(x,
              `log_GFP` = log10(`eGFP.A`),
              `log_RFP` = log10(`mCherry.A`),
              `TFT_ratio` = log(`mCherry.A`/`eGFP.A`, base = 2),
              `PSV_ratio` = log(`eGFP.A`/`mCherry.A`, base = 2))}
all.set <- lapply(all.set, fsApply, PSV.TFT.transform)


##-----
## [x]
## get the total number of cells for each flowFrame
## nrow is passed as an optional arg to fsApply here
total.cells <- lapply(all.set, fsApply, nrow)

a <- all.set[[3]][[2]]
xyplot(`SSC.A` ~ `FSC.A`, data = a, smooth = F,
       filter = curv2Filter(x = "FSC.A", y = "SSC.A", bwFac = 3, gridsize = c(250,250)))
b <- split(x = a, f = new.c)
sapply(b, function(x){nrow(exprs(x))})

new.c <- curv2Filter(x = "FSC.A", y = "SSC.A", bwFac = 3, gridsize = c(250,250))
a <- all.set[[1]][[1]]
b <- split(x = a, f = new.c)
nrow(exprs(a[[1]]))
nrow(exprs(b[[1]]))
nrow(b[[4]])
## this is a nice 'initial' gate; combine this w/ simple FSC-based approach 
myf <- curv2Filter(x = "FSC.A", y = "SSC.A", bwFac = 3, gridsize = c(250,250))
## doesn't work b/c mode is always detector max value 
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

Mode(exprs(a[,1]))
a <- all.set[[3]][[1]]
densityplot(exprs(a[,1]))
a.d <- density(exprs(a[,1]))
str(exprs(a))

a.d

fam <- median(exprs(a$FSC.A))
fmy <- c(fam - (fam * 0.25), fam + (fam * 0.25))
sam <- median(exprs(a$SSC.A))
smy <- c(sam - (sam * 0.25), sam + (sam * 0.25))

rect <- rectangleGate(filterId = "FSC/SSC filter",
                      "SSC.A" = smy, "FSC.A" = fmy)
xyplot(`SSC.A` ~ `FSC.A`, data = a, smooth = F, filter = rect)
plot(exprs(a$FSC.A), exprs(a$SSC.A), pch = 19, cex = 0.02)
abline(v = fam)
abline(v = mean(exprs(a$FSC.A)), col = "green")
summary(a)

nf <- norm2Filter(x = "FSC.A", y = "SSC.A", method = "covMcd")

xyplot(`SSC.A` ~ `FSC.A`, data = a, smooth = F, filter = nf)

b <- all.set[[1]][[2]]

xyplot(`SSC.A` ~ `FSC.A`, data = b, smooth = F, filter = nf)

new_set <- read.flowSet(path = "~/Desktop/data/flow/2020.03.11_norm2Filter_testing/fcs")
a <- new_set[[10]]
nf <- norm2Filter(x = "FSC-A", y = "SSC-A", method = "cov.rob", scale.factor = 1)
xyplot(`SSC.A` ~ `FSC.A`, data = a, smooth = F, filter = nf, stat = TRUE, abs = TRUE, pos = 0.5)
ellipsoidGate

a_s <- split(new_set[[10]], km)
a_s


xyplot(`SSC.A` ~ `FSC.A`, data = a, smooth = F, filter = km)
xyplot(`SSC-A` ~ `FSC-A`, data = a_s[[1]], smooth = F)
xyplot(`SSC-A` ~ `FSC-A`, data = new_set[[10]], smooth = F)
a_s <- a_s[[1]]

kmeansFilter
