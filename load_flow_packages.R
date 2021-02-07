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
packages <- c("colorspace", "lattice",
              "ggvis", "dygraphs",
              "DescTools", "viridis",
              "latticeExtra", "session")
sapply(X = packages, FUN = package_installer)
sapply(X = packages, FUN = require, character.only = TRUE)


## -----
## bioconductor
bioc_package_installer <- function(x){if(!requireNamespace(x))
                                          BiocManager::install(x, INSTALL_opts = '--no-lock')}
bioc_packages <-  c("flowCore", "flowViz",
                    "flowUtils", "flowStats",
                    "flowFP", "geneplotter",
                    "Rgraphviz", "ggcyto")
sapply(X = bioc_packages, FUN = bioc_package_installer)
sapply(X = bioc_packages, FUN = require, character.only = TRUE)


## -----
## required for merging flowsets into a single flowframe
source(file = "https://raw.githubusercontent.com/mac230/flow_scripts/master/set2frame.R")

