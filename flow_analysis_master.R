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
packages <- c("colorspace", "lattice", "ggvis", "dygraphs", "DescTools")
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
needed.dirs <- c("/fcs", "/results", "/tables", "/scripts")
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
by.thr.tft    <- "BY.*Thr.*.fcs"
rm.thr.tft    <- "RM.*Thr.*.fcs"
doa10.thr.tft <- "doa10.*Thr.*.fcs"


##############
## USER INPUT:
##############
## for later use in plots
strain.names <- c("no reporter", "BY Thr TFT", "RM Thr TFT", "doa10 Thr TFT")



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
                    by.thr.tft,
                    rm.thr.tft,
                    doa10.thr.tft)

all.set     <- lapply(all.strains, function(x){read.flowSet(files = NULL, path = ".", pattern = x, alter.names = T, min.limit = 1)})
## str(all.set[[1]]@phenoData@data$name)

## make sure we use the correct parameters 
rfp <- "mCherry.A"
stopifnot(
    length(
        grep(pattern = rfp,
             x = colnames(exprs(all.set[[1]][[1]])))) > 0)

##################
## END USER INPUT:
##################


## -----
## <<Color_Setup>>
## linking colors to strain names in R
## I think I should be able to make something
## akin to an lisp association list where
## there is a strain name and associated color
col.untagged <- c(color = gray(0.7),   name = "no reporter")
col.by       <- c(color = "#7A9BCCFF", name = ".*BY.*")
col.rm       <- c(color = "#CC7AAAFF", name = ".*RM.*")
col.rpn4     <- c(color = "#CCAB7AFF", name = ".*rpn4.*")
col.ubr1     <- c(color = "#88CCBBFF", name = ".*ubr1.*")
col.doa10    <- c(color = "#A3CC7AFF", name = ".*doa10.*")
cols.list    <- list(col.untagged, col.by, col.rm, col.rpn4, col.ubr1, col.doa10)

col.out <- sapply(X = cols.list, FUN = function(x){
                      grepl(pattern = x["name"], x = strain.names )
                  })
col.out <- as.logical(unlist(sapply(1:ncol(col.out), FUN = function(x){
                      max(col.out[, x])
                      })))
all.cols <- unlist(sapply(X = cols.list[col.out], FUN = function(x){identity(x["color"])}))


## output a dummmy plot to assess strain/color mapping
setwd(results.dir)
pdf(file = "color_mapping.pdf", height = 7, width = 7, bg = "transparent")
barplot(rep(4, length(strain.names)), col = all.cols, ylim = c(0, 5.5))
box()
legend(x = "topleft", legend = strain.names, lty = 1, lwd = 7.5, col = all.cols, bg = "white")
legend(x = "topright", y = NA,
       legend = unlist(lapply(X = cols.list, FUN = function(x){identity(x)["name"]})),
       col = unlist(lapply(X = cols.list, FUN = function(x){identity(x)["color"]})),
       lty = 1, lwd = 7.5,  bg = "white")
dev.off()


##-----
## [x]
## write strain/replicate groupings to a table for inspection
## view w/ 'column -t -s "," ./tables/strain_replicate_groupings.txt'
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
## <<TFT_Transformation>>
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
              `PSV_ratio` = log(`eGFP.A`/`mCherry.A`, base = 2),
              ## 'no log' TFT ratio
              `nl_TFT_ratio` = (`mCherry.A`/`eGFP.A`)
              )}
all.set <- lapply(all.set, fsApply, PSV.TFT.transform)


##-----
## [x]
## get the total number of cells for each flowFrame
## nrow is passed as an optional arg to fsApply here
total.cells <- lapply(all.set, fsApply, nrow)


##-----
## <<Cell_Gating>>
## [x]
## 02.27.2019 try this w/ curv2Filter w/ a big bandwidth setting to grab the
## main cloud of cells we take only cells in 'area 1' (the gate), not 'rest'
## (the cells outside the gate)
## <<FSC_SSC_Gate>>
initial.split <- function(x) {
    split(x, f = curv2Filter(x = "FSC.A", y = "SSC.A", bwFac = 7, gridsize = c(250,250)),
          population = "area 1", flowSet = TRUE, codeflowSet = TRUE)
}
## this object is a list of flowsets of the cells in the initial gate (area 1)
## each flowset in this list contains only 1 flowFrame
initial.split.all <- lapply(all.set, fsApply, initial.split)


##-----
## [x]
## plot the cells w/ their associated filter gate
setwd(results.dir)
dir.create(path = paste0(results.dir, "/cell_gate_plots"))
cell.gate.dir <- paste0(results.dir, "/cell_gate_plots")
setwd(cell.gate.dir)

xy.initial.pdf <- function(x){
    pdf(file = paste0("initial_", x@description$"TUBE NAME", ".pdf"), height = 7, width = 7)
    print(xyplot(`SSC.A` ~ `FSC.A`, data = x,
                 filter = curv2Filter(x = "FSC.A", y = "SSC.A", bwFac = 10, gridsize = c(250,250)),
                 smooth = F))
    dev.off()
}
lapply(all.set, fsApply, xy.initial.pdf)


##-----
## [x] - 2020.04.19 - no longer using due to fsc gating approach below
## plot the results of the pre-filter plus curv2Filter gating
## start by undoing the complicated list structure the filter operation creates
## this yields a list of flowSets
## initial.curv.split <- unlist(initial.split.all)
##setwd(cell.gate.dir)
##xy.initial.curv.pdf <- function(x) {
##    pdf(file = paste0("curv_", x@description$"TUBE NAME", "_.pdf"), height = 7, width = 7)
##    print(xyplot(`SSC.A` ~ `FSC.A`, data = x,
##                 filter = curv2Filter(x = "FSC.A", y = "SSC.A", bwFac = 2, gridsize = c(250,250)),
##                 smooth = F))
##    dev.off()
##}
##lapply(initial.curv.split, fsApply, xy.initial.curv.pdf)


## -----
## <<FSC_Gate>>
## a function to gate the cells to include only haploids.
## we identify these as a sharp peak in the lower end of
## the fsc density plot.  I take 10% above and below the
## max density value
fsc.gate.generator <- function(x){
    fsc.dens  <- density(exprs(x[, 1]))
    ## return the index of the maximum y value of the density estimate
    fsc.max   <- fsc.dens[[1]][which.max(fsc.dens[[2]])]
    fsc.upper <- (fsc.max * 0.10) + fsc.max
    fsc.lower <- fsc.max - (fsc.max * 0.10)
    fsc.gate  <- c(fsc.lower, fsc.upper)
}

curv.split <- function(x){
    split(x, f = rectangleGate("FSC.A" = fsc.gate.generator(x)),
          population = "defaultRectangleGate+",
          flowSet = T, codeflowSet = T)}
curv.set <- lapply(all.set, fsApply, curv.split)


##-----
## [x]
## plot the results of the pre-filter plus curv2Filter gating
## start by undoing the complicated list structure the filter operation creates
## this yields a list of flowSets
setwd(cell.gate.dir)
xy.fsc.curv.pdf <- function(x){
    pdf(file = paste0("curv_", x@description$"TUBE NAME", "_.pdf"), height = 7, width = 7)
    print(xyplot(`SSC.A` ~ `FSC.A`, data = x, main = x@description$"TUBE NAME",
                 filter = rectangleGate("FSC.A" = fsc.gate.generator(x)),
                 smooth = F))
    dev.off()
    }
lapply(all.set, fsApply, xy.fsc.curv.pdf)

## plot the fsc density and associated gate
## I use 'Map' here to color the plots by strain 
fsc.density.curv.pdf <- function(data, color){
    fsApply(data, function(x){
                pdf(file = paste0("fsc_density_", x@description$"TUBE NAME", "_.pdf"),
                    height = 7, width = 7)
                plot(density(exprs(x[, 1])),
                     xlab = colnames(exprs(x))[1],
                     main = x@description$"TUBE NAME",
                     col = color, lwd = 2.5)
                abline(v = fsc.gate.generator(x),
                       col = gray(0.4),
                       lty = 3, lwd = 2.5)
                dev.off()
            }
            )}

Map(f = fsc.density.curv.pdf, all.set, all.cols)


gated.xy.fsc.curv.pdf <- function(x){
    pdf(file = paste0("xy_sub_population_curv_", x@description$"TUBE NAME", "_.pdf"), height = 7, width = 7)
    print(xyplot(`SSC.A` ~ `FSC.A`, data = x, main = x@description$"TUBE NAME",
                 filter = rectangleGate("FSC.A" = fsc.gate.generator(x)),
                 smooth = F))
    dev.off()
}
lapply(curv.set, function(x){lapply(x, fsApply, gated.xy.fsc.curv.pdf)})


## -----
## <<Gate-Based_Subsetting>>
## [x]
## the output of the filtering operation is a list
## curv.set:
## 1. curv.set -> list of flowsets
##    curv.set[[1]] -> flowSet w/ 4 experiments

##    2. curv.set[[1]][[1]] -> flowSet
##       flowSet of the different curv gates

##       3. curv.set[[1]][[1]][[1]] -> the actual flowframe of each gate
##                                     ("defaultRectangleGate+")

## use unlist to get a simpler list structure
## the result is a list of flowsets


## now convert the list of fsc-gated flowsets into a list of list of
## Start by unlisting the original dataset.
## Then loop over the unlisted data and grab the flowframe of interest.
curv.set <- unlist(curv.set)
c.final.frame  <- list()

for(i in 1:length(curv.set)){
    c.final.frame[[i]] <- curv.set[[i]][[1]]
}


##-----
## [x]
## list of flowsets that result from the first filter
## i.final.frame is a list of flowFrames
i.set <- unlist(initial.split.all)
i.final.frame <- list()
for(i in 1:length(i.set)){
i.final.frame[[i]] <- i.set[[i]][[1]]
}


##-----
## [x]
## list of flowsets of the ungated cells
## u.final.frame is a list of flowFrames
u.set <- unlist(lapply(all.set, fsApply, list))
u.final.frame <- list()
for(i in 1:length(u.set)){
u.final.frame[[i]] <- u.set[[i]]
}


##-----
## [x]
## get the names of the original samples and use them as regex
## 'grepl' returns logical output; w/ 'lapply', test each regex on the list of file/tube names
tube.names <- unlist(lapply(c.final.frame, function(x) {print(x@description$GUID.original)}))
strain.regex.logical <- lapply(all.strains, function(x){grepl(x = tube.names, pattern = x)})


##-----
## [x]
## now use the logical vectors for grouping
all.data <- list(curv = c.final.frame, initial = i.final.frame, ungated = u.final.frame)
names(all.data)


## 1. all.data is a list of 3 lists, each of these 3 lists is a list of 24 ungrouped flowframes
## 2. lapply all.data to pass a list of 24 flowframes to a function
## 3. make that function an lapply to pass each of 24 flowframes to some test
## 4. end result should be a list of 3, with each of these 3 lists having 5 lists (strains)
## 3 levels to all.data -> all.data[[I-gated_set]][[II-strain_replicates]][[III-individual flowframe]]
all.data <- lapply(all.data,
                   function(x){
                       lapply(seq_along(strain.regex.logical),
                              function(y){
                                  x[strain.regex.logical[[y]]]
                              })
                   })


##-----
## [x]
## now we do logicle transform for plotting 
logicle.trans <- logicleTransform(transformationId = "logicle fluor transform")
logicle.func <- function(x){
    transform(x,
              `log_GFP` = logicle.trans(`eGFP.A`),
              `log_RFP` = logicle.trans(`mCherry.A`))}

all.data.l <- lapply(1:length(all.data), function(x){
                         lapply(all.data[[x]], function(y){
                                    lapply(y, function(q){
                                               logicle.func(q)}
                                           )}
                                )})


##-----
## [x]
## logicle plots of individual samples
setwd(results.dir)
dir.create(path = paste0(results.dir, "/logicle_2D_plots_individual"))
logicle.dir <- paste0(results.dir, "/logicle_2D_plots_individual")
setwd(logicle.dir)

lapply(1:length(all.data.l), function(x){
           lapply(all.data.l[[x]], function(y){
                      lapply(y, function(q){
                                 pdf(file = paste0(names(all.data[x]),
                                                   "_",
                                                   q@description$'TUBE NAME',
                                                   ".pdf"),
                                     height = 7,
                                     width = 7)
                                 print(xyplot(`log_RFP` ~ `log_GFP`,
                                              data = q,
                                              smooth = F,
                                              strip = paste0(names(all.data[x]),
                                                             " ",
                                                             gsub("_",
                                                                  " ",
                                                                  q@description$'TUBE NAME')),
                                              prepanel=function(){return(list(xlim = c(0, 4), ylim = c(0, 4)))}))
                      dev.off()
                      })
                  })
           })


##-----
## [x]
## per gate per group logicle plots
## set up lists for each gate
setwd(results.dir)
dir.create(path = paste0(results.dir, "/logicle_2D_plots_groups"))
logicle.groups.dir <- paste0(results.dir, "/logicle_2D_plots_groups")
setwd(logicle.groups.dir)


## set up lists for each gate
c.frame.l   <- vector(mode = "list", length = length(strain.names))
i.frame.l   <- vector(mode = "list", length = length(strain.names))
u.frame.l   <- vector(mode = "list", length = length(strain.names))
all.frame.l <- list(curv = c.frame.l, initial = i.frame.l, ungated = u.frame.l)


## convert individual flowsets to single flowframes
all.frame.l <- lapply(1:length(all.frame.l), function(x){
                          lapply(1:length(all.frame.l[[x]]), function(y){
                                     all.frame.l[[x]][[y]] <- set2Frame(as(unlist(all.data.l[[x]][[y]]), "flowSet"))
                                 })
                      })


## name the frames that will comprise each gate
## 'gsub' removes spaces 
names(all.frame.l) <- names(all.data)
for(i in 1:length(all.frame.l)){
    names(all.frame.l[[i]]) <- gsub(pattern = " ", replacement = "_", x = strain.names)
}


## plot
## 'gsub' here ensures no spaces in file names
lapply(1:length(all.frame.l), function(x){
           lapply(1:length(all.frame.l[[x]]), function(y){
                      pdf(file = paste0(names(all.frame.l[x]),
                                        "_",
                                        gsub(" ", "_",
                                             names(all.frame.l[[x]][y])), ".pdf"),
                          height = 7,
                          width = 7)
                      print(xyplot(`log_RFP` ~ `log_GFP`,
                                   data = all.frame.l[[x]][[y]],
                                   smooth = F,
                                   strip = paste0(names(all.data[x]),
                                                  " ",
                                                  gsub("_", " ", names(all.frame.l[[x]][y]))),
                                   prepanel=function(){
                                       return(list(xlim = c(0, 4),
                                                   ylim = c(0, 4)))
                                   }))
                      dev.off()
                  })
       })


##-----
## <<Cell_Count_Table_and_Plots>>
## [x]
## check that final data result has expected number of cells
## write the number of cells from each gating step to a table
setwd(tables.dir)

c.counts <- unlist(lapply(all.data[[1]], function(x){
                              lapply(x, function(y){
                                         nrow(y)
                                     })
                          }))

i.counts <- unlist(lapply(all.data[[2]], function(x){
                              lapply(x, function(y){
                                         nrow(y)
                                     })
                          }))

u.counts <- unlist(lapply(all.data[[3]], function(x){
                              lapply(x, function(y){
                                         nrow(y)
                                     })
                          }))

table.names <- unlist(lapply(all.data[[1]], function(x){
                                 lapply(x, function(y){
                                            print(y@description$GUID.original)
                                        })
                             }))

table_factor <- vector(mode = "list", length = length(all.set))
for(i in 1:length(all.set)){
    table_factor[[i]] <- rep(i, times = length(all.set[[i]]))
    }
table_factor <- unlist(table_factor)

all.cell.counts <- list(curv_gate = c.counts, initial_gate = i.counts, ungated = u.counts, table.factor = table_factor, names = table.names)
write.table(x = all.cell.counts, file = "cell_counts_by_gate.txt", append = F, sep = ",", quote = F, row.names = F)

cell.data <- read.csv(file = "cell_counts_by_gate.txt", header = T)
cell.data$table.factor <- factor(x = cell.data$table.factor, levels = unique(cell.data$table.factor), labels = strain.names)


##-----
## [x]
## box plots of cell counts x group
setwd(results.dir)
dir.create(path = paste0(results.dir, "/cell_count_plots"))
cell.count.dir <- paste0(results.dir, "/cell_count_plots")
setwd(cell.count.dir)
lapply(1:length(all.data), function(x){
           pdf(file = paste0("cell_count_boxplot_", names(all.data[x]), ".pdf"), height = 7, width = 7, bg = "transparent")
           par(cex.axis = 0.8)
           boxplot(cell.data[, x] ~ cell.data$table.factor, ylab = "Cell Count", col = gray(0.9))
           dev.off()})


## single pdf of all cell counts
par(mfrow = c(3,1))
pdf(file = "all_cell_counts_boxplot_.pdf", height = 7, width = 7, bg = "transparent")
lapply(1:length(all.data), function(x){
           par(cex.axis = 0.8)
           boxplot(cell.data[, x] ~ cell.data$table.factor, ylab = "Cell Count", col = gray(0.9))
       })
dev.off()
par(mfrow = c(1,1))


##-----
## <<Check_Gate_Groupings>>
## [x]
## get tube names for each set of samples
## unlist step results in a list of lists -> out[[I - gate]][[II - strain]]
out <- lapply(all.data, function(x){
                  lapply(x, function(y){
                             unlist(lapply(y, function(q){
                                        print(strsplit(x = q@description$GUID.original, split = "_[0-9]{3}.fcs"))
                                    }))
                         })
              })


##-----
## [x]
## write output to a file that gives replicate grouping across frames of 'out'
setwd(tables.dir)
for(i in 1:length(out)){
    for(j in 1:length(i)){
        for(k in 1:length(j)){
            cat(paste0(out[[i]], "_", names(out[i])), file = "gates_by_groups_list.txt", append = T, sep = "\n")
        }
    }
}


##-----
## <<Data_Frame_Conversion>>
## [x]
## convert to data frames
all.data.e <- lapply(all.data, function(x){
                         lapply(x, function(y){
                                    lapply(y, function(q){
                                               q <- as.data.frame(exprs(q))
                                           })
                                })
                     })


##-----
## [x]
## save the data as an R object for later loading if needed
setwd(base.dir)
dir.create(path = paste0(base.dir, "/R_objects"))
objs.dir <- paste0(base.dir, "/R_objects")
setwd(objs.dir)
## flowsets
saveRDS(all.data, file = "curv_initial_ungated_flowsets")
## exprs matrix
saveRDS(all.data.e, file = "curv_initial_ungated_flowsets_exprs")
## test <- readRDS(file = "curv_initial_ungated_sets")


##-----
## <<Between_Groups_Processing>>
## [x]
## between groups stuff
## start by grouping replicates in a list of lists
all.groups.e <- lapply(all.data.e, function(x){
                           lapply(x, function(y){
                                      x <- do.call("rbind", y)
                                  })
                       })

## name the backgrounds in each dataset
## remember, have to 'for' loop this; can't use 'lapply'
for(i in seq_along(all.groups.e)){
           names(all.groups.e[[i]]) <- strain.names
       }

## write out the total number of cells for each background
setwd(tables.dir)
cell.counts <- vector(mode = "list", length = length(all.groups.e))
names(cell.counts) <- names(all.groups.e)

cell.counts <- lapply(all.groups.e, function(x){
                          unlist(lapply(1:length(x), function(y){
                                            nrow(x[[y]])
                                        }))
                      })

## assign names
for(i in 1:length(cell.counts)){
    names(cell.counts[[i]]) <- strain.names
}

## create a new table to dump output into 
cat(paste0("Cell_Count, Strain, Gate,", "\n"), file = "total_cells_per_group_gate.txt")

## loop over the data to write the table 
for(i in seq_along(cell.counts)){
    for(h in seq_along(cell.counts[[i]])){
        cat(paste0(cell.counts[[i]][[h]], ", ", 
                   gsub(" ", "_", names(cell.counts[[i]][h])), ", ",
                   names(cell.counts[i])),
                   file = "total_cells_per_group_gate.txt",
                   append = T,
                   sep = "\n")
        }
    }


##-----
## <<Between_Groups_Plots>>
## [x]
## need to add names to individual elements
setwd(results.dir)
dir.create(path = paste0(results.dir, "/between_groups_plots"))
between.groups.dir <- paste0(results.dir, "/between_groups_plots")
setwd(between.groups.dir)

## set up names and limits for parameters - these may have to change
## these are for density plots, so x axis is what's below
## 1 - "FSC.A", 2 - "SSC.A", 3 - "eGFP.A", 4 - "mCherry.A", 5 - "Time"
## 6 - "log_GFP", 7 - "log_RFP", 8 - "TFT_ratio", 9 - "PSV_ratio", 10 - "nl_TFT_ratio"
x.lab   <- gsub(pattern = "_", replacement = " ", names(all.groups.e[[1]][[1]]))
x.min   <- c(0, 0, 0, 0, 0, 2, 2, -6, -2, 0)
x.max   <- c(2.5e5, 2e5, 2e4, 2e4, 2.5e3, 5, 5, 2, 6, 1)
leg.pos <- c(rep("topright", 5), "topleft", rep("topright", 3))

## this can help catch errors related to the correct number of parameters above 
stopifnot(length(exprs(all.set[[1]][[1]][1, ])) == length(x.min))

## ylim takes the min and max of the density each parameter
## i.e., the values will always be low
y.lim <- vector(mode = "list", length = 3)
for(h in 1:length(y.lim)){
    for(i in 1:length(x.lab)){
        y.lim[[h]][i] <- max(
            unlist(lapply(all.groups.e[[h]],
                          function(x){1.1 * max(density(x[, i])$y)
                                           })
                   ))}}


Map(f = function(gate, name){
        for(i in 1:length(x.lab)){
            ## get y limit for density plot
            lapply(1:length(all.groups.e[[1]]),
                   function(x){
                       ## plot
                       pdf(file = paste0(name, "_",
                                         gsub(" ", "_", x.lab[i]),
                                         "_between_groups", ".pdf"),
                           height = 7, width = 7, bg = "transparent")
                       plot(density(all.groups.e[[gate]][[x]][, i], adjust = 0.75),
                            col = "white",
                            ylim = c(0, y.lim[[gate]][i]),
                            xlim = c(x.min[i], x.max[i]),
                            xlab = x.lab[i], main = name)

                       ## map across gates and do so on a per-strain basis
                       Map(f = function(x, y){
                               lines(density(all.groups.e[[gate]][[x]][, i],
                                             adjust = 0.75),
                                     col = y, lwd = 2)},
                           x = 1:length(all.groups.e[[gate]]), y = all.cols)

                       ## legend
                       legend(x = "topleft", legend = strain.names,
                              lty = 1, lwd = 5,
                              col = all.cols, bg = "white")
                       dev.off()

                   })
        }}, gate = seq_along(all.groups.e), name = names(all.groups.e))


##-----
## <<Between_Groups_Replicates>>
## [x]
## between groups w/ replicates
setwd(results.dir)
dir.create(path = paste0(results.dir, "/between_groups_w_replicates_plots"))
between.groups.replicates.dir <- paste0(results.dir, "/between_groups_w_replicates_plots")
setwd(between.groups.replicates.dir)

Map(f = function(gate, name){
        for(i in 1:length(x.lab)){

            ## get y limit for density plot
            lapply(1:length(all.groups.e[[1]]),
                   function(x){
                       ## plot
                       pdf(file = paste0(name, "_",
                                         gsub(" ", "_", x.lab[i]),
                                         "_between_groups_w_replicates.pdf"),
                           height = 7,
                           width = 7,
                           bg = "transparent")
                       
                       plot(density(all.groups.e[[gate]][[x]][, i], adjust = 0.75),
                            col = "white",
                            ylim = c(0, y.lim[[gate]][i]),
                            xlim = c(x.min[i], x.max[i]),
                            xlab = x.lab[i],
                            main = name)

                       ## map across gates and do so on a per-strain basis
                       Map(f = function(x, y){
                        lines(density(all.groups.e[[gate]][[x]][, i],
                                      adjust = 0.75), col = y, lwd = 3)
                            lapply(1:length(all.data.e[[gate]][[x]]),
                                   function(n){
                                       lines(density(all.data.e[[gate]][[x]][[n]][, i],
                                                     adjust = 0.75),
                                             col = y,
                                             lwd = 0.8)})},
                            x = seq_along(all.groups.e[[gate]]), y = all.cols)

                       ## legend
                       legend(x = "topleft",
                              legend = strain.names,
                              lty = 1, lwd = 5,
                              col = all.cols, bg = "transparent")
                       dev.off()
                   }
        )}}, gate = 1:length(all.groups.e), name = names(all.groups.e))


## -----
## <<Between_Gates_Plots>>
## [x]
## between gates plots
setwd(results.dir)
dir.create(path = paste0(results.dir, "/between_gates_plots"))
between.gates.dir <- paste0(results.dir, "/between_gates_plots")
setwd(between.gates.dir)

Map(f = function(strain, name){
        for(i in 1:length(x.lab)){

            ## get y limit for density plot
            y.lim <- c(0, max(unlist(lapply(1:3,function(x)
                                            {max(density(all.groups.e[[x]][[strain]][, i])$y)
                                            }))))

            ## plot 
            pdf(file = paste0(gsub(" ", "_", x.lab[i]), "_",
                              gsub(" ", "_", name),
                              "_between_gates.pdf"),
                height = 7,
                width = 7,
                bg = "transparent")
            plot(density(all.groups.e[[1]][[strain]][, i],
                         adjust = 0.75),
                 col = gray(0.8),
                 ylim = y.lim,
                 xlim = c(x.min[i], x.max[i]),
                 xlab = x.lab[i],
                 main = name)

            ## map across gates and do so on a per-strain basis
            Map(f = function(x, y){
                    lines(density(all.groups.e[[x]][[strain]][, i],
                                  adjust = 0.75),
                          col = all.cols[strain],
                          lty = y,
                          lwd = 2)},
                x = 1:3, y = 1:3)

            ## legend
            legend(x = "topleft",
                   legend = c("curv", "initial", "ungated"),
                   lty = 1:3,
                   lwd = 5,
                   col = all.cols[strain],
                   bg = "white",
                   seg.len = 3)
            dev.off()

        }
    }, strain = seq_along(strain.names), name = strain.names)


## <<Replicate_Plots>>
##-----
## [x]
## replicate plots
## this gives a nice layout w/ replicates separately and grouped 
setwd(results.dir)
dir.create(path = paste0(results.dir, "/replicate_plots"))
replicates.dir <- paste0(results.dir, "/replicate_plots")
setwd(replicates.dir)

lapply(1:length(all.groups.e), function(x){
           lapply(1:length(all.groups.e[[x]]), function(y){
                      for(i in 1:length(x.lab)){
                          pdf(file = paste0(gsub(" ", "_", names(all.groups.e[x])),
                                            "_",
                                            gsub(" ", "_", x.lab[i]),
                                            "_",
                                            gsub(" ", "_", strain.names[y]),
                                            "_replicate_plot.pdf"),
                              height = 7, width = 7,
                              bg = "transparent")

                          y.i.lim <- c(0, 1.2 * max(density(all.groups.e[[x]][[y]][, i])$y))
                              plot(density(all.groups.e[[x]][[y]][, i], adjust = 0.75),
                                   xlab = x.lab[i],
                                   main = strain.names[y],
                                   xlim = c(x.min[i], x.max[i]),
                                   ylim = y.i.lim,
                                   col = "white")

                              cols <- rainbow(n = length(all.data.e[[x]][[y]]),
                                              s = 0.7,
                                              v = 0.7,
                                              alpha = 1,
                                              end = 0.8)

                              ## the individual replicates line
                          lapply(seq_along(all.data.e[[x]][[y]]),
                                 function(r){lines(density(all.data.e[[x]][[y]][[r]][, i],
                                                           adjust = 0.75),
                                                   col = cols[r])})
                              
                              names <- fsApply(all.set[[y]], function(x){
                                                   paste0(x@description$"TUBE NAME")})

                              ## the overall group line
                              lines(density(all.groups.e[[x]][[y]][, i], adjust = 0.75),
                                    col = gray(0.7), lwd = 2)
                              
                              legend(x = "topleft",
                                     legend = c(names, "all replicates"),
                                     lty = 1,
                                     lwd = 5,
                                     col = c(cols, gray(0.7)),
                                     bty = "",
                                     bg = "transparent")

                          dev.off()
                      }
                  })
                  }
           )


##-----
## <<2D_Scatter_Plots>>
## [x]
setwd(results.dir)
dir.create(path = paste0(results.dir, "/2D_scatter_plots"))
scatterplot.dir <- paste0(results.dir, "/2D_scatter_plots")
setwd(scatterplot.dir)

lapply(1:length(all.groups.e), function(x){
           lapply(1:length(all.groups.e[[x]]), function(y){

                      pdf(file = paste0(names(all.groups.e[x]),
                                        "_",
                                        strain.names[y],
                                        "_GFP_mCh_scatter.pdf"),
                          height = 7, width = 7, bg = "transparent")

                      plot(all.groups.e[[x]][[y]]$log_GFP,
                           all.groups.e[[x]][[y]]$log_RFP,
                           xlim = c(1, 5),
                           ylim = c(1, 5),
                           pch = 19,
                           cex = 0.1,
                           col = all.cols[y],
                           xlab = "log10 sfGFP", ylab = "log 10 mCherry",
                           main = strain.names[y])
                      dev.off()
                  })
       })


##-----
## <<Between_Groups_Boxplots>>
## [x]
## between groups boxplots 
## merge exprs datasets into a single dataframe,
## use nrow to create a factor for each strain

## create a factor for each gate corresponding to the individual strains
## e.g. nrow(all.groups.e[[1]][[1]]) = number of cells of one strain type
strain.factor <- vector(mode = "list", length = length(all.groups.e))
strain.factor <- lapply(all.groups.e, function(x){
                            unlist(lapply(1:length(x), function(y){
                                              rep(x = y, times = nrow(x[[y]]))
                                          }))
                                   })

## now use the 'strain.names' object to label the new factor 
strain.factor <- lapply(1:length(strain.factor), function(x){
                            factor(x = strain.factor[[x]],
                                   levels = unique(strain.factor[[x]]),
                                   labels = strain.names)
                            })

## so, here we merge everything into a single dataset
## at this point, 'strain.factor' is still separate from
## our actual data, so merge it in at each gate level
## the new object 'all.groups.df' remains a list w/ 3
## levels, corresponding to our 3 levels of gating 
all.groups.df <- vector(mode = "list", length = length(all.groups.e))
all.groups.df <- lapply(all.groups.e, function(x){
                            do.call("rbind", x)
})

## add the strain factor to the dataframes
for(i in 1:length(all.groups.df)){
    all.groups.df[[i]]$strain <- strain.factor[[i]]
}

## we'll loop over columns for plotting, so make a list to exclude the 'strain' factor column
fac.test <- lapply(all.groups.df, function(x){
                       unlist(lapply(1:ncol(all.groups.df[[1]]), function(y){
                              !is.factor(x[, y])
                              }))
                       })

fac.test <- lapply(fac.test, function(x){
                       x <- x[x > 0]
                   })

## create pdfs
setwd(results.dir)
dir.create(path = paste0(results.dir, "/between_groups_boxplots"))
boxplots.dir <- paste0(results.dir, "/between_groups_boxplots")
setwd(boxplots.dir)

invisible(lapply(1:length(all.groups.df), function(x){
                     lapply(1:length(fac.test[[x]]), function(y){
                                pdf(file = paste0(names(all.groups.df[x]),
                                                  "_",
                                                  names(all.groups.df[[x]][y]),
                                                  "_strain_boxplot.pdf"),
                                    height = 7,
                                    width = 7,
                                    bg = "transparent")
                                par(cex.axis = 0.8) 
                                boxplot(all.groups.df[[x]][, y] ~ all.groups.df[[x]]$strain,
                                        main = names(all.groups.df[x]),
                                        names = strain.names,
                                        ylab = gsub("_", " ", names(all.groups.df[[x]][y])),
                                        col = gray(0.95),
                                        pars = list(outpch = 19,
                                                    outcol = "#112255AA",
                                                    outcex = 0.4))
                                dev.off()
                            })
                 }))


##-----
## <<Replicate_Mean_Stripcharts>>
## [x]
## start by creating a list of lists identical in size
## to the individual replicates list
all.data.means <- vector(mode = "list", length = length(all.data.e))

## now get the mean of each replicate, typically a sample > 10k cells
## the result of this operation is a list of list that's the same
## structure as 'all.data.e', but instead of, e.g., 10,000 rows, each
## individual sample now has only a single value (the mean) for each
## parameter 
all.data.means <- lapply(1:length(all.data.e), function(x){
                             all.data.means[[x]] <- lapply(1:length(all.data.e[[x]]), function(y){
lapply(1:length(all.data.e[[x]][[y]]), function(i){
           sapply(X = 1:ncol(all.data.e[[x]][[y]][[i]]),
                  FUN = function(m){mean(all.data.e[[x]][[y]][[i]][, m])})
       })
})
})

## name the gates, strains, and replicates
names(all.data.means) <- names(all.data.e)
## i = gates
for(i in seq_along(all.data.means)){
    names(all.data.means[[i]]) <- names(all.groups.e[[i]])
    ## h = strains
    for(h in seq_along(all.data.means[[i]])){
        ## g = replicates        
        for(g in seq_along(all.data.means[[i]][[h]])){
            rep.names <- vector()
            rep.names[g] <- all.set[[h]][[g]]@description$"TUBE NAME"
            names(all.data.means[[i]][[h]]) <- rep.names
            names(all.data.means[[i]][[h]][[g]]) <- names(all.groups.e[[1]][[1]])
        }}}


## create a similar list that gets the mean of the replicate means
## you'll use this later for plotting and statistics
## in the first step, bind each replicate into a single list that represents
## a strain
all.data.mean.lines <- lapply(all.data.means, function(x){
                                  lapply(1:length(x), function(y){
                                             colMeans(do.call("rbind", x[[y]]))
                                         })
                              })

for(i in seq_along(all.data.mean.lines)){
    names(all.data.mean.lines[[i]]) <- strain.names
}


## now bind the list of strains together via 'rbind'
all.data.mean.lines <- lapply(all.data.mean.lines, function(x){
                                  as.data.frame(do.call("rbind", x))
                              })

all.data.means <- lapply(all.data.means, function(x){
                             lapply(x, function(y){
                                        as.data.frame(do.call("rbind", y))
                                    })
                         })

all.data.means <- lapply(all.data.means, function(x){
                             do.call("rbind", x)
                         })


strain.rep.factor <- vector(mode = "list",
                            length = length(all.data.e))
strain.rep.factor <- lapply(1:length(all.data.e), function(x){
                                strain.rep.factor[[x]] <- factor(x = unlist(lapply(1:length(all.data.e[[x]]), function(y){
rep(y, times = length(all.data.e[[x]][[y]]))
})), labels = strain.names)
})


for(i in 1:length(all.data.means)){
    all.data.means[[i]]$strain <- strain.rep.factor[[i]]
}


## now plot across gates/parameters
## set up some dirs first
setwd(results.dir)
dir.create(path = paste0(results.dir, "/between_groups_mean_stripcharts"))
mean.stripcharts.dir <- paste0(results.dir, "/between_groups_mean_stripcharts")
setwd(mean.stripcharts.dir)

lapply(1:length(all.data.means), function(x){
           lapply(1:(ncol(all.data.means[[x]]) - 1), function(y){
                      pdf(file = paste0(names(all.data.means[x]), "_",
                                        names(all.data.means[[x]][y]),
                                        "_mean_betw_groups_strip.pdf"),
                          height = 7,
                          width = 7,
                          bg = "transparent")
                      ## first plot is 'dummy' plot so I can add lines
                      ## for the means.  After the lines are added,
                      ## overplot the same data again so individual
                      ## points aren't hidden behind the mean lines
                      stripchart(all.data.means[[x]][, y] ~ all.data.means[[x]]$strain,
                                 vertical = T,
                                 pch = NA,
                                 cex.axis = 0.8,
                                 main = names(all.data.means[x]),
                                 ylab = gsub("_", " ", names(all.data.means[[x]][y])))
                      lapply(1:length(strain.names),
                             function(q){
                                 lines(x = c(q - 0.25, q + 0.25),
                                       y = rep(all.data.mean.lines[[x]][q, y], 2),
                                       lwd = 2.5,
                                       col = "black")})
                      stripchart(all.data.means[[x]][, y] ~ all.data.means[[x]]$strain,
                                 vertical = T,
                                 cex.axis = 0.8,
                                 pch = 21,
                                 lwd = 1.25,
                                 col = gray(0.2),
                                 bg = gray(0.8),
                                 add = T,
                                 method = "jitter",
                                 jitter = 0.1,
                                 cex = 1.25)
                      dev.off()
                             })
                  })


## -----
## <<Replicate_Median_Stripcharts>>
## [x]
## start by creating a list of lists identical in size
## to the individual replicates list
all.data.medians <- vector(mode = "list", length = length(all.data.e))

## now get the mean of each replicate, typically a sample > 10k cells
## the result of this operation is a list of list that's the same
## structure as 'all.data.e', but instead of, e.g., 10,000 rows, each
## individual sample now has only a single value (the mean) for each
## parameter 
all.data.medians <- lapply(1:length(all.data.e), function(x){
                             all.data.means[[x]] <- lapply(1:length(all.data.e[[x]]), function(y){
lapply(1:length(all.data.e[[x]][[y]]), function(i){
           sapply(X = 1:ncol(all.data.e[[x]][[y]][[i]]),
                  FUN = function(m){median(all.data.e[[x]][[y]][[i]][, m])})
       })
})
})

## name the gates, strains, and replicates
names(all.data.medians) <- names(all.data.e)
## i = gates
for(i in seq_along(all.data.medians)){
    names(all.data.medians[[i]]) <- names(all.groups.e[[i]])
    ## h = strains
    for(h in seq_along(all.data.medians[[i]])){
        ## g = replicates        
        for(g in seq_along(all.data.medians[[i]][[h]])){
            rep.names <- vector()
            rep.names[g] <- all.set[[h]][[g]]@description$"TUBE NAME"
            names(all.data.medians[[i]][[h]]) <- rep.names
            names(all.data.medians[[i]][[h]][[g]]) <- names(all.groups.e[[1]][[1]])
        }}}


## create a similar list that gets the mean of the replicate medians
## you'll use this later for plotting and statistics
## in the first step, bind each replicate into a single list that represents
## a strain
all.data.median.lines <- lapply(all.data.medians, function(x){
                                  lapply(1:length(x), function(y){
                                             colMeans(do.call("rbind", x[[y]]))
                                         })
                              })

for(i in seq_along(all.data.median.lines)){
    names(all.data.median.lines[[i]]) <- strain.names
}


## now bind the list of strains together via 'rbind'
all.data.median.lines <- lapply(all.data.median.lines, function(x){
                                  as.data.frame(do.call("rbind", x))
                              })

all.data.medians <- lapply(all.data.medians, function(x){
                             lapply(x, function(y){
                                        as.data.frame(do.call("rbind", y))
                                    })
                         })

all.data.medians <- lapply(all.data.medians, function(x){
                             do.call("rbind", x)
                         })


strain.rep.factor <- vector(mode = "list",
                            length = length(all.data.e))
strain.rep.factor <- lapply(1:length(all.data.e), function(x){
                                strain.rep.factor[[x]] <- factor(x = unlist(lapply(1:length(all.data.e[[x]]), function(y){
rep(y, times = length(all.data.e[[x]][[y]]))
})), labels = strain.names)
})


for(i in 1:length(all.data.medians)){
    all.data.medians[[i]]$strain <- strain.rep.factor[[i]]
}




## now plot across gates/parameters
## set up some dirs first
setwd(results.dir)
dir.create(path = paste0(results.dir, "/between_groups_median_stripcharts"))
median.stripcharts.dir <- paste0(results.dir, "/between_groups_median_stripcharts")
setwd(median.stripcharts.dir)


lapply(1:length(all.data.medians), function(x){
           lapply(1:(ncol(all.data.medians[[x]]) - 1), function(y){
                      pdf(file = paste0(names(all.data.medians[x]), "_",
                                        names(all.data.medians[[x]][y]),
                                        "_median_betw_groups_strip.pdf"),
                          height = 7,
                          width = 7,
                          bg = "transparent")
                      ## 'dummy' plot that we'll overplot on 
                      stripchart(all.data.medians[[x]][, y] ~ all.data.medians[[x]]$strain,
                                 vertical = T,
                                 pch = NA,
                                 cex.axis = 0.8,
                                 main = names(all.data.medians[x]),
                                 ylab = gsub("_", " ", names(all.data.medians[[x]][y])))
                      ## add mean lines
                      lapply(1:length(strain.names), function(q){
                                 lines(x = c(q - 0.25, q + 0.25),
                                       y = rep(all.data.median.lines[[x]][q, y], 2),
                                       lwd = 2.5,
                                       col = gray(0))})
                      ## overplot medians
                      stripchart(all.data.medians[[x]][, y] ~ all.data.medians[[x]]$strain,
                                 vertical = T,
                                 cex.axis = 0.8,
                                 pch = 21,
                                 lwd = 1.25,
                                 col = gray(0.2),
                                 bg = gray(0.8),
                                 add = T,
                                 method = "jitter",
                                 jitter = 0.1,
                                 cex = 1.25)
                      dev.off()
                             })
                  })


## -----
## <<Summary_Tables>>
## [x]
## write the data we use for statistics to tables
## I format the decimal places to 3 to make reading
## the output easier.  This is accomplished via
## 'sprintf', which converts to character, so make
## separate objects for writing table output 
setwd(tables.dir)
all.data.means.round <- all.data.means
for(i in seq_along(all.data.means.round)){
        for(h in seq_along(all.data.means.round[[i]])){
            if(!is.factor(all.data.means.round[[i]][, h]))
                all.data.means.round[[i]][, h] <- sprintf("%.3f", all.data.means.round[[i]][, h])
}}


## group means table
Map(f = function(x, name){
        write.table(x = x,
                    file = paste0(name, "_group_means_table.csv"),
                    append = F,
                    quote = F,
                    sep = ",",
                    row.names = T)},
        x = all.data.mean.lines,
        name = names(all.data.mean.lines)
)

lapply(seq_along(all.data.means.round), function(x){
           write.table(x = all.data.means.round[[x]], 
                       file = paste0(names(all.data.means.round[x]),
                                   "_sample_means_table.csv"),
                       ## necessary to prevent empty column
                       ## where rownames inserted 
                       col.names = NA,
                       row.names = T,
                       ## don't enclose everything in double quotes
                       quote = F,
                       sep = ",")}
       )


## same operation but for medians 
all.data.medians.round <- all.data.medians
setwd(tables.dir)
all.data.medians.round <- all.data.medians
for(i in seq_along(all.data.medians.round)){
        for(h in seq_along(all.data.medians.round[[i]])){
            if(!is.factor(all.data.medians.round[[i]][, h]))
                all.data.medians.round[[i]][, h] <- sprintf("%.3f", all.data.medians.round[[i]][, h])
}}


lapply(seq_along(all.data.medians.round), function(x){
           write.table(x = all.data.medians.round[[x]],
                       file = paste0(names(all.data.medians.round[x]),
                                     "_sample_medians_table.csv"),
                       ## necessary to prevent empty column
                       ## where rownames inserted 
                       col.names = NA,
                       row.names = T,
                       ## don't enclose everything in double quotes
                       quote = F,
                       sep = ",")}
       )

## group medians table
Map(f = function(x, name){
        write.table(x = x,
                    file = paste0(name, "_group_medians_table.csv"),
                    append = F,
                    quote = F,
                    sep = ",",
                    row.names = T)},
        x = all.data.median.lines,
        name = names(all.data.median.lines)
)




##-----
## <<Flow_Cytometry_Statistics>>
## [x]
## set up directories
setwd(base.dir)
dir.create(path = paste0(base.dir, "/statistical_analysis"))
stats.dir <- paste0(base.dir, "/statistical_analysis")
setwd(stats.dir)

lapply(1:length(all.data.means), function(x){
           lapply(1:(ncol(all.data.means[[x]])), function(y){
                      if(!is.factor(all.data.means[[x]][, y])){
                          mean.aov <- aov(all.data.means[[x]][, y] ~
                                              all.data.means[[x]]$strain)
                          s.mean   <- capture.output(summary(mean.aov))
                          ptest    <- capture.output(PostHocTest(x = mean.aov,
                                                                 method = "lsd",
                                                                 conf.level = 0.95))
                      cat(c("\n\n-----",
                            paste0(names(all.data.means[[x]][y])),
                            "-----\n",
                            s.mean, ptest),
                          file = paste0(names(all.data.means[x]),
                                        "_means_statistics.txt"),
                          sep = "\n",
                          append = T)
                      }
                  })
       })


lapply(1:length(all.data.medians), function(x){
           lapply(1:(ncol(all.data.medians[[x]])-1), function(y){
                      if(!is.factor(all.data.means[[x]][, y])){
                          median.aov <- aov(all.data.medians[[x]][, y] ~
                                                all.data.medians[[x]]$strain)
                          s.median   <- capture.output(summary(median.aov))
                          ptest      <- capture.output(PostHocTest(x = median.aov,
                                                                   method = "lsd",
                                                                   conf.level = 0.95))
                      cat(c("\n\n-----",
                            paste0(names(all.data.medians[[x]][y])),
                            "-----\n",
                            s.median, ptest),
                          file = paste0(names(all.data.medians[x]),
                                        "_medians_statistics.txt"),
                          sep = "\n",
                          append = T)
                          }
                  })
       })


##-----
## <<Cell_Count_Stripcharts>> 
## [x]
## last thing is the replicate cell count plots
## put these into the stripchart format I used above
## this uses the object 'strain.rep.factor' that isn't 
## created until late in the code, so this piece of analysis
## ends up here. 
setwd(cell.count.dir)
individual.cell.counts <- lapply(all.data.e, function(x){
                                     lapply(x, function(y){
                                                lapply(y, function(q){
                                                           nrow(q)
                                                       })
                                            })
                                 })

individual.cell.counts <- lapply(individual.cell.counts, function(x){
                                     unlist(lapply(x, function(y){
                                                       unlist(y)
                                                   }))
                                 })

individual.cell.dataframe <- Map(f = function(x, y){data.frame(count = x, strain = y)},
                                 x = individual.cell.counts, y = strain.rep.factor)

individual.mean.lines <- lapply(individual.cell.dataframe, function(x){
                                    unlist(lapply(strain.names, function(y){
                                                      mean(x[x$strain == y, 1])
                                                  }))
                                })

lapply(1:length(individual.cell.dataframe), function(x){
           pdf(file = paste0(names(individual.cell.dataframe[x]),
                             "_individual_cell_count_strip.pdf"),
               height = 7, width = 7, bg = "transparent")
           ## dummy chart
           stripchart(count ~ strain, data = individual.cell.dataframe[[x]],
                      vertical = T, cex.axis = 0.8, pch = NA, lwd = 1.25,
                      method = "jitter", jitter = 0.1, ylab = "Cell Count", 
                      main = names(individual.cell.dataframe[x]))
           ## plot the mean lines
           lapply(1:length(strain.names), function(q){
                                 lines(x = c(q - 0.25, q + 0.25),
                                       y = rep(individual.mean.lines[[x]][q], 2),
                                       lwd = 2.5, col = gray(0))})
           ## overplot the actual values
           stripchart(count ~ strain,
                      data = individual.cell.dataframe[[x]],
                      vertical = T,
                      cex.axis = 0.8,
                      pch = 21,
                      lwd = 1.25,
                      col = gray(0.2),
                      bg = gray(0.8),
                      cex = 1.25,
                      method = "jitter",
                      jitter = 0.1, add = T)
           dev.off()
       })
