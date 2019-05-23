##-----
## [x]
## directory and package setup
## use a wider console to make reading output easier
options(width = 200)                    
packages = c("colorspace", "flowCore", "flowViz", "flowUtils", "flowStats", "flowFP", "geneplotter", "ggcyto", "DescTools")
lapply(packages, require, character.only = TRUE)
source(file = "~/Desktop/emacs/R/functions/set2frame.R")




##-----
## [x]
## user-specified options - these will change for each analysis depending on strains/reporters
##############
## USER INPUT:
############## 
base.dir       <- "~/Desktop/data/flow/2019.05.22_Thr_dsRed_unfused_TFT_flow"
setwd(base.dir)
dir.create(path = paste0(base.dir, "/fcs"))
dir.create(path = paste0(base.dir, "/results"))
dir.create(path = paste0(base.dir, "/tables"))
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
no.reporter         <- ".*untagged.*fcs"
by.Thr.dsRed.tft    <- "BY.*Thr_dsRed.*.fcs"
rm.Thr.dsRed.tft    <- "RM.*Thr_dsRed.*.fcs"
ubr1.Thr.dsRed.tft  <- "ubr1.*Thr.*dsRed.*.fcs"
doa10.Thr.dsRed.tft <- "doa10.*Thr.*dsRed.*.fcs"


##############
## USER INPUT:
##############
## for later use in plots
strain.names <- c("untagged", "BY_Thr_dsRed_TFT", "RM_Thr_dsRed_TFT", "ubr1_Thr_dsRed_TFT", "doa10_Thr_dsRed_TFT")





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
work.dir
all.strains <- list(no.reporter, by.Thr.dsRed.tft, rm.Thr.dsRed.tft, ubr1.Thr.dsRed.tft, doa10.Thr.dsRed.tft) 
all.set     <- lapply(all.strains, function(x){read.flowSet(files = NULL, path = ".", pattern = x, alter.names = T, min.limit = 1)})
str(all.set[[1]])




#################
## END USER INPUT:
#################




##-----
## [x]
## write strain/replicate groupings to a table for inspection
setwd(tables.dir)
replicates.out <- lapply(all.set, function(x){x@phenoData@data$name})
lapply(replicates.out, function(x){cat(c(x, "\n"), file = "strain_replicate_groupings.txt", append = T, sep = ", ")})




##-----
## [x]
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
              `lGFP` = log10(`eGFP.A`),
              `lRFP` = log10(`mCherry.A`),
              `TFTR` = log(`mCherry.A`/`eGFP.A`, base = 2),
              `PSVR` = log(`eGFP.A`/`mCherry.A`, base = 2))}
all.set <- lapply(all.set, fsApply, PSV.TFT.transform)




##-----
## [x]
## get the total number of cells for each flowFrame
## nrow is passed as an optional arg to fsApply here
total.cells <- lapply(all.set, fsApply, nrow)




##-----
## [x]
## 02.27.2019 try this w/ curv2Filter w/ a big bandwidth setting to grab the main cloud of cells
## we take only cells in 'area 1' (the gate), not 'rest' (the cells outside the gate)
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
## [x]
## plot the results of the pre-filter plus curv2Filter gating
## start by undoing the complicated list structure the filter operation creates
## this yields a list of flowSets
initial.curv.split <- unlist(initial.split.all)
setwd(cell.gate.dir)
xy.initial.curv.pdf <- function(x) {
    pdf(file = paste0("curv_", x@description$"TUBE NAME", "_.pdf"), height = 7, width = 7)
    print(xyplot(`SSC.A` ~ `FSC.A`, data = x,
                 filter = curv2Filter(x = "FSC.A", y = "SSC.A", bwFac = 2, gridsize = c(250,250)),
                 smooth = F))
    dev.off()
}
lapply(initial.curv.split, fsApply, xy.initial.curv.pdf)




##-----
## [x]
## now automatically gate the cells by splitting w/ curv2Filter
curv.split <- function(x){
    split(x, f = curv2Filter(x = "FSC.A", y = "SSC.A", bwFac = 2, gridsize = c(250,250)),
          flowSet = TRUE, codeflowSet = TRUE)
}
curv.set <- lapply(initial.curv.split, fsApply, curv.split)

## the output of the filtering operation is a list
## curv.set:
## 1. curv.set -> list of flowsets
##    curv.set[[1]] -> flowSet w/ 4 experiments

##    2. curv.set[[1]][[1]] -> flowSet
##       flowSet of the different curv gates         

##       3. curv.set[[1]][[1]][[1]] -> the actual flowframe of each gate ("rest" in the preceding)
##          so there could be, e.g., 3 elements, "rest", "area 1", "area 2"
##          this is what needs to be ordered based on n. cells

## use unlist to get a simpler list structure
## the result is a list of flowsets
curv.set    <- unlist(curv.set)
order.curvs <- function(x) order(fsApply(x, nrow), decreasing = T)
## this will produce a list w/ the population of interest for each experiment
## the second object is a list of the gates in numeric order
curv.order        <- lapply(curv.set, order.curvs)
curv.orig         <- lapply(curv.order, function(x) {seq(from = 1, to = length(x), by = 1)})




##-----
## [x]
## now create a logical list for subsetting - take the pop w/ 2nd most cells (1st most = all ungated)
curv.logical <- list()
for(i in 1:length(curv.order)){
    curv.logical[[i]] <- ifelse(curv.orig[[i]] == curv.order[[i]][2], T, F)
}




##-----
## [x]
## now use a for() loop to subset so we keep only 1 gate
## curv.final.frame is a list of flowframes
c.final       <- list()
c.final.frame <- list()
for (i in 1:length(curv.set)){
    c.final[[i]] <- subset(curv.set[[i]], curv.logical[[i]])
    c.final.frame[[i]] <- c.final[[i]][[1]]
}
## now do the same for the initial gate and ungated cells




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
## plot the final gating result
setwd(cell.gate.dir)
final.gating.result.pdf <- function(x) {
    pdf(file = paste0("final_gate_", x@description$"TUBE NAME", ".pdf"), height = 7, width = 7)
    print(xyplot(`SSC.A` ~ `FSC.A`, data = x, smooth = F))
    dev.off()
}
lapply(c.final.frame, final.gating.result.pdf)
setwd(base.dir)




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
              `lGFP` = logicle.trans(`eGFP.A`),
              `lRFP` = logicle.trans(`mCherry.A`))}

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
                      pdf(file = paste0(names(all.data[x]), "_", q@description$'TUBE NAME', ".pdf"), height = 7, width = 7)
                      print(xyplot(`lRFP` ~ `lGFP`, data = q,  smooth = F,
                             strip = paste0(names(all.data[x]), "_", q@description$'TUBE NAME'),
                             prepanel=function()
                             {return(list(xlim = c(0, 4), ylim = c(0, 4)))}))
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
names(all.frame.l) <- names(all.data)
for(i in 1:length(all.frame.l)){
    names(all.frame.l[[i]]) <- strain.names
}


## plot
lapply(1:length(all.frame.l), function(x){
           lapply(1:length(all.frame.l[[x]]), function(y){
                      pdf(file = paste0(names(all.frame.l[x]), "_", names(all.frame.l[[x]][y]), ".pdf"), height = 7, width = 7)
                      print(xyplot(`lRFP` ~ `lGFP`, data = all.frame.l[[x]][[y]],  smooth = F,
                                   strip = paste0(names(all.data[x]), "_", names(all.frame.l[[x]][y])),
                                   prepanel=function()
                                   {return(list(xlim = c(0, 4), ylim = c(0, 4)))}))
                      dev.off()
                  })})




##-----
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




##-----
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
## [x]
## between groups stuff
## start by grouping replicates in a list of lists
all.groups.e <- lapply(all.data.e, function(x){
                           lapply(x, function(y){
                                      x <- do.call("rbind", y)
                                  })
                       })




##-----
## [x]
## between groups plots
## need to add names to individual elements
setwd(results.dir)
dir.create(path = paste0(results.dir, "/between_groups_plots"))
between.groups.dir <- paste0(results.dir, "/between_groups_plots")
setwd(between.groups.dir)

## set up names and limits for parameters - these may have to change 
## "FSC.A"  "SSC.A"  "eGFP.A" "mCherry.A" "Time"   "lgfp"   "lmch"   "tftr"   "psvr"  
x.lab   <- names(all.groups.e[[1]][[1]])
x.min   <- c(0, 0, 0, 0, 0, 2, 2, -6, -2)
x.max   <- c(2.5e5, 2e5, 2e4, 2e4, 2.5e3, 5, 5, 2, 6)
leg.pos <- c(rep("topright", 5), "topleft", rep("topright", 3))

y.lim <- vector(mode = "list", length = 3)
for(h in 1:length(y.lim)){
    for(i in 1:length(x.lab)){
        y.lim[[h]][i] <- max(unlist(lapply(all.groups.e[[h]], function(x){
                                               max(density(x[, i])$y)
                                           })
                                    ))}}

Map(f = function(gate, name){
        for(i in 1:length(x.lab)){

            ## get y limit for density plot
            lapply(1:length(all.groups.e[[1]]), function(x){

            ## plot
            pdf(file = paste0(name, "_", x.lab[i], "_", "_between_groups", ".pdf"), height = 7, width = 7, bg = "transparent")
            plot(density(all.groups.e[[gate]][[x]][, i], adjust = 2), col = "white", ylim = c(0, y.lim[[gate]][i]), xlim = c(x.min[i], x.max[i]), xlab = x.lab[i], main = name)

            ## set up colors for the different gates
            cols <- rainbow(n = length(all.groups.e[[gate]]), s = 0.5, v = 0.8, alpha = 0.9, end = 0.7)

            ## map across gates and do so on a per-strain basis
            Map(f = function(x, y){lines(density(all.groups.e[[gate]][[x]][, i], adjust = 2), col = y)}, x = 1:length(all.groups.e[[gate]]), y = cols)

            ## legend
            legend(x = leg.pos[i], legend = strain.names, lty = 1, lwd = 5, col = cols, bg = "white")
            dev.off()

            })
        }}, gate = 1:length(all.groups.e), name = names(all.groups.e))




##-----
## [x]
## between groups w/ replicates
setwd(results.dir)
dir.create(path = paste0(results.dir, "/between_groups_w_replicates_plots"))
between.groups.replicates.dir <- paste0(results.dir, "/between_groups_w_replicates_plots")
setwd(between.groups.replicates.dir)


Map(f = function(gate, name){
        for(i in 1:length(x.lab)){

            ## get y limit for density plot
            lapply(1:length(all.groups.e[[1]]), function(x){

            ## plot
            pdf(file = paste0(name, "_", x.lab[i], "_", "_between_groups_w_replicates", ".pdf"), height = 7, width = 7, bg = "transparent")
            plot(density(all.groups.e[[gate]][[x]][, i], adjust = 2), col = "white", ylim = c(0, y.lim[[gate]][i]), xlim = c(x.min[i], x.max[i]), xlab = x.lab[i], main = name)

            ## set up colors for the different gates
            cols <- rainbow(n = length(all.groups.e[[gate]]), s = 0.5, v = 0.8, alpha = 0.9, end = 0.7)
            light.cols <- rainbow(n = length(all.groups.e[[gate]]), s = 0.7, v = 0.9, alpha = 0.5, end = 0.7)    

                ## map across gates and do so on a per-strain basis
                Map(f = function(x, y){
                        lines(density(all.groups.e[[gate]][[x]][, i], adjust = 2), col = y, lwd = 3)
                        lapply(1:length(all.data.e[[gate]][[x]]), function(n){
                                   lines(density(all.data.e[[gate]][[x]][[n]][, i], adjust = 2), col = light.cols[x], lwd = 0.8)})},
                            x = 1:length(all.groups.e[[gate]]), y = cols)

            ## legend
            legend(x = leg.pos[i], legend = strain.names, lty = 1, lwd = 5, col = cols, bg = "white")
            dev.off()

            }
        )}}, gate = 1:length(all.groups.e), name = names(all.groups.e))




##-----
## [x]
## between gates plots
setwd(results.dir)
dir.create(path = paste0(results.dir, "/between_gates_plots"))
between.gates.dir <- paste0(results.dir, "/between_gates_plots")
setwd(between.gates.dir)

Map(f = function(strain, name){
        for(i in 1:length(x.lab)){

            ## get y limit for density plot
            y.lim <- c(0, max(unlist(lapply(1:3, function(x){max(density(all.groups.e[[x]][[strain]][, i])$y)}))))

            ## plot 
            pdf(file = paste0(x.lab[i], "_", name, "_between_gates", ".pdf"), height = 7, width = 7, bg = "transparent")
            plot(density(all.groups.e[[1]][[strain]][, i], adjust = 2), col = gray(0.8), ylim = y.lim, xlim = c(x.min[i], x.max[i]), xlab = x.lab[i], main = name)

            ## set up colors for the different gates
            cols <- rainbow(n = length(all.groups.e), s = 0.5, v = 0.8, alpha = 0.9, end = 0.5)

            ## map across gates and do so on a per-strain basis
            Map(f = function(x, y){lines(density(all.groups.e[[x]][[strain]][, i], adjust = 2), col = y)}, x = 1:3, y = cols)

            ## legend
            legend(x = leg.pos[i], legend = c("curv", "initial", "ungated"), lty = 1, lwd = 5, col = cols, bg = "white")
            dev.off()

        }
    }, strain = 1:length(strain.names), name = strain.names)




##-----
## [x]
## replicate plots
## this gives a nice layout w/ replicates separately and grouped 
setwd(results.dir)
dir.create(path = paste0(results.dir, "/replicate_plots"))
replicates.dir <- paste0(results.dir, "/replicate_plots")
setwd(replicates.dir)

y.lim <- vector(mode = "list", length = 3)
for(h in 1:length(y.lim)){
    for(i in 1:length(x.lab)){
        y.lim[[h]][i] <- max(unlist(lapply(all.groups.e[[h]], function(x){
                                               max(density(x[, i])$y)
                                           })
                                    ))}}

lapply(1:length(all.groups.e), function(x){
           lapply(1:length(all.groups.e[[x]]), function(y){
                      for(i in 1:length(x.lab)){
                          pdf(file = paste0(names(all.groups.e[x]), "_", x.lab[i], "_", strain.names[y], "_", "replicate_plot.pdf"),
                              height = 7, width = 7, bg = "transparent")
                          plot(density(all.groups.e[[x]][[y]][, i], adjust = 2), xlab = x.lab[i], main = strain.names[y],
                               xlim = c(x.min[i], x.max[i]), ylim = c(0, y.lim[[x]][i]), col = "white")
                          lapply(all.data.e[[x]][[y]], function(r){lines(density(r[, i], adjust = 2), col = gray(0.8))})
                          lines(density(all.groups.e[[x]][[y]][, i], adjust = 2), col = "#AA44AA", lwd = 2)
                          dev.off()
                      }
                  })
       })




##-----
## [x]
## 2D scatter plots
setwd(results.dir)
dir.create(path = paste0(results.dir, "/2D_scatter_plots"))
scatterplot.dir <- paste0(results.dir, "/2D_scatter_plots")
setwd(scatterplot.dir)

lapply(1:length(all.groups.e), function(x){
           lapply(1:length(all.groups.e[[x]]), function(y){
                      pdf(file = paste0(names(all.groups.e[x]), "_", strain.names[y], "_GFP_mCh_scatter.pdf"), height = 7, width = 7, bg = "transparent")
                      plot(all.groups.e[[x]][[y]]$lgfp, all.groups.e[[x]][[y]]$lmch,
                           xlim = c(1, 5), ylim = c(1, 5), pch = 19, cex = 0.1,
                           col = gray(0.5, alpha = 0.5), xlab = "log10 sfGFP", ylab = "log 10 mCherry",
                           main = strain.names[y])
                      dev.off()
                  })
       })




##-----
## [x]
## between groups boxplots 
## merge exprs datasets into a single dataframe, use nrow to create a factor for each strain
## create a vector of cell counts for each gate
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


## create a factor for each gate corresponding to the individual strains
strain.factor <- vector(mode = "list", length = length(all.groups.e))

strain.factor <- lapply(all.groups.e, function(x){
                            unlist(lapply(1:length(x), function(y){
                                              rep(x = y, times = nrow(x[[y]]))
                                          }))
                                   })

strain.factor <- lapply(1:length(strain.factor), function(x){
                            factor(x = strain.factor[[x]], levels = unique(strain.factor[[x]]), labels = strain.names)
                            })

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
                                pdf(file = paste0(names(all.groups.df[x]), "_", names(all.groups.df[[x]][y]), "_", "strain_boxplot.pdf"),
                                    height = 7, width = 7, bg = "transparent")
                                par(cex.axis = 0.8) 
                                boxplot(all.groups.df[[x]][, y] ~ all.groups.df[[x]]$strain,
                                        main = names(all.groups.df[x]),
                                        xlab = names(all.groups.df[[x]][y]),
                                        col = gray(0.95),
                                        pars = list(outpch = 19, outcol = "#112255AA", outcex = 0.4))
                                dev.off()
                            })
                 }))





##-----
## [x]
## replicate mean stripcharts
## start by creating a list of lists identical in size to the individual replicates list
all.data.means <- vector(mode = "list", length = length(all.data.e))
all.data.means <- lapply(1:length(all.data.e), function(x){
                             all.data.means[[x]] <- lapply(1:length(all.data.e[[x]]), function(y){
                                                               all.data.means[[x]][[y]] <- vector(mode = "list", length = length(all.data.e[[x]][[y]]))
                                                       })
                         })


## now get the mean of each replicate, which is a sample of 1000's of cells
all.data.means <- lapply(1:length(all.data.e), function(x){
                             all.data.means[[x]] <- lapply(1:length(all.data.e[[x]]), function(y){
                                                                        lapply(1:length(all.data.e[[x]][[y]]), function(i){
                                                                                   sapply(X = 1:ncol(all.data.e[[x]][[y]][[i]]),
                                                                                          FUN = function(m){mean(all.data.e[[x]][[y]][[i]][, m])})
                                                                               })
                                                           })
                         })


## create a similar list that gets the mean of the replicate means
## you'll use this later for plotting and statistics
## in the first step, bind each replicate into a single list that represents a strain
all.data.mean.lines <- lapply(all.data.means, function(x){
                                  lapply(1:length(x), function(y){
                                             colMeans(do.call("rbind", x[[y]]))
                                         })})

## now bind the list of strains together via 'rbind'
all.data.mean.lines <- lapply(all.data.mean.lines, function(x){
                                  as.data.frame(do.call("rbind", x))
})



all.data.means <- lapply(all.data.means, function(x){
                             lapply(x, function(y){
                                        as.data.frame(do.call("rbind", y))
                                    })
                         })

names(all.data.means) <- names(all.data.e)

for(i in 1:length(all.data.means)){
    names(all.data.means[[i]]) <- strain.names
}

all.data.means <- lapply(all.data.means, function(x){
                             do.call("rbind", x)
                         })

params <- names(all.data.e[[1]][[1]][[1]])
for(i in 1:length(all.data.means)){
    names(all.data.means[[i]]) <- params
}

strain.rep.factor <- vector(mode = "list", length = length(all.data.e))
strain.rep.factor <- lapply(1:length(all.data.e), function(x){
                                strain.rep.factor[[x]] <- factor(x =  unlist(lapply(1:length(all.data.e[[x]]), function(y){
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
                          height = 7, width = 7, bg = "transparent")
                      stripchart(all.data.means[[x]][, y] ~ all.data.means[[x]]$strain,
                                 vertical = T, pch = NA,
                                 cex.axis = 0.8, main = names(all.data.means[x]),
                                 ylab = names(all.data.means[[x]][y]))
                      lapply(1:length(strain.names), function(q){
                                 lines(x = c(q - 0.25, q + 0.25),
                                       y = rep(all.data.mean.lines[[x]][q, y], 2),
                                       lwd = 2.5, col = gray(0))})
                                 stripchart(all.data.means[[x]][, y] ~ all.data.means[[x]]$strain,
                                            vertical = T, cex.axis = 0.8, pch = 21, lwd = 1.25,
                                            col = gray(0.2), bg = gray(0.8), add = T,
                                            method = "jitter", jitter = 0.1, cex = 1.25)
                      dev.off()
                             })
                  })




##-----
## [x]
## replicate median stripcharts
## start by creating a list of lists identical in size to the individual replicates list
all.data.medians <- vector(mode = "list", length = length(all.data.e))
all.data.medians <- lapply(1:length(all.data.e), function(x){
                             all.data.medians[[x]] <- lapply(1:length(all.data.e[[x]]), function(y){
                                                                 all.data.medians[[x]][[y]] <- vector(mode = "list", length = length(all.data.e[[x]][[y]]))
                                                       })
                         })


## now get the median of each replicate, which is a sample of 1000's of cells
all.data.medians <- lapply(1:length(all.data.e), function(x){
                               all.data.medians[[x]] <- lapply(1:length(all.data.e[[x]]), function(y){
                                                                 lapply(1:length(all.data.e[[x]][[y]]), function(i){
                                                                            sapply(X = 1:ncol(all.data.e[[x]][[y]][[i]]),
                                                                                   FUN = function(m){median(all.data.e[[x]][[y]][[i]][, m])})
                                                                               })
                                                           })
                         })


## create a similar list that gets the mean of the replicate means
## you'll use this later for plotting and statistics
## in the first step, bind each replicate into a single list that represents a strain
all.data.median.lines <- lapply(all.data.medians, function(x){
                                  lapply(1:length(x), function(y){
                                             colMeans(do.call("rbind", x[[y]]))
                                         })})

## now bind the list of strains together via 'rbind'
all.data.median.lines <- lapply(all.data.median.lines, function(x){
                                  as.data.frame(do.call("rbind", x))
})



all.data.medians <- lapply(all.data.medians, function(x){
                             lapply(x, function(y){
                                        as.data.frame(do.call("rbind", y))
                                    })
                         })

names(all.data.medians) <- names(all.data.e)

for(i in 1:length(all.data.medians)){
    names(all.data.medians[[i]]) <- strain.names
}

all.data.medians <- lapply(all.data.medians, function(x){
                               do.call("rbind", x)
                         })

params <- names(all.data.e[[1]][[1]][[1]])
for(i in 1:length(all.data.medians)){
    names(all.data.medians[[i]]) <- params
}

strain.rep.factor <- vector(mode = "list", length = length(all.data.e))
strain.rep.factor <- lapply(1:length(all.data.e), function(x){
                                strain.rep.factor[[x]] <- factor(x =  unlist(lapply(1:length(all.data.e[[x]]), function(y){
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
                          height = 7, width = 7, bg = "transparent")
                      stripchart(all.data.medians[[x]][, y] ~ all.data.medians[[x]]$strain,
                                 vertical = T, pch = NA,
                                 cex.axis = 0.8, main = names(all.data.medians[x]),
                                 ylab = names(all.data.medians[[x]][y]))
                      lapply(1:length(strain.names), function(q){
                                 lines(x = c(q - 0.25, q + 0.25),
                                       y = rep(all.data.median.lines[[x]][q, y], 2),
                                       lwd = 2.5, col = gray(0))})
                                 stripchart(all.data.medians[[x]][, y] ~ all.data.medians[[x]]$strain,
                                            vertical = T, cex.axis = 0.8, pch = 21, lwd = 1.25,
                                            col = gray(0.2), bg = gray(0.8), add = T,
                                            method = "jitter", jitter = 0.1, cex = 1.25)
                      dev.off()
                             })
                  })




##-----
## [x]
## statistics
## set up directories


setwd(base.dir)
dir.create(path = paste0(base.dir, "/statistical_analysis"))
stats.dir <- paste0(base.dir, "/statistical_analysis")
setwd(stats.dir)

names(all.data.means[[1]])
all.data.means[[1]][1:10, ]



params <- c(params[1:(length(params)-1)], "nl_tftr", "strain")

for(i in 1:length(all.data.means)){
    all.data.means[[i]]$psvr <- all.data.means[[i]]$mCherry.A/all.data.means[[i]]$eGFP.A
    names(all.data.means[[i]]) <- params
}

all.data.means
lapply(1:length(all.data.means), function(x){
           lapply(1:(ncol(all.data.means[[x]])-1), function(y){
                      mean.aov <- aov(all.data.means[[x]][, y] ~ all.data.means[[x]]$strain)
                      s.mean   <- capture.output(summary(mean.aov))
                      ptest    <- capture.output(PostHocTest(x = mean.aov, method = "lsd", conf.level = 0.95))
                      cat(c("\n\n-----",
                            paste0(names(all.data.means[[x]][y])),
                            "-----\n",
                            s.mean, ptest), file = paste0(names(all.data.means[x]), "_means_statistics.txt") , sep = "\n",
                          append = T)
                  })
       })




for(i in 1:length(all.data.means)){
    all.data.medians[[i]]$psvr <- all.data.medians[[i]]$mCherry.A/all.data.medians[[i]]$eGFP.A
    names(all.data.medians[[i]]) <- params
}


lapply(1:length(all.data.medians), function(x){
           lapply(1:(ncol(all.data.medians[[x]])-1), function(y){
                      median.aov <- aov(all.data.medians[[x]][, y] ~ all.data.medians[[x]]$strain)
                      s.median   <- capture.output(summary(median.aov))
                      ptest      <- capture.output(PostHocTest(x = median.aov, method = "lsd", conf.level = 0.95))
                      cat(c("\n\n-----",
                            paste0(names(all.data.medians[[x]][y])),
                            "-----\n",
                            s.median, ptest), file = paste0(names(all.data.medians[x]), "_medians_statistics.txt") , sep = "\n",
                          append = T)
                  })
       })




##-----
## [x]
## last thing is the replicate cell count plots
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

individual.cell.dataframe <- Map(f = function(x, y){data.frame(count = x, strain = y)}, x = individual.cell.counts, y = strain.rep.factor)

individual.mean.lines <- lapply(individual.cell.dataframe, function(x){
                                    unlist(lapply(strain.names, function(y){
                                                      mean(x[x$strain == y, 1])
                                                  }))
                                })

lapply(1:length(individual.cell.dataframe), function(x){

           pdf(file = paste0(names(individual.cell.dataframe[x]), "_individual_cell_count_strip.pdf"),
               height = 7, width = 7, bg = "transparent")
           
           stripchart(count ~ strain, data = individual.cell.dataframe[[x]],
                      vertical = T, cex.axis = 0.8, pch = NA, lwd = 1.25,
                      method = "jitter", jitter = 0.1, ylab = "Cell Count", 
                      main = names(individual.cell.dataframe[x]))
           
           lapply(1:length(strain.names), function(q){
                                 lines(x = c(q - 0.25, q + 0.25),
                                       y = rep(individual.mean.lines[[x]][q], 2),
                                       lwd = 2.5, col = gray(0))})

           stripchart(count ~ strain, data = individual.cell.dataframe[[x]],
                      vertical = T, cex.axis = 0.8, pch = 21, lwd = 1.25,
                      col = gray(0.2), bg = gray(0.8), cex = 1.25,
                      method = "jitter", jitter = 0.1, add = T)

           dev.off()

       })

