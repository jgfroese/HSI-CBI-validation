########################################################
## HSI-CBI - MODEL VALIDATION WITH PRESENCE-ONLY DATA ##
########################################################

##############################
## Load required R packages ##
##############################

require(ecospat)
require(data.table)
require(zoo)

######################
## Data preparation ##
######################

## For spatial data preparation in package {raster} or alternative GIS software see "README.md"

## Step 1: compute model-predicted values (e.g. HSI) at species presence records
# read exported .TXT file for each model (m1,m2,...) / validation presences (p1,p2,...) combination as data frame
pred.m1.p1 = read.csv("[...]/your-combined-model1-presences1.txt")
pred.m1.p2 = read.csv("[...]/your-combined-model1-presences2.txt")
...
pred.m2.p1 = read.csv("[...]/your-combined-model2-presences1.txt")
pred.m2.p2 = read.csv("[...]/your-combined-model2-presences2.txt")
...

# combine all data frames in list
pred.list <- list (pred.m1.p1, pred.m1.p2, ..., pred.m2.p1, pred.m2.p2, ...)

# homogenise 4 column names (if different in exported .TXT files)
for (i in seq_along(pred.list)) {
  setnames(pred.list[[i]], c("ID", "Pixelcount", "HSI", "Presences"))
}

# calculate total number of presence records per HSI value (multiple presences may fall within one pixel)
n.list <- your.list.length # the number of data frames in your list
pred.sum.list <- vector("list", n.list)
for (i in seq_along(pred.list)) {
  pred.sum.list[[i]] <- aggregate(cbind(Pixelcount*Presences)~HSI, data = pred.list[[i]], sum)
}

# add descriptive column names and
for (i in seq_along(pred.sum.list)) {
  setnames(pred.sum.list[[i]], c("HSI", "Presences"))
}

# convert into list of vectors (= HSI at species presence records)
pred.v.list = vector("list", n.list)
for (i in seq_along(pred.sum.list)) {
  pred.v = vector()
  for (j in 1:length(pred.sum.list[[i]]$Presences)) {
    for (k in 1:pred.sum.list[[i]][j, 2]) {
      pred.v <-append(pred.v, pred.sum.list[[i]][j, 1])
    }
  }
  pred.v.list[[i]] <- append(pred.v.list[[i]], pred.v)
}

## Step 2: compute model-predicted values (e.g. HSI) across validation background
# read exported .TXT file for each model (m1,m2,...) / validation background (b1,b2,...) combination as data frame
exp.m1.b1 = read.csv("[...]/your-combined-model1-background1.txt")
exp.m1.b2 = read.csv("[...]/your-combined-model1-background2.txt")
...
exp.m2.b1 = read.csv("[...]/your-combined-model2-background1.txt")
exp.m2.b2 = read.csv("[...]/your-combined-model2-background2.txt")
...

# combine all data frames in list
exp.list <- list (exp.m1.b1, exp.m1.b2, ..., exp.m2.b1, exp.m2.b2, ...)

# homogenise 4 column names (if different in exported .TXT files)
for (i in seq_along(exp.list)) {
  setnames(exp.list[[i]], c("ID", "HSI", "Pixelcount"))
}

# convert into list of vectors (= HSI across validation background [potentially very large!])
exp.v.list = vector("list", n.list)
for (i in seq_along(exp.list)) {
  exp.v = vector()
  for (j in 1:length(exp.list[[i]]$Pixelcount)) {
    for (k in 1:exp.list[[i]][j, 2]) {
      exp.v <-append(exp.v, exp.list[[i]][j, 1])
    }
  }
  exp.v.list[[i]] <- append(exp.v.list[[i]], exp.v)
}


###################
## Data analysis ##
###################

## Step 3: apply function "boyce {ecospat}" to each model (m1,m2,...) / validation data (v1,v2,...) combination with parameters:
exp.v.list[[i]] # HSI across each validation data background (m1.b1, m1.b2, ... m2.b1, m2.b2, ...)
pred.v.list[[i]] # HSI at each validation presences (m1.p1, m1.p2, ... m2.p1, m2.p2, ...)
nclass = 0 # defaults to moving window (continuous, classification-independent) computation with arguments
window.w = 10 # moving window width (i.e. 10 adjacent HSI values are considered in each computation)
res = 100 # resolution factor (i.e. 100 computations across model-predicted value range)
PEplot = F # no PEplot is generated (customised plotting below)
boyce.m1.v1 <- ecospat.boyce(exp.v.list[[1]], pred.v.list[[1]], nclass, window.w, res, PEplot)
boyce.m1.v2 <- ecospat.boyce(exp.v.list[[2]], pred.v.list[[2]], nclass, window.w, res, PEplot)
...
boyce.m2.v1 <- ecospat.boyce(exp.v.list[[3]], pred.v.list[[3]], nclass, window.w, res, PEplot)
boyce.m2.v2 <- ecospat.boyce(exp.v.list[[4]], pred.v.list[[4]], nclass, window.w, res, PEplot)

## Step 4: investigate results of CBI analysis 
# combine all results in list
boyce.list <- list (boyce.m1.v1, boyce.m1.v2, ..., boyce.m2.v1, boyce.m2.v2, ...)

# print out CBI ($Spearman.cor) and PEmax (highest computed predicted-to-expected ratio)
CBI.list = vector("list", n.list); PEmax.list = vector("list", n.list);
for (i in seq_along(boyce.list)) {
  CBI.list[[i]] <- append(CBI.list[[i]],boyce.list[[i]]$Spearman.cor);
  PEmax.list[[i]] <- append(PEmax.ist[[i]],max(boyce.list[[i]]$F.ratio, na.rm = TRUE))
}
CBI.list; PEmax.list

## Step 5: investigate model specificity: compute proportion of validation background with HSI above threshold
t.HSI <- your.HSI.threshold # the user-defined HSI threshold (can be defined following investigation of PE-rato plots below)
HSI.t <- vector("list", n.list) # replace "t" with user-defined threshold to maintain transparency
for (i in seq_along(HSI.t)) {
  HSI.t[[i]] <- aggregate(Pixelcount~HSI > 59.99, data = exp.list[[i]], sum) / sum(exp.list[[i]]$Pixelcount)
}
HSI.t.list <- vector("list", n.list)
for (i in seq_along(HSI.t.list)) {
  HSI.t.list[[i]] <- append(HSI.t.list[[i]], (round(HSI.t[[i]][2, "Pixelcount"] * 100, digits = 0)))
}
HSI.t.list 


######################
## Validation plots ##
######################

## Option 1: plot P/E ratios of each each model (m1,m2,...) / validation data (v1,v2,...) combination separately
graphics.off()
for (i in seq_along(boyce.list)) {
  plot(boyce.list[[i]]$HS, boyce.list[[i]]$F.ratio, type = "n", 
       xlab = "Habitat suitability index (HSI)", ylab = "Predicted-to-expected (P/E) ratio",
       xlim = c(10, 90), ylim = c(0, 5), main = paste("boyce.list", i)) # adjust xlim/ylim to your data # type = "n" omits points
  lines(boyce.list[[i]]$HS, boyce.list[[i]]$F.ratio, col = "black") # draws line instead of points
  av20 <- rollmean (na.fill(boyce.list[[i]]$F.ratio, 0), 20, fill = "extend") # smoothed average across 20 PE ratio points
  lines(boyce.list[[i]]$HS, av20, col = "red", lwd = 2) # trend line (from smoothed average)
}

## Option 2: compare P/E ratios between model predictions or validation data sets by plotting from a nested list list
# (a) if comparing different validation data sets for the same model prediction, nest lists as follows:
boyce.m1 <- list (boyce.m1.v1, boyce.m1.v2, ...)
boyce.m2 <- list (boyce.m2.v1, boyce.m2.v2, ...)
plot.m.list <- list (model_1 = boyce.m1, model_2 = boyce.m2, ...)

# (b)if comparing different model predictions for the same validation data sets, nest lists as follows:
boyce.v1 <- list (boyce.m1.v1, boyce.m2.v1, ...)
boyce.v2 <- list (boyce.m1.v2, boyce.m2.v2, ...)
plot.v.list <- list (validation_1 = boyce.v1, validation_2 = boyce.v2, ...) # give each result in the list a name

# create plot (adjust graphical pars for your data and length of plot.list, here shown for (a) plot.m.list with length = 4)
graphics.off()
par(mfrow = c(2, 3), mar = c(2, 2, 2, 1), oma = c(3, 3, 0, 0), cex = 0.6)
for (i in seq_along(plot.m.list)) {
  plot(plot.m.list[[i]][[1]]$HS, plot.m.list[[i]][[1]]$F.ratio, type = "n", xlab = '', ylab = '', 
       xlim = c(10, 90), ylim = c(0, 5), main = paste(names(boyce.list[i]))) # adjust xlim/ylim to your data
  lines(boyce.list[[i]][[1]]$HS, boyce.list[[i]][[1]]$F.ratio, col = "black", lty = 1) # one line per validation data set
  lines(boyce.list[[i]][[2]]$HS, boyce.list[[i]][[2]]$F.ratio, col = "black", lty = 2)
  av20 <- rollmean((na.fill(boyce.list[[i]][[1]]$F.ratio, 0) + na.fill(boyce.list[[i]][[2]]$F.ratio, 0)) / 2, 20,
       fill = "extend") # smoothed average across all validation data sets (20 PE ratio points)
  lines((boyce.list[[i]][[1]]$HS + boyce.list[[i]][[2]]$HS) / 2, av20, col = "red",lwd = 2) # trend line
}
plot.new()
legend("center", c("name of validation data 1", "name of validation data 1", "Smoothed average"), 
       bty = "n", lty = c(1,2,1), col = c("black", "black", "red"), lwd = c(1, 1, 3), cex = 1.5)
mtext("Habitat suitability index (HSI)", side = 1, outer = TRUE, cex = 1.2, line = 1.3)
mtext("Predicted-to-expected (P/E) ratio", side = 2, outer = TRUE, cex = 1.2, line = 1.1)