library(geomorph)
# assuming symDat_3D$coords is p × k × n (e.g., 19 × 2 × 7400)

## Symmetrical Datasets ####
symDat_2D <- readRDS("./data/symGaits_2D_GDF.RDS") 
symDat_3D <- symDat_2D
symDat_3D$coords <- arrayspecs(symDat_2D$coords, 19, 2)

# Make Base PCA Alignment & Plots

dimnames(symDat_3D$coords)[[3]] <- as.character(seq(1:7400))  # rows (optional)
dimnames(symDat_3D$coords)[[2]] <- c("X","Y")            # columns
#dimnames(symDat_3D$coords)[[1]] <- as.character(seq_len(nrow(out_3D)))  # rows (optional)

#


pca_sym <- gm.prcomp(symDat_3D$coords)
PCA_out <- plot(pca_sym)
# PCA you already did:
# pca_sym <- gm.prcomp(symDat_3D$coords)

# choose a PC:
pc <- 1

# indices of min/max scorers on that PC
i_min <- which.min(pca_sym$x[, pc])
i_max <- which.max(pca_sym$x[, pc])

# corresponding shapes (landmarks) from your array
shape_min <- symDat_3D$coords[ , , i_min]
shape_max <- symDat_3D$coords[ , , i_max]

# Procrustes mean shape across all
ref <- mshape(symDat_3D$coords)

par(mfrow = c(1,3), mar = c(2,2,3,1))
plotRefToTarget(ref, shape_min, method = "TPS", mag = 2)
title(sprintf("PC%d empirical MIN (score=%.3f)", pc, pca_sym$x[i_min, pc]))

plotRefToTarget(ref, ref, method = "TPS", mag = 2)
title("Mean shape")

plotRefToTarget(ref, shape_max, method = "TPS", mag = 2)
title(sprintf("PC%d empirical MAX (score=%.3f)", pc, pca_sym$x[i_max, pc]))


install.packages('gifski')
library(gifski)


Y   <- symDat_3D$coords
ref <- mshape(Y)
p   <- dim(Y)[1]; k <- dim(Y)[2]

pc      <- 3
i_min   <- which.min(pca_sym$x[, pc])
i_max   <- which.max(pca_sym$x[, pc])
shape_min <- Y[,, i_min]
shape_max <- Y[,, i_max]

steps <- 41
w_seq <- seq(0, 1, length.out = steps)   # 0 = min, 1 = max

png_dir   <- tempdir()
png_files <- file.path(png_dir, sprintf("specinterp_%03d.png", seq_along(w_seq)))

for (i in seq_along(w_seq)) {
  w      <- w_seq[i]
  target <- (1 - w) * shape_min + w * shape_max
  
  png(png_files[i], width=900, height=900, res=150); par(mar=c(2,2,3,1))
  ok <- TRUE
  if (k == 2) {
    ok <- tryCatch({ plotRefToTarget(ref, target, method="TPS", mag=1); TRUE },
                   error=function(e) FALSE)
  }
  if (k == 3 || !ok) {
    plotRefToTarget(ref, target, method="points", mag=1, label=FALSE)
  }
  title(sprintf("Specimen morph  w=%.2f  (min→max on PC%d)", w, pc))
  dev.off()
}

gifski(png_files, gif_file = "Specimen_min_to_max_PC3.gif", width=900, height=900, delay=0.06)



#make a sequence of three images, min mean max, but remove the wireframe. 
#Make the points bold. Then join the points by lines. 
#There should be one line for points 1-9 (spine) 
#a second line for pts 10-14 (forelimb) 
#and one line for hindlimb (15-19)


# assumes: symDat_3D$coords (p x 2 x n), pca_sym, pc already defined
stopifnot(dim(symDat_3D$coords)[2] == 2)

pc <- 2
Y   <- symDat_3D$coords
ref <- mshape(Y)
i_min <- which.min(pca_sym$x[, pc])
i_max <- which.max(pca_sym$x[, pc])

shape_min <- Y[,, i_min]
shape_max <- Y[,, i_max]
shape_mean <- ref

# landmark groups
spine    <- 1:9
forelimb <- 10:14
hindlimb <- 15:19

# consistent axes across all three panels
allX <- c(shape_min[,1], shape_mean[,1], shape_max[,1])
allY <- c(shape_min[,2], shape_mean[,2], shape_max[,2])
pad  <- 0.05 * max(diff(range(allX)), diff(range(allY)))
xlim <- range(allX) + c(-pad, pad)
ylim <- range(allY) + c(-pad, pad)

draw_shape <- function(shp, title_txt) {
  plot(NA, xlim = xlim, ylim = ylim, asp = 1, xlab = "", ylab = "",
       axes = FALSE, main = title_txt)
  # lines first (so points sit on top)
  lines(shp[spine,    1], shp[spine,    2], lwd = 3)
  lines(shp[forelimb, 1], shp[forelimb, 2], lwd = 3)
  lines(shp[hindlimb, 1], shp[hindlimb, 2], lwd = 3)
  # bold points
  points(shp[,1], shp[,2], pch = 21, bg = "black", col = "black", cex = 1.4)
}

op <- par(mfrow = c(1,3), mar = c(1.5,1.5,3,0.5))
on.exit(par(op), add = TRUE)

draw_shape(shape_min,  sprintf("PC%d MIN (score = %.3f)", pc, pca_sym$x[i_min, pc]))
draw_shape(shape_mean, "Mean shape")
draw_shape(shape_max,  sprintf("PC%d MAX (score = %.3f)", pc, pca_sym$x[i_max, pc]))




