
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

# Visualize deformations from the mean
par(mfrow = c(1,2), mar = c(2,2,2,2))
plotRefToTarget(ref, shape_min, method = "TPS", mag = 2,
                main = paste0("Empirical min on PC", pc))
plotRefToTarget(ref, shape_max, method = "TPS", mag = 2,
                main = paste0("Empirical max on PC", pc))
par(mfrow = c(1,1))

mtext("Thin-plate spline deformation from Procrustes PC1", side = 3, outer = TRUE, cex = 0.9)
par(op)





# choose a PC:
pc <- 2

# indices of min/max scorers on that PC
i_min <- which.min(pca_sym$x[, pc])
i_max <- which.max(pca_sym$x[, pc])

# corresponding shapes (landmarks) from your array
shape_min <- symDat_3D$coords[ , , i_min]
shape_max <- symDat_3D$coords[ , , i_max]

# Procrustes mean shape across all
ref <- mshape(symDat_3D$coords)

# Visualize deformations from the mean
par(mfrow = c(1,2), mar = c(2,2,2,2))
plotRefToTarget(ref, shape_min, method = "TPS", mag = 2,
                main = paste0("Empirical min on PC", pc))
plotRefToTarget(ref, shape_max, method = "TPS", mag = 2,
                main = paste0("Empirical max on PC", pc))
par(mfrow = c(1,1))

mtext("Thin-plate spline deformation from Procrustes PC2", side = 3, outer = TRUE, cex = 0.9)
par(op)


library(dplyr)
library(tidyr)
library(ggplot2)


# PC scores straight from the PCA
pc_df <- as.data.frame(PCA_out$PC.points)
names(pc_df) <- c("PC1","PC2")

# Bind your metadata (must be the same length/order as the coords used in PCA)
pc_df$cluster <- symDat_2D$cluster
pc_df$gait    <- symDat_2D$gait
pc_df$species <- symDat_2D$species
pc_df$clade   <- symDat_2D$clade
pc_df$mass    <- symDat_2D$mass


levels(pc_df$cluster)

# Without transparency (left)
p1 <- ggplot(data=pc_df, aes(x=PC1, group=cluster, fill=cluster)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme()
p1

# Without transparency (left)
p2 <- ggplot(data=pc_df, aes(x=PC2, group=cluster, fill=cluster)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme()
p2

