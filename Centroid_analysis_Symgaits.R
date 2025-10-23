library(geomorph)
library(tidyverse)



setwd('G:/Project files/josh geomorph')

## Symmetrical Datasets ####
symDat_2D <- readRDS("./data/symGaits_2D_GDF.RDS") 
dim(symDat_2D$coords)
dimnames(symDat_2D$coords)

symDat_3D <- symDat_2D
symDat_3D$coords <- arrayspecs(symDat_2D$coords, 19, 2)

dimnames(symDat_3D$coords)[[3]] <- as.character(seq(1:7400))  # rows (optional)
dimnames(symDat_3D$coords)[[2]] <- c("X","Y")            # columns
dimnames(symDat_3D$coords)[[1]] <- as.character(seq_len(nrow(out_3D)))  # rows (optional)

# Make Base PCA Alignment & Plots

pca_sym <- gm.prcomp(symDat_3D$coords)
PCA_out <- plot(pca_sym)
plot(PCA_out$PC.points, asp=1)

pca_sym$A

# Get the PC Proportions of Variance
pca_sym_symmary <- summary(pca_sym)
pca_sym_symmary$PC.summary


# 
# # pick the stride number you want (1 to 74)
# stride_num <- 3  
# 
# frames_per_stride <- 100
# start_idx <- (stride_num - 1) * frames_per_stride + 1
# end_idx   <- stride_num * frames_per_stride
# 
# this_stride <- as.data.frame(PCA_out$PC.points[start_idx:end_idx, ])
# colnames(this_stride) <- c('X','Y')
# 
# # Example 1: plot X vs Y trajectory for that stride
# plot(this_stride$X, this_stride$Y, pch=20,
#      xlab = "X", ylab = "Y",
#      main = paste("Stride", stride_num, "trajectory"))
# 
# 
# #all strides 
# 
# frames_per_stride <- 100
# n_strides <- 74
# 
# plot(NULL, xlim = range(PCA_out$PC.points[,1]),
#      ylim = range(PCA_out$PC.points[,2]),
#      xlab = "X", ylab = "Y",
#      main = "All stride trajectories",
#      asp = 1)
# 
# for (stride_num in 1:n_strides) {
#   start_idx <- (stride_num - 1) * frames_per_stride + 1
#   end_idx   <- stride_num * frames_per_stride
#   
#   this_stride <- as.data.frame(PCA_out$PC.points[start_idx:end_idx, ])
#   colnames(this_stride) <- c('X','Y')
#   
#   lines(this_stride$X, this_stride$Y, col = rgb(0,0,0,0.3))  # semi-transparent lines
#   points(this_stride$X, this_stride$Y, pch = 20, cex = 0.5, col = rgb(0,0,1,0.3))
# }
# 

##plot centroids

frames_per_stride <- 100
n_strides <- 74

# empty plot with correct ranges
plot(NULL, xlim = range(PCA_out$PC.points[,1]),
     ylim = range(PCA_out$PC.points[,2]),
     xlab = "X", ylab = "Y",
     main = "Stride trajectories and centroids",
     asp = 1)

centroids <- matrix(NA, nrow = n_strides, ncol = 2)

for (stride_num in 1:n_strides) {
  start_idx <- (stride_num - 1) * frames_per_stride + 1
  end_idx   <- stride_num * frames_per_stride
  
  this_stride <- as.data.frame(PCA_out$PC.points[start_idx:end_idx, ])
  colnames(this_stride) <- c("X","Y")
  
  # plot the stride path
  lines(this_stride$X, this_stride$Y, col = rgb(0,0,0,0.2))
  
  # compute centroid
  centroids[stride_num, ] <- colMeans(this_stride)
}

# plot centroids as red points
points(centroids[,1], centroids[,2], pch = 19, col = "red", cex = 1.2)

sng <- seq(1,7400,100)


length(symDat_2D$file[sng])


bigdat<-data.frame(PC1=centroids[,1], 
                   PC2=centroids[,2], 
                   clade=symDat_2D$clade[sng], 
                   species=symDat_2D$species[sng], 
                   cluster=symDat_2D$cluster[sng],
                   gait=symDat_2D$gait[sng],
                   mass=symDat_2D$mass[sng])

write.csv(bigdat, 'bigdat_Symmetrical_v2.csv')


bigdat <- bigdat[-which(bigdat$species == 'Thylacinus cynocephalus'),]

table(bigdat$gait)

library(vegan)
adonis2(cbind(bigdat$PC1, bigdat$PC2) ~ clade + species + cluster + gait + mass,
        data = bigdat, method = "euclidean", by = "terms")

adonis2(cbind(bigdat$PC1, bigdat$PC2) ~ cluster + mass + gait,
        data = bigdat, method = "euclidean", by = "terms")



library(lme4)
library(car)
library(effectsize)
library(performance)

# PC1 model
#m1 <- lmer(PC1 ~ cluster + mass + gait + (1|clade/species), data = bigdat)
m1 <- lmer(PC1 ~ cluster + mass + gait + (1|species), data = bigdat)

# PC2 model
#m2 <- lmer(PC2 ~ cluster + mass + gait + (1|clade/species), data = bigdat)
m2 <- lmer(PC2 ~ cluster + mass + gait + (1|species), data = bigdat)

# Type II tests for fixed effects
car::Anova(m1, type = 2)
car::Anova(m2, type = 2)

# Partial eta² for effect sizes
eta_squared(m1, partial = TRUE)
eta_squared(m2, partial = TRUE)

# Marginal vs conditional R² (variance explained by fixed vs fixed+random)
r2_nakagawa(m1)
r2_nakagawa(m2)


library(emmeans)

# Optional (slightly better small-sample dfs for lmer):
# install.packages("pbkrtest")
emm_options(lmer.df = if (requireNamespace("pbkrtest", quietly = TRUE)) "kenward-roger" else "satterthwaite")

# Set how emmeans averages over factors:
## - weights="equal" (default): equal-weight over levels of gait
## - weights="proportional": weight by sample size per gait level (often preferable if unbalanced)
wts <- "proportional"

# Choose a reference mass (mean) for the marginal means
mbar <- mean(bigdat$mass, na.rm = TRUE)

# ---- PC1: pairwise cluster comparisons ----
emm1 <- emmeans(m1, ~ cluster, weights = wts, at = list(mass = mbar))
pc1_pairs <- pairs(emm1, adjust = "tukey")           # Tukey all-pairs tests
pc1_pairs

emm1 <- emmeans(m1, ~ gait, weights = wts, at = list(mass = mbar))
pc1_pairs <- pairs(emm1, adjust = "tukey")           # Tukey all-pairs tests
pc1_pairs


# ---- PC2: pairwise cluster comparisons ----
emm2 <- emmeans(m2, ~ cluster, weights = wts, at = list(mass = mbar))
pc2_pairs <- pairs(emm2, adjust = "tukey")
pc2_pairs

emm2 <- emmeans(m2, ~ gait, weights = wts, at = list(mass = mbar))
pc2_pairs <- pairs(emm2, adjust = "tukey")
pc2_pairs


boxplot(PC1 ~ cluster, bigdat)
boxplot(PC2 ~ cluster, bigdat)

