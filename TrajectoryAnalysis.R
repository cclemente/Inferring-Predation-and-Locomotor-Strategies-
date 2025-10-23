# Load Packages ####
library(geomorph)
library(gtsummary)
library(viridisLite)
library(broom)
library(scales)
library(openxlsx)
source("./functions/tidy_trajectory_summaries_function.R")
source("./functions/species_average_gdf_function.R")
source("./functions/plot_pc_shape_extremes_from_TA_function.R")

# Load Data ####

## Symmetrical Datasets ####
symDat_2D <- readRDS("./data/symGaits_2D_GDF.RDS") 
dim(symDat_2D$coords)
symDat_2D_walkOnly <- readRDS("./data/symGaits_2D_GDF_walkOnly.RDS")
symDat_2D_walkOnly_noThy <- readRDS("./data/symGaits_2D_GDF_walkOnly_noThy.RDS")

symDat_2D_noThy <- readRDS("./data/symGaits_2D_GDF_noThy.RDS")

## Asymetrical Datasets ####
aSymDat_2D <- readRDS("./data/aSymGaits_2D_GDF.RDS")

## Inspect Data ####
summary(symDat_2D)
summary(symDat_2D_walkOnly)
summary(symDat_2D_noThy)
summary(aSymDat_2D)

unique(aSymDat_2D$cluster)
unique(aSymDat_2D$gait)

# Check to see that $coords is 2D
dim(symDat_2D$coords)
dim(symDat_2D_walkOnly$coords)
dim(symDat_2D_noThy$coords)
dim(aSymDat_2D$coords)

##  Build Color Palattes ####
unique(symDat_2D$cluster)

col_full <- c("#703017","#f36a3c","#179172" ,"#fbd301","#6d57cf","#c5b351")
col_full_noWalk <- c("#f36a3c","#179172", "#fbd301", "#6d57cf", "#c5b358")
col_full_clade <- c("#440154FF","#31688EFF","#35B779FF","#FDE725FF")

cols <- viridis(4)
col_fade <- c("#70301777","#f36a3c77","#17917277", "#fbd30177", "#6d57cf77", "#c5b35877")
col_fade_clade <- c("#AA88AAFF","#98B3C7FF","#9ADBC0FF","#FEF38FFF")

col_fade_noWalk <- c("#f36a3c77","#17917277","#c5b35877","#6d57cf77","#fbd30177")


# Multivariate Regression


test <- procD.lm(coords ~ cluster * clade * frame , data = symDat_2D)
test
summary(test)


# View Raw PCA Info ####

symm_PCA_space <- gm.prcomp(symDat_2D$coords, GLS = TRUE)
sym_PCA_sym <- summary(symm_PCA_space)

basefit <- gm.prcomp(symDat_2D$coords)
plot(basefit)$
baseFitSum <- summary(basefit)
baseFitSum

basePlot <- plot(basefit)

# General Models - Exploratory ####



## Model considering shape Change ####
# Try Polynomial

dat <- symDat_2D
dat$frame <- as.numeric(dat$frame)
# center/scale frame to improve numerical stability
dat$frame_sc <- scale(dat$frame, center = TRUE, scale = TRUE)

fit_poly <- procD.lm(
  coords ~ poly(frame_sc, 3) + clade + cluster +
    clade:poly(frame_sc, 3) + cluster:poly(frame_sc, 3) + clade:cluster,
  data = dat,
  SS.type = "III",
  iter = 999
)

anova(fit_poly, effect.type = "F")

# Try ordered approach
# 0) ensure frame is ordered 0..100 (or the set you actually have)
dat$frameF <- factor(dat$frame, levels = sort(unique(dat$frame)), ordered = FALSE)
dat$frameF

# 1) Fit Procrustes (RRPP) MANOVA with interactions that encode trajectories
#    Type III SS makes the clade/trait tests marginal (i.e., adjusted for the others)
fit22 <- procD.lm(
  coords ~ frameF + clade + cluster + clade:frameF + cluster:frameF + clade:cluster,
  data = dat,
  SS.type = "III",
  iter = 999
)

# 2) ANOVA table (permutation tests)
tab <- anova(fit22, effect.type = "F")
tab

# 3) Partial R^2 for each term (variance explained, adjusted for residuals)
SS <- tab$`SS`
residSS <- SS[which(rownames(tab) == "Residuals")]
partial_R2 <- SS[rownames(tab) != "Residuals"] / (SS[rownames(tab) != "Residuals"] + residSS)
data.frame(term = rownames(tab)[rownames(tab) != "Residuals"],
           partial_R2 = partial_R2,
           p_value = tab$`Pr(>F)`[rownames(tab) != "Residuals"])

summary(symDat_2D)

unique(symDat_2D$mass)
unique(symDat_2D$speed)
unique(symDat_2D$speed)
unique(symDat_2D$cluster)
unique(symDat_2D$clade)
unique(symDat_2D$frame)

fit1 <- lm.rrpp(coords ~ cluster + clade * frame, data = symDat_2D,
                      iter = 500, print.progress = TRUE)
summary(fit1)


# Run Trajectory Analysis ####

# 1) Q: How different are motion patterns in different hunting clusters?
#         a) Run trajectory analysis by hunting cluster on all data
# 2) Q: How sensitive are these results by differences in gait?
#         b) Run trajectory analysis on walking gait only
# 3) Q: How sensitive are these results to phylogeny?
#         c) Run trajectory analysis on clades
#         d) Run trajectory analysis on cluster with effect of clade removed
#         e) Run trajectory analysis on cluster using species level-analysis


## Symmetrical Data ####
### 1a) Run Hunting Cluster TAs ####


# Fit by Hunting
fit.symDat.cluster  <- lm.rrpp(coords ~ cluster * frame, 
                  data = symDat_2D, iter = 500, print.progress = TRUE)

summary(fit.symDat.cluster)

fit.symDat.cluster

TA.symDat.cluster <- trajectory.analysis(fit.symDat.cluster, groups = symDat_2D$cluster, 
                             traj.pts = symDat_2D$frame, print.progress = TRUE)
TA.symDat.cluster$pca


saveRDS(TA.symDat.cluster, "outputData/symDat_cluster_TA.RDS")

#### Make Stats & Plots ####
TA.symDat.cluster <- readRDS("outputData/symDat_cluster_TA.RDS")
summary(TA.symDat.cluster)

plot.symDat.cluster <- plot(TA.symDat.cluster, pch = 25, 
                           bg = col_fade, cex = 0.4,
                           xlim = c(-0.5, 0.4),
                           ylim = c(-0.3, 0.3))

plot.symDat.cluster

add.trajectories(plot.symDat.cluster, traj.pch = 21, 
                 traj.bg = col_full, 
                 start.bg = "black", end.bg = "white")

SD <- summary(TA.symDat.cluster, attribute = "SD")
MD <- summary(TA.symDat.cluster, attribute = "MD")
TC <- summary(TA.symDat.cluster, attribute = "TC", angle.type = "deg")

ta_Stats <- cbind(SD$summary.table, MD$summary.table, TC$summary.table)

# Print out the % Variation Explained by 
fitvals <- TA.symDat.cluster$fit$LM$fitted # fitted Procrustes coords
gp <- gm.prcomp(fitvals)
balls <- summary(gp) # % variation explained
balls$PC.summary

ta_Stats_list <- list("TA_Stats" = ta_Stats, "PC_contributions" = balls$PC.summary)
write.xlsx(ta_Stats_list, "./statsOutput/symDat_TA__byCluster_summaryStats.xlsx", rowNames = TRUE)

# 
TA.symDat.cluster$fit$LM$fitted

# Plot the PC Shape Extremes
plot_pc_shape_extremes_from_TA(symDat_2D, TA.symDat.cluster)

#### Pairwise: Which is Thylacine Closest To ? ####
TC$pairwise.tables

### 2b) Repeat Analysis ~ Cluster with ONLY walking ####

symDat_2D_walkOnly$frame
symDat_2D_walkOnly$cluster

# Fit by Hunting
fit.symDat_onlyWalk.cluster  <- lm.rrpp(coords ~ cluster * frame, 
                               data = symDat_2D_walkOnly, iter = 500, print.progress = TRUE)

summary(fit.symDat_onlyWalk.cluster)

TA.symDat_onlyWalk.cluster <- trajectory.analysis(fit.symDat_onlyWalk.cluster, groups = symDat_2D_walkOnly$cluster, 
                                         traj.pts = symDat_2D_walkOnly$frame, print.progress = TRUE)

saveRDS(TA.symDat_onlyWalk.cluster, "outputData/symDat_onlyWalk_cluster_TA.RDS")
#### Make Stats and Plots ####
TA.symDat_onlyWalk.cluster <- readRDS("outputData/symDat_onlyWalk_cluster_TA.RDS")

summary(TA.symDat_onlyWalk.cluster)

plot.symDat_onlyWalk.cluster <- plot(TA.symDat_onlyWalk.cluster, pch = 25, 
                            bg = col_fade_noWalk, cex = 0.4,
                            xlim = c(-0.5, 0.4),
                            ylim = c(-0.3, 0.3))

add.trajectories(plot.symDat_onlyWalk.cluster, traj.pch = 21, 
                 traj.bg = col_full_noWalk, 
                 start.bg = "black", end.bg = "white")

SD <- summary(TA.symDat_onlyWalk.cluster, attribute = "SD")
MD <- summary(TA.symDat_onlyWalk.cluster, attribute = "MD")
TC <- summary(TA.symDat_onlyWalk.cluster, attribute = "TC", angle.type = "deg")
ta_Stats <- cbind(SD$summary.table, MD$summary.table, TC$summary.table)

# Print out the % Variation Explained by 
fitvals <- TA.symDat_onlyWalk.cluster$fit$LM$fitted # fitted Procrustes coords
gp <- gm.prcomp(fitvals)
balls <- summary(gp) # % variation explained
balls$PC.summary

ta_Stats_list <- list("TA_Stats" = ta_Stats, "PC_contributions" = balls$PC.summary)
write.xlsx(ta_Stats_list, "./statsOutput/symDat_onlyWalk_TA__byCluster_summaryStats.xlsx", rowNames = TRUE)

#### Pairwise: Which is Thylacine Closest To ? ####
TC$pairwise.tables

### 3c) Run Trajectory Analysis on Clades ####

### Here use the data without Thylacine.
## Or better to set it to Dasyurid??...?

symDat_2D_noThy$frame
symDat_2D_noThy$cluster
unique(symDat_2D_noThy$clade)

# Fit by Hunting
fit.symDat.clade  <- lm.rrpp(coords ~ clade * frame, 
                               data = symDat_2D_noThy, iter = 500, print.progress = TRUE)

summary(fit.symDat.clade)


TA.symDat.clade <- trajectory.analysis(fit.symDat.clade, groups = symDat_2D_noThy$clade, 
                                         traj.pts = symDat_2D_noThy$frame, print.progress = TRUE)

saveRDS(TA.symDat.clade, "outputData/symDat_clade_TA.RDS")

#### Make Stats and Plots ####
TA.symDat.clade <- readRDS("outputData/symDat_clade_TA.RDS")

summary(TA.symDat.clade)

TA.symDat.clade$n.trajectories

plot.symDat.clade <- plot(TA.symDat.clade, pch = 25, 
                                     bg = col_fade_clade, cex = 0.4)

add.trajectories(plot.symDat.clade, traj.pch = 21, 
                 traj.bg = col_full_clade, 
                 start.bg = "black", end.bg = "white")

SD <- summary(TA.symDat.clade, attribute = "SD")
MD <- summary(TA.symDat.clade, attribute = "MD")
TC <- summary(TA.symDat.clade, attribute = "TC", angle.type = "deg")
ta_Stats <- cbind(SD$summary.table, MD$summary.table, TC$summary.table)

# Print out the % Variation Explained by 
fitvals <- TA.symDat.clade$fit$LM$fitted # fitted Procrustes coords
gp <- gm.prcomp(fitvals)
balls <- summary(gp) # % variation explained
balls$PC.summary

ta_Stats_list <- list("TA_Stats" = ta_Stats, "PC_contributions" = balls$PC.summary)
write.xlsx(ta_Stats_list, "./statsOutput/symDat_TA__byClade_summaryStats.xlsx", rowNames = TRUE)

### 3d) The ~ Cluster with the effect of Clade Removed ####

### Calculate Clade-Free Residuals of Shape ####

# Start with full data
symDat_2D

# Assume: dat has columns clade, cluster, frame
# coords is an n × q matrix (q = p*k), here q = 38
coords_mat <- as.matrix(symDat_2D$coords)
q <- ncol(coords_mat)

# Infer k (2D vs 3D). Here q %% 2 == 0 and q %% 3 != 0, so k = 2.
k <- if (q %% 3 == 0) 3 else 2
p <- q / k
stopifnot(p == floor(p))  # sanity

# Remove effect of Clade

fit_clade       <- lm.rrpp(coords ~ clade, data = symDat_2D, iter = 999)

# Get residuals (these are shapes with clade effects removed)
res_mat <- residuals(fit_clade, type = "response")


# Convert original and residuals to arrays
coords_arr <- arrayspecs(coords_mat, p = p, k = k)  # p × k × n
res_arr    <- arrayspecs(res_mat,   p = p, k = k)   # p × k × n

# Add back grand mean so we’re in shape space
grand_mean <- mshape(coords_arr)                     # p × k
coords_adj <- sweep(res_arr, MARGIN = 1:2, grand_mean, `+`)  # p × k × n

symDat_2D$coords

####
symDat_2D_noClade <- symDat_2D
symDat_2D_noClade$coords <- coords_adj

# set to 2D array again...
symDat_2D_noClade$coords <- two.d.array(symDat_2D_noClade$coords, sep = ".")
dim(symDat_2D_noClade$coords)

# New Fit
fit.symDat.cladeFree  <- lm.rrpp(coords ~ cluster * frame, data = symDat_2D_noClade, iter = 500)

summary(fit.symDat.cladeFree)

TA.symDat_cladeFree.cluster <- trajectory.analysis(fit.symDat.cladeFree, groups = symDat_2D_noClade$cluster, 
                                       traj.pts = symDat_2D_noClade$frame, print.progress = TRUE)

saveRDS(TA.symDat_cladeFree.cluster, "outputData/symDat_cladeFree_cluster_TA.RDS")

#### Make Stats and Plots ####

TA.symDat_cladeFree.cluster <- readRDS("outputData/symDat_cladeFree_cluster_TA.RDS")

summary(TA.symDat_cladeFree.cluster)

plot.symDat.clade <- plot(TA.symDat_cladeFree.cluster, pch = 25, 
                          bg = col_fade, cex = 0.4)

add.trajectories(plot.symDat.clade, traj.pch = 21, 
                 traj.bg = col_full, 
                 start.bg = "black", end.bg = "white")

SD <- summary(TA.symDat_cladeFree.cluster, attribute = "SD")
MD <- summary(TA.symDat_cladeFree.cluster, attribute = "MD")
TC <- summary(TA.symDat_cladeFree.cluster, attribute = "TC", angle.type = "deg")

ta_Stats <- cbind(SD$summary.table, MD$summary.table, TC$summary.table)

# Print out the % Variation Explained by 
fitvals <- TA.symDat_cladeFree.cluster$fit$LM$fitted # fitted Procrustes coords
gp <- gm.prcomp(fitvals)
balls <- summary(gp) # % variation explained
balls$PC.summary

ta_Stats_list <- list("TA_Stats" = ta_Stats, "PC_contributions" = balls$PC.summary)
write.xlsx(ta_Stats_list, "./statsOutput/symDat_cladeFree_TA__byCluster_summaryStats.xlsx", rowNames = TRUE)

#### Pairwise: Which is Thylacine Closest To ? ####
TC$pairwise.tables

### 3e) Using Species-Average Analysis ####

# Use the walk only data
summary(symDat_2D_walkOnly_noThy)

# Need to DE-Factor the carry_fields or the function will barf
symDat_2D_walkOnly_noThy$cluster <- as.character(symDat_2D_walkOnly_noThy$cluster)
symDat_2D_walkOnly_noThy$clade <- as.character(symDat_2D_walkOnly_noThy$clade)

# Do species level means #
symDat_walkOnly_spAvg <- species_average_gdf(symDat_2D_walkOnly_noThy, 
                                    align = "global",
                                    carry_fields = c("cluster", "clade"))


unique(symDat_walkOnly_spAvg$species)
unique(symDat_walkOnly_spAvg$cluster)
unique(symDat_walkOnly_spAvg$clade)


dim(symDat_walkOnly_spAvg$coords)

# Set $coords to 2D array 
symDat_2D_walkOnly_spAvg <- symDat_walkOnly_spAvg
symDat_2D_walkOnly_spAvg$coords <-  two.d.array(symDat_walkOnly_spAvg$coords)

# Re-Factor the metadata
symDat_2D_walkOnly_spAvg$frame <- as.factor(symDat_2D_walkOnly_spAvg$frame)
symDat_2D_walkOnly_spAvg$clade <- as.factor(symDat_2D_walkOnly_spAvg$clade)
symDat_2D_walkOnly_spAvg$cluster <- as.factor(symDat_2D_walkOnly_spAvg$cluster)

# Check the output
unique(symDat_2D_walkOnly_spAvg$frame)
unique(symDat_2D_walkOnly_spAvg$clade)
unique(symDat_2D_walkOnly_spAvg$cluster)

summary(symDat_2D_walkOnly_spAvg)
dim(symDat_2D_walkOnly_spAvg$coords)



# Fit by Hunting
fit.symDat_onlyWalk_spAvg.cluster  <- lm.rrpp(coords ~ cluster * frame, 
                                        data = symDat_2D_walkOnly_spAvg, iter = 500, print.progress = TRUE)

summary(fit.symDat_onlyWalk_spAvg.cluster)

TA.symDat_onlyWalk_spAvg.cluster <- trajectory.analysis(fit.symDat_onlyWalk_spAvg.cluster, groups = symDat_2D_walkOnly_spAvg$cluster, 
                                                  traj.pts = symDat_2D_walkOnly_spAvg$frame, print.progress = TRUE)

saveRDS(TA.symDat_onlyWalk_spAvg.cluster, "outputData/symDat_onlyWalk_spAvg_cluster_TA.RDS")

#### Make Stats and Plots ####
TA.symDat_onlyWalk_spAvg.cluster <- readRDS("outputData/symDat_onlyWalk_spAvg_cluster_TA.RDS")

summary(TA.symDat_onlyWalk_spAvg.cluster)

plot.symDat_onlyWalk_spAvg.cluster <- plot(TA.symDat_onlyWalk_spAvg.cluster, pch = 25, 
                                     bg = col_fade_noWalk, cex = 0.4,
                                     xlim = c(-0.5, 0.4),
                                     ylim = c(-0.3, 0.3))

add.trajectories(plot.symDat_onlyWalk_spAvg.cluster, traj.pch = 21, 
                 traj.bg = col_full_noWalk[1:4], 
                 start.bg = "black", end.bg = "white")

SD <- summary(TA.symDat_onlyWalk_spAvg.cluster, attribute = "SD")
MD <- summary(TA.symDat_onlyWalk_spAvg.cluster, attribute = "MD")
TC <- summary(TA.symDat_onlyWalk_spAvg.cluster, attribute = "TC", angle.type = "deg")

ta_Stats <- cbind(SD$summary.table, MD$summary.table, TC$summary.table)

# Print out the % Variation Explained by 
fitvals <- TA.symDat_onlyWalk_spAvg.cluster$fit$LM$fitted # fitted Procrustes coords
gp <- gm.prcomp(fitvals)
balls <- summary(gp) # % variation explained
balls$PC.summary

ta_Stats_list <- list("TA_Stats" = ta_Stats, "PC_contributions" = balls$PC.summary)

write.xlsx(ta_Stats, "./statsOutput/symDat_onlyWalk_spAvg_TA__byClade_summaryStats.xlsx", rowNames = TRUE)

#### Pairwise: Which is Thylacine Closest To ? ####
TC$pairwise.tables

## Assymetrical Data ####
### 1a) Run Hunting Cluster TAs ####

# Fit by Hunting
fit.symDat.cluster  <- lm.rrpp(coords ~ cluster * frame, 
                               data = symDat_2D, iter = 500, print.progress = TRUE)

summary(fit.symDat.cluster)

fit.symDat.cluster

TA.symDat.cluster <- trajectory.analysis(fit.symDat.cluster, groups = symDat_2D$cluster, 
                                         traj.pts = symDat_2D$frame, print.progress = TRUE)
TA.symDat.cluster$pca


saveRDS(TA.symDat.cluster, "outputData/symDat_cluster_TA.RDS")

TA.symDat.cluster <- readRDS("outputData/symDat_cluster_TA.RDS")
summary(TA.symDat.cluster)

plot.symDat.cluster <- plot(TA.symDat.cluster, pch = 25, 
                            bg = col_fade, cex = 0.4,
                            xlim = c(-0.5, 0.4),
                            ylim = c(-0.3, 0.3))

plot.symDat.cluster

add.trajectories(plot.symDat.cluster, traj.pch = 21, 
                 traj.bg = col_full, 
                 start.bg = "black", end.bg = "white")

SD <- summary(TA.symDat.cluster, attribute = "SD")
MD <- summary(TA.symDat.cluster, attribute = "MD")
TC <- summary(TA.symDat.cluster, attribute = "TC", angle.type = "deg")

ta_Stats <- cbind(SD$summary.table, MD$summary.table, TC$summary.table)

# Print out the % Variation Explained by 
fitvals <- TA.symDat.cluster$fit$LM$fitted # fitted Procrustes coords
gp <- gm.prcomp(fitvals)
balls <- summary(gp) # % variation explained
balls$PC.summary

ta_Stats_list <- list("TA_Stats" = ta_Stats, "PC_contributions" = balls$PC.summary)
write.xlsx(ta_Stats_list, "./statsOutput/symDat_TA_byClade_summaryStats.xlsx", rowNames = TRUE)

# 
TA.symDat.cluster$fit$LM$fitted

# Plot the PC Shape Extremes
plot_pc_shape_extremes_from_TA(symDat_2D, TA.symDat.cluster)