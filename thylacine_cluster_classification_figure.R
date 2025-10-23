
install.packages('patchwork')

# ----- Packages -----
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(scales)
library(cowplot)  # or cowplot if you prefer

# ----- CONFIG -----
target_species <- "Thylacinus cynocephalus"

# If you need to restrict to symmetrical gaits, do it here, e.g.:
# bigdat <- subset(bigdat, gait_group == "sym")  # <-- adjust to your column

# Map numeric clusters to names
cluster_labels <- c(
  "Anteater",                  # 1
  "Opportunistic Grappler",    # 2
  "Large Grappler",            # 3
  "Opportunistic Pouncer",     # 4
  "Social Hunters"             # 5
)

# Map numeric clusters to names
cluster_labels <- c(
  "Anteaters",                  # 1
  "OppGrap",    # 2
  "LrgGrapp",            # 3
  "OppPounce",     # 4
  "SocHunt"             # 5
)


# ----- TRAIN/TEST SPLIT -----
train <- bigdat %>%
  filter(species != target_species,
         !is.na(PC1), !is.na(PC2),
         cluster %in% 1:5) %>%                 # drop '6' unknown
  mutate(cluster = factor(cluster, levels = 1:5, labels = cluster_labels))

test  <- bigdat %>%
  filter(species == target_species,
         !is.na(PC1), !is.na(PC2)) %>%
  select(PC1, PC2)

stopifnot(nrow(test) >= 1, nrow(train) >= 5)

# Thylacine centroid
thyl_mu <- data.frame(PC1 = mean(test$PC1), PC2 = mean(test$PC2))

# Cluster centroids (training set)
centroids <- train %>%
  group_by(cluster) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

# Pooled covariance for Mahalanobis (fallback to diagonal if singular)
S <- try(cov(train[, c("PC1","PC2")]), silent = TRUE)
if (inherits(S, "try-error") || !is.finite(det(S)) || det(S) <= .Machine$double.eps) {
  S <- diag(apply(train[, c("PC1","PC2")], 2, var, na.rm = TRUE))
}

# Euclidean vs Mahalanobis nearest cluster to thylacine centroid
MU <- as.matrix(centroids[, c("PC1","PC2")])
x  <- as.numeric(thyl_mu)

euclid_d <- rowSums((MU - matrix(x, nrow(MU), 2, byrow = TRUE))^2)
euclid_pick <- centroids$cluster[which.min(euclid_d)]

maha_d <- apply(MU, 1, function(m) mahalanobis(thyl_mu, center = m, cov = S))
maha_pick <- centroids$cluster[which.min(maha_d)]

# ----- PANEL A: Morphospace scatter with ellipses, centroids, and thylacine -----
# Only draw ellipses for clusters with >= 3 points (stat_ellipse needs enough data)
library(dplyr)
# clusters with enough points for ellipses
ellipse_ok <- train %>% count(cluster) %>% filter(n >= 3)

pA <- ggplot(train, aes(PC1, PC2, color = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_point(data = centroids, aes(PC1, PC2, color = cluster),
             inherit.aes = FALSE, shape = 21, fill = NA, size = 4, stroke = 1.1, show.legend = FALSE) +
  (if (nrow(ellipse_ok) > 0)
    stat_ellipse(data = dplyr::semi_join(train, ellipse_ok, by = "cluster"),
                 type = "norm", level = 0.95, linewidth = 0.6, alpha = 0.4)
   else NULL) +
  geom_point(data = test, aes(PC1, PC2), inherit.aes = FALSE,
             shape = 17, size = 3, color = "black") +
  geom_point(data = thyl_mu, aes(PC1, PC2), inherit.aes = FALSE,
             shape = 8, size = 3.6, color = "black") +
  geom_segment(data = centroids %>% filter(cluster == euclid_pick),
               aes(x = thyl_mu$PC1, y = thyl_mu$PC2, xend = PC1, yend = PC2),
               inherit.aes = FALSE, linetype = "dashed", linewidth = 0.7) +
  geom_segment(data = centroids %>% filter(cluster == maha_pick),
               aes(x = thyl_mu$PC1, y = thyl_mu$PC2, xend = PC1, yend = PC2),
               inherit.aes = FALSE, linetype = "solid", linewidth = 0.8) +
  ggrepel::geom_text_repel(
    data = centroids,
    aes(x = PC1, y = PC2, label = cluster),   # <-- add x,y here
    inherit.aes = FALSE, size = 3, box.padding = 0.3, seed = 1
  ) +
  scale_color_brewer(palette = "Dark2", name = "Cluster") +
  labs(
    title = "A. Symmetrical-gait morphospace (PC1–PC2)",
    subtitle = paste0("Dashed = Euclidean nearest (", euclid_pick,
                      "); Solid = Mahalanobis nearest (", maha_pick, ")"),
    x = "PC1", y = "PC2"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")





# ----- PANEL B: Bootstrap assignment stability -----
set.seed(42)
B <- 1000
lev <- levels(train$cluster)
eu <- setNames(integer(length(lev)), lev)
mh <- setNames(integer(length(lev)), lev)

for (b in 1:B) {
  boot <- do.call(rbind, lapply(split(train, train$cluster), function(df) {
    idx <- sample.int(nrow(df), size = nrow(df), replace = TRUE)
    df[idx, , drop = FALSE]
  }))
  boot_cent <- boot %>% group_by(cluster) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
  MUb <- as.matrix(boot_cent[, c("PC1","PC2")])
  # Euclidean
  ed <- rowSums((MUb - matrix(x, nrow(MUb), 2, byrow = TRUE))^2)
  eu[as.character(boot_cent$cluster[which.min(ed)])] <- eu[as.character(boot_cent$cluster[which.min(ed)])] + 1
  # Mahalanobis (pooled, robustified)
  S_b <- try(cov(boot[, c("PC1","PC2")]), silent = TRUE)
  if (inherits(S_b, "try-error") || !is.finite(det(S_b)) || det(S_b) <= .Machine$double.eps) {
    S_b <- diag(apply(boot[, c("PC1","PC2")], 2, var, na.rm = TRUE))
  }
  md <- apply(MUb, 1, function(m) mahalanobis(thyl_mu, center = m, cov = S_b))
  mh[as.character(boot_cent$cluster[which.min(md)])] <- mh[as.character(boot_cent$cluster[which.min(md)])] + 1
}
prop_eu <- eu / B
prop_mh <- mh / B

boot_df <- tibble(
  cluster = factor(names(prop_mh), levels = lev),
  Mahalanobis = as.numeric(prop_mh),
  Euclidean   = as.numeric(prop_eu)
) |>
  pivot_longer(cols = c(Mahalanobis, Euclidean), names_to = "Method", values_to = "Prop") |>
  arrange(desc(Method), desc(Prop))

pB <- ggplot(boot_df, aes(x = reorder(cluster, Prop), y = Prop, fill = Method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.7)) +
  scale_fill_manual(values = c(Mahalanobis = "#4E79A7", Euclidean = "#A0CBE8")) +
  labs(title = "B. Bootstrap assignment stability (B = 1,000)",
       x = NULL, y = "Win proportion",
       subtitle = "Mahalanobis vs Euclidean nearest-centroid of thylacine centroid") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# ----- PANEL C: Mahalanobis percentile (centrality) -----
# MD distribution of training points to *their own* cluster centroid
md_train <- split(train, train$cluster) |>
  lapply(function(df) {
    mu <- colMeans(df[, c("PC1","PC2")], na.rm = TRUE)
    mahalanobis(df[, c("PC1","PC2")], center = mu, cov = S)
  })

# Thylacine centroid MD to each cluster's centroid
md_thyl <- sapply(levels(train$cluster), function(k) {
  mu_k <- as.numeric(centroids[centroids$cluster == k, c("PC1","PC2")])
  mahalanobis(thyl_mu, center = mu_k, cov = S)
})

# Percentile of thylacine MD within each cluster’s MD distribution
perc_df <- tibble(
  cluster = factor(names(md_thyl), levels = lev),
  percentile = mapply(function(k, val) mean(md_train[[k]] <= as.numeric(val), na.rm = TRUE),
                      names(md_thyl), md_thyl)
)

pC <- ggplot(perc_df, aes(x = cluster, y = percentile)) +
  geom_hline(yintercept = 0.5, linetype = "dotted") +
  geom_point(size = 3, color = "#4E79A7") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "C. Centrality of thylacine centroid within each cluster",
       x = NULL, y = "Mahalanobis percentile (lower = more central)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))




# ----- Assemble & Save -----



# 1) Pull a single legend (from pA) and hide legends on the panels
leg <- get_legend(pA + theme(legend.position = "bottom") +
                    guides(color = guide_legend(nrow = 1, byrow = TRUE)))

pA_nl <- pA + theme(legend.position = "none")
pB_nl <- pB + theme(legend.position = "none")
pC_nl <- pC + theme(legend.position = "none")

# 2) Bottom row (B | C)
rowBC <- plot_grid(
  pB_nl, pC_nl,
  ncol = 2, rel_widths = c(1, 1),
  labels = c("B", "C"), label_size = 12,
  align = "hv", axis = "tblr"
)

# 3) Stack A over (B|C)
fig_core <- plot_grid(
  pA_nl, rowBC,
  ncol = 1, rel_heights = c(2, 1),
  labels = c("A", ""), label_size = 12,
  align = "v"
)

# 4) Add the shared legend at the bottom
fig <- plot_grid(
  fig_core, leg,
  ncol = 1, rel_heights = c(1, 0.10)
)

fig




###


## ====== Build a clean trajectory table (CLUSTER, not clade) ======
# If your cluster labels are numeric 1..5, map to names here:
cluster_labels <- c(
  `1` = "Anteater",
  `2` = "Opportunistic Grappler",
  `3` = "Large Grappler",
  `4` = "Opportunistic Pouncer",
  `5` = "Social Hunters"
)

traj <- merge(
  setNames(aggx, c("frame","cluster","PC1")),
  setNames(aggy, c("frame","cluster","PC2")),
  by = c("frame","cluster"), all = FALSE
)
traj <- traj[complete.cases(traj), ]

# Harmonize cluster names
traj$cluster <- as.character(traj$cluster)
traj$cluster_name <- ifelse(traj$cluster %in% names(cluster_labels),
                            cluster_labels[traj$cluster], traj$cluster)

# order by frame within cluster
traj <- traj[order(traj$cluster_name, traj$frame), ]

## ====== Colors ======
cls  <- sort(unique(traj$cluster_name))
cols <- setNames(hcl.colors(length(cls), "Dark 3"), cls)
cols_pt <- lapply(cols, adjustcolor, alpha.f = 0.55)  # transparent points

## ====== Optional: thylacine observations & centroid ======
# Try to get the two thylacine points from previous 'test' data; fallback to clade_gdf
if (exists("test") && all(c("PC1","PC2") %in% names(test))) {
  thyl_pts <- as.matrix(test[, c("PC1","PC2")])
} else if (exists("clade_gdf") && exists("sng")) {
  thyl_pts <- cbind(PC1 = clade_gdf$Xav[sng], PC2 = clade_gdf$Yav[sng])
} else {
  thyl_pts <- NULL
}
thyl_mu <- if (!is.null(thyl_pts)) c(PC1 = mean(thyl_pts[,1], na.rm=TRUE),
                                     PC2 = mean(thyl_pts[,2], na.rm=TRUE)) else NULL

## ====== Cluster centroids for distance lines ======
# Prefer centroids from 'train' (non-thylacine dataset); else use traj means
if (exists("train") && all(c("PC1","PC2","cluster") %in% names(train))) {
  trn <- train
  # map train cluster to names if numeric
  trn$cluster <- as.character(trn$cluster)
  trn$cluster_name <- ifelse(trn$cluster %in% names(cluster_labels),
                             cluster_labels[trn$cluster], trn$cluster)
  centroids <- aggregate(cbind(PC1,PC2) ~ cluster_name, data = trn, mean)
  names(centroids)[1] <- "cluster"
} else {
  centroids <- aggregate(cbind(PC1,PC2) ~ cluster_name, data = traj, mean)
  names(centroids)[1] <- "cluster"
}
centroids <- centroids[centroids$cluster %in% cls, , drop=FALSE]

## Covariance for Mahalanobis (pooled); fallback to diagonal if singular
S <- if (exists("train")) try(cov(train[, c("PC1","PC2")]), silent=TRUE) else try(cov(traj[, c("PC1","PC2")]), silent=TRUE)
if (inherits(S, "try-error") || !is.finite(det(S)) || det(S) <= .Machine$double.eps) {
  base_mat <- if (exists("train")) train[, c("PC1","PC2")] else traj[, c("PC1","PC2")]
  S <- diag(apply(base_mat, 2, function(v) var(v, na.rm=TRUE)))
}

## ====== Compute nearest clusters from thylacine centroid (if available) ======
euclid_pick <- maha_pick <- NULL
if (!is.null(thyl_mu)) {
  MU <- as.matrix(centroids[, c("PC1","PC2")])
  x  <- as.numeric(thyl_mu)
  euclid_d <- rowSums((MU - matrix(x, nrow(MU), 2, byrow=TRUE))^2)
  maha_d   <- apply(MU, 1, function(m) mahalanobis(thyl_mu, center = m, cov = S))
  euclid_pick <- centroids$cluster[which.min(euclid_d)]
  maha_pick   <- centroids$cluster[which.min(maha_d)]
}

## ====== Plot canvas (range includes traj, centroids, and thylacine) ======
x_rng <- range(c(traj$PC1, centroids$PC1, if (!is.null(thyl_mu)) thyl_mu["PC1"]), na.rm=TRUE)
y_rng <- range(c(traj$PC2, centroids$PC2, if (!is.null(thyl_mu)) thyl_mu["PC2"]), na.rm=TRUE)

plot(x_rng, y_rng, type="n", xlab="PC1", ylab="PC2",
     main = "Average stride path in PC space by cluster")

## ====== Draw per-frame points + paths per cluster ======
for (cl in cls) {
  sub <- traj[ traj$cluster_name == cl, ]
  o   <- order(sub$frame)
  # per-frame mean points (semi-transparent)
  points(sub$PC1[o], sub$PC2[o], pch = 16, cex = 0.6, col = cols_pt[[cl]])
  # path over the stride
  lines(sub$PC1[o], sub$PC2[o], col = cols[cl], lwd = 2)
}

## ====== Draw cluster centroids ======
for (cl in centroids$cluster) {
  with(centroids[centroids$cluster==cl,],
       points(PC1, PC2, pch = 21, cex = 1.6, col = cols[cl], bg = NA, lwd = 1.2))
}

## ====== Thylacine points + centroid + distance lines ======
if (!is.null(thyl_pts)) {
  points(thyl_pts[,1], thyl_pts[,2], pch = 17, cex = 1.1, col = "black")
}
if (!is.null(thyl_mu)) {
  points(thyl_mu["PC1"], thyl_mu["PC2"], pch = 8, cex = 1.3, col = "black", lwd = 1.2)
  # lines to nearest centroids
  if (!is.null(euclid_pick)) {
    with(centroids[centroids$cluster==euclid_pick,],
         segments(thyl_mu["PC1"], thyl_mu["PC2"], PC1, PC2, lty = 2, lwd = 1.1))  # dashed (Euclid)
  }
  if (!is.null(maha_pick)) {
    with(centroids[centroids$cluster==maha_pick,],
         segments(thyl_mu["PC1"], thyl_mu["PC2"], PC1, PC2, lty = 1, lwd = 1.6))  # solid (Mahalanobis)
  }
}

## ====== Non-overlapping labels at mid-frame along each path ======
labpos <- do.call(rbind, lapply(cls, function(cl){
  sub <- traj[ traj$cluster_name == cl, ]
  sub <- sub[order(sub$frame), ]
  mid <- sub[ which.min(abs(sub$frame - median(sub$frame, na.rm=TRUE))) , ]
  data.frame(cluster = cl, x = mid$PC1, y = mid$PC2)
}))

if (requireNamespace("maptools", quietly = TRUE)) {
  maptools::pointLabel(labpos$x, labpos$y, labels = labpos$cluster,
                       cex = 0.95, col = cols[labpos$cluster], doPlot = TRUE)
} else {
  text(labpos$x, labpos$y, labels = labpos$cluster,
       cex = 0.95, col = cols[labpos$cluster], pos = 3, xpd = NA)
}

legend("topright", legend = cls, col = cols, lwd = 2, bty = "n", cex = 0.9, title = "Cluster",
       inset = 0.01)



## --- Overlay ALL stride points (PC1/PC2), colored by cluster ---

# 1) Factor = cluster IDs for each stride
Factor <- clade_gdf$cluster  # numeric 1..5 (6 = unknown), or a factor

# 2) Which rows to plot? If you use an index 'sng', keep it; else plot all
idx <- if (exists("sng")) {
  sng
} else {
  rep(TRUE, length(Factor))
}

# 3) Extract stride-level PC scores
x_str <- clade_gdf$Xav[idx]
y_str <- clade_gdf$Yav[idx]
fac_str <- as.character(Factor[idx])

# 4) Map stride clusters to the same names/colors used for the paths
#    (cluster_labels and cols were defined earlier)
name_str <- ifelse(fac_str %in% names(cluster_labels),
                   cluster_labels[fac_str],
                   fac_str)

# 5) Drop NAs / bad rows
ok <- is.finite(x_str) & is.finite(y_str) & !is.na(name_str) & name_str %in% names(cols)
x_str <- x_str[ok]; y_str <- y_str[ok]; name_str <- name_str[ok]

# 6) Build a per-point color vector (semi-transparent so the cloud doesn’t overwhelm)
col_str <- adjustcolor(unname(cols[name_str]), alpha.f = 0.5)

# 7) Plot the stride cloud
points(x_str, y_str, pch = 16, cex = 1, col = col_str)

# (Optional) if you want slightly larger, less transparent points:
# points(x_str, y_str, pch = 16, cex = 0.9, col = adjustcolor(unname(cols[name_str]), alpha.f = 0.4))


