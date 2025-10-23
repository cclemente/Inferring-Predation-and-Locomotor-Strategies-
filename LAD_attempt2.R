
library(class)

## Assume you still have:
## train: non-thylacine rows with cluster != 6, columns PC1, PC2, cluster
## test_X: the two thylacine rows with PC1, PC2

bigdat <- read.csv('bigdat_Symmetrical_v2.csv')

target <- "Thylacinus cynocephalus"

# Train on known clusters (drop thylacine + any unknown '6')
train <- subset(bigdat, species != target & !is.na(PC1) & !is.na(PC2) & cluster != 6)
train$cluster <- droplevels(factor(train$cluster))

# Test (thylacine) – only PCs needed
test_X <- subset(bigdat, species == target, select = c(PC1, PC2))



# 1) Thylacine centroid (species-level mean in PC space)
thyl_mu <- data.frame(PC1 = mean(test_X$PC1), PC2 = mean(test_X$PC2))

# 2) Nearest centroid (Euclidean)
cent <- aggregate(cbind(PC1,PC2) ~ cluster, data = train, mean)
MU   <- as.matrix(cent[, c("PC1","PC2")])
x    <- as.numeric(thyl_mu)
euclid_d   <- rowSums((MU - matrix(x, nrow(MU), 2, byrow=TRUE))^2)
euclid_pick <- cent$cluster[which.min(euclid_d)]

# 3) Mahalanobis (uses pooled covariance for stability)
covmat <- cov(train[, c("PC1","PC2")])
if (det(covmat) <= .Machine$double.eps) {
  covmat <- diag(apply(train[, c("PC1","PC2")], 2, var, na.rm=TRUE))
}
maha_d     <- apply(MU, 1, function(m) mahalanobis(thyl_mu, center = m, cov = covmat))
maha_pick  <- cent$cluster[which.min(maha_d)]

# 4) kNN on the centroid (nonparametric check)
k <- max(1, 2*floor(sqrt(nrow(train))/2)+1)  # odd k ~ sqrt(n)
knn_cls_mu <- knn(train[, c("PC1","PC2")], thyl_mu, cl = droplevels(train$cluster), k = k, prob = TRUE)
knn_p_mu   <- attr(knn_cls_mu, "prob")

# 5) Print a concise summary
cat("\nThylacine centroid (PC1, PC2):", round(thyl_mu$PC1, 3), ",", round(thyl_mu$PC2, 3), "\n")
cat("Nearest-centroid (Euclidean):", as.character(euclid_pick), "  distances:", round(euclid_d, 4), "\n")
cat("Nearest-centroid (Mahalanobis):", as.character(maha_pick), "  distances:", round(maha_d, 4), "\n")
cat("kNN on centroid: class", as.character(knn_cls_mu), " (p =", round(knn_p_mu, 3), ", k =", k, ")\n")

# Optional: tie-aware call if distances are nearly equal (within 10% of min)
tie_thresh <- 0.10
close_euclid <- cent$cluster[(euclid_d - min(euclid_d))/min(euclid_d) <= tie_thresh]
close_maha   <- cent$cluster[(maha_d   - min(maha_d))/min(maha_d)   <= tie_thresh]
cat("Tie-aware (Euclid) close clusters:", paste(close_euclid, collapse="/"), "\n")
cat("Tie-aware (Mahalanobis) close clusters:", paste(close_maha, collapse="/"), "\n")





set.seed(42)
B  <- 1000
lev <- levels(train$cluster)
eu <- setNames(integer(length(lev)), lev)
mh <- setNames(integer(length(lev)), lev)

for (b in 1:B) {
  # Stratified bootstrap: resample within each cluster to the same size
  boot <- do.call(rbind, lapply(split(train, train$cluster), function(df) {
    idx <- sample.int(nrow(df), size = nrow(df), replace = TRUE)
    df[idx, , drop = FALSE]
  }))
  boot$cluster <- droplevels(boot$cluster)
  
  # Recompute cluster centroids
  cent_b <- aggregate(cbind(PC1, PC2) ~ cluster, data = boot, mean)
  MUb    <- as.matrix(cent_b[, c("PC1","PC2")])
  x      <- as.numeric(thyl_mu)
  
  ## Euclidean winner
  ed <- rowSums((MUb - matrix(x, nrow(MUb), 2, byrow = TRUE))^2)
  eu[as.character(cent_b$cluster[which.min(ed)])] <-
    eu[as.character(cent_b$cluster[which.min(ed)])] + 1
  
  ## Mahalanobis winner (pooled covariance; fall back if singular)
  S <- try(cov(boot[, c("PC1","PC2")]), silent = TRUE)
  if (inherits(S, "try-error") || !is.finite(det(S)) || det(S) <= .Machine$double.eps) {
    S <- diag(apply(boot[, c("PC1","PC2")], 2, var, na.rm = TRUE))
  }
  md <- apply(MUb, 1, function(m) mahalanobis(thyl_mu, center = m, cov = S))
  mh[as.character(cent_b$cluster[which.min(md)])] <-
    mh[as.character(cent_b$cluster[which.min(md)])] + 1
}

# Proportions of wins across bootstraps
prop_eu <- round(eu / B, 3)
prop_mh <- round(mh / B, 3)

prop_eu
prop_mh



#figures








