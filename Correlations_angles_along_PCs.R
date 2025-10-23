
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
#PCA_out <- plot(pca_sym)

#helper functions

# --- Oriented 2D angle between (A->B) and (C->D); base range 0..360
angle_between_lines <- function(x1,y1,x2,y2, x3,y3,x4,y4) {
  v1 <- c(x2 - x1, y2 - y1)
  v2 <- c(x4 - x3, y4 - y3)
  n1 <- sqrt(sum(v1^2)); n2 <- sqrt(sum(v2^2))
  if (n1 == 0 || n2 == 0) return(NA_real_)
  dot <- sum(v1 * v2)
  cr  <- v1[1]*v2[2] - v1[2]*v2[1]     # z-component cross
  ang <- atan2(cr, dot) * 180/pi       # (-180,180]
  (ang + 360) %% 360                   # -> [0,360)
}

# --- unwrap a degree series to be continuous (allow >360, <0, etc.)
unwrap_deg <- function(a) {
  out <- a
  for (i in 2:length(out)) {
    d <- out[i] - out[i-1]
    if (d > 180)  out[i:length(out)] <- out[i:length(out)] - 360
    if (d < -180) out[i:length(out)] <- out[i:length(out)] + 360
  }
  out
}

# ---  per-shape angles (ABSOLUTE angles, no folding)
angle_vec <- function(shape, PTS) {
  apply(PTS, 1, function(r)
    angle_between_lines(shape[r[1],1],shape[r[1],2],
                        shape[r[2],1],shape[r[2],2],
                        shape[r[3],1],shape[r[3],2],
                        shape[r[4],1],shape[r[4],2])
  ) |> as.numeric()
}

# --- angles across a sequence of shapes; unwrap each angle column
# shapes: p x 2 x n   (e.g., your symDat_3D$coords or a reconstructed path)
angle_vec_series <- function(shapes, PTS, unwrap = TRUE) {
  n  <- dim(shapes)[3]
  m  <- nrow(PTS)
  M  <- matrix(NA_real_, n, m)
  for (i in seq_len(n)) M[i,] <- angle_vec(shapes[,,i], PTS)
  if (unwrap) for (j in seq_len(m)) M[,j] <- unwrap_deg(M[,j])
  colnames(M) <- paste0("ang", seq_len(m))
  M
}


# --- Forelimb & hindlimb ratio for one 2D shape (p x 2)
limb_ratios <- function(shape) {
  # FL: 10-11-12-13-14 vs straight 10->14
  fl_poly <- seg_len(shape[10,],shape[11,]) + seg_len(shape[11,],shape[12,]) +
    seg_len(shape[12,],shape[13,]) + seg_len(shape[13,],shape[14,])
  fl_chord <- seg_len(shape[10,],shape[14,])
  # HL: 15-16-17-18-19 vs straight 15->19
  hl_poly <- seg_len(shape[15,],shape[16,]) + seg_len(shape[16,],shape[17,]) +
    seg_len(shape[17,],shape[18,]) + seg_len(shape[18,],shape[19,])
  hl_chord <- seg_len(shape[15,],shape[19,])
  c(FLRatio = fl_poly / fl_chord, HLRatio = hl_poly / hl_chord)
}


# --- Simple limb “polyline length” helper
seg_len <- function(a,b) sqrt(sum((a - b)^2))

#Bug fixed from your code: for the min ratios you accidentally used column 4 in the denominator. Above uses only cols 1–2 for min shapes (and 3–4 for max if you’re comparing two shapes side-by-side; see Approach B below).

PTS <- matrix(
  c(1, 3, 4, 5, 180,  # HeSp 1
    7, 8, 4, 5, 180,  # Spine1 2
    7, 8, 5, 6, 180,  # Spine2 3
    5, 7,15,16, 180,  # FeSpTheta 4
    5, 7,10,11,   0,  # HuSpTheta 5
    15,16,16,17, 180,  # KneeTheta 6
    16,17,17,18, 180,  # AnkleTheta 7
    10,11,11,12, 180,  # ElbowTheta 8
    11,12,12,13, 180), # WristTheta 9
  nrow = 9, byrow = TRUE
)


ang_names <- c("HeSp_1","Spine1_2","Spine2_3","FeSp_4","HuSp_5",
               "KneeTheta_6","AnkleTheta_7","ElbowTheta_8","WristTheta_9")

#Approach A — Empirical slope: angle ~ PC across all specimens
#This tells you, for each PC axis, how each angle changes per unit PC (with SE/p/R²). It uses aligned array directly


Y <- symDat_3D$coords            # p x 2 x n
pca <- pca_sym
n  <- dim(Y)[3]

# angles (0..360, unwrapped across frames)
A <- angle_vec_series(Y, PTS, unwrap = TRUE)
colnames(A) <- ang_names

# ratios per specimen (unchanged)
R <- matrix(NA_real_, nrow = n, ncol = 2); colnames(R) <- c("FLRatio","HLRatio")
for (i in seq_len(n)) R[i,] <- limb_ratios(Y[,,i])

head(A)

# For a chosen PC (e.g., PC1), fit angle ~ PC1 and ratios ~ PC1
# choose a PC:
pc <- 2
df <- data.frame(PC = pca$x[, pc], A, R)

get_fit_row <- function(y, pcname = "PC") {
  fit <- lm(y ~ PC, data = df)
  s   <- summary(fit)
  c(slope = coef(fit)[2],
    SE    = s$coefficients[2,2],
    t     = s$coefficients[2,3],
    p     = s$coefficients[2,4],
    R2    = s$r.squared)
}

res_angles <- t(apply(df[ , ang_names, drop = FALSE], 2, get_fit_row))
res_ratios <- t(apply(df[ , c("FLRatio","HLRatio")], 2, get_fit_row))

res_PC1 <- rbind(res_angles, res_ratios)
round(res_PC1, 4)


# ang_names <- c("HeSp_1","Spine1_2","Spine2_3","FeSp_4","HuSp_5",
#                "KneeTheta_6","AnkleTheta_7","ElbowTheta_8","WristTheta_9")


# assumes: ang_names, A (n × m angles), and pca_sym already exist
pc <- 2
df <- data.frame(PC = pca_sym$x[, pc], A)
measure_names <- colnames(A)  # or c(ang_names, "FLRatio","HLRatio") if you add ratios

op <- par(mfrow = c(3, 3), mar = c(3,3,2,1), mgp = c(1.8,0.6,0))
on.exit(par(op), add = TRUE)

for (nm in measure_names) {
  plot(df$PC, df[[nm]], pch = 20, cex = 0.6,
       xlab = sprintf("PC%d score", pc), ylab = nm)
  abline(lm(df[[nm]] ~ df$PC), lwd = 2, col='red')
}


#fix wrist angles if required. (only run once!)
#then rerun lines above. 
ind<-which(A[,9]>500)
A[ind,9]<-A[ind,9]-360



par(mfrow = c(1, 2))
    df <- data.frame(PC = pca$x[, pc], A, R)
    
    plot(df$PC, df$FLRatio, pch = 20, cex = 0.6,
         xlab = sprintf("PC%d score", pc), ylab = 'FLratio')
    abline(lm(df$FLRatio ~ df$PC), lwd = 2, col='red')
    
    plot(df$PC, df$HLRatio, pch = 20, cex = 0.6,
         xlab = sprintf("PC%d score", pc), ylab = 'HLRatio')
    abline(lm(df$HLRatio ~ df$PC), lwd = 2, col='red')


# install.packages(c("officer","flextable"))  # run once
library(officer)
library(flextable)

# --- prepare table from your model results ---
tab <- as.data.frame(res_PC1)                 # works for matrix or data.frame
tab$Measure <- rownames(tab)
# keep a sensible column order; adjust if your names differ
expected <- c("slope","SE","t","p","R2")
have <- intersect(expected, names(tab))
tab <- tab[, c("Measure", have), drop = FALSE]

# round numeric cols (except p which we format separately)
num_cols <- setdiff(intersect(names(tab), c("slope","SE","t","R2")), "p")
for (nm in num_cols) tab[[nm]] <- round(tab[[nm]], 4)

# pretty p-values + significance stars
fmt_p <- function(p) ifelse(p < 1e-4, "<0.0001", sprintf("%.4f", p))
stars <- function(p) ifelse(p < .001, "***",
                            ifelse(p < .01, "**",
                                   ifelse(p < .05, "*", "")))
if ("p" %in% names(tab)) {
  tab$p_sig <- paste0(fmt_p(tab$p), " ", stars(tab$p))
  tab$p <- NULL
  names(tab)[names(tab) == "p_sig"] <- "p"
}

# optional: nicer header labels
names(tab) <- sub("^SE$", "SE", names(tab))
names(tab) <- sub("^R2$", "R²", names(tab))

# --- build a Word table ---
ft <- flextable(tab)
ft <- autofit(ft)
ft <- align(ft, j = setdiff(names(tab), "Measure"), align = "center", part = "all")
ft <- bold(ft, part = "header")
ft <- valign(ft, valign = "center", part = "all")

# optional caption
cap <- "Angle/ratio change per unit PC1 score (linear model: y ~ PC1). p: two-sided; stars denote significance (* < .05, ** < .01, *** < .001)."

doc <- read_docx()
doc <- body_add_par(doc, "PC1 angle/ration slopes", style = "heading 2")
doc <- body_add_par(doc, cap, style = "Normal")
doc <- body_add_flextable(doc, ft)
print(doc, target = "PC1_angle_slopes_v2.docx")




cap <- "Angle/ratio change per unit PC2 score (linear model: y ~ PC2). p: two-sided; stars denote significance (* < .05, ** < .01, *** < .001)."

doc <- read_docx()
doc <- body_add_par(doc, "PC2 angle/ration slopes", style = "heading 2")
doc <- body_add_par(doc, cap, style = "Normal")
doc <- body_add_flextable(doc, ft)
print(doc, target = "PC2_angle_slopes_v2.docx")




