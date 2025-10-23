## ---- packages ----
pkgs <- c("dplyr","stringr","vegan","lme4","car","effectsize",
          "performance","flextable","officer","broom","broom.mixed")
invisible(lapply(pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)))
lapply(pkgs, library, character.only = TRUE)

## ---- helpers ----
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "<0.001", sprintf("%.3f", p)))
strip_rownames <- function(x) { x$Term <- rownames(x); rownames(x) <- NULL; x }

permanova_tables <- function(dat, gait_label) {
  # Full model (note: cluster may be aliased by clade/species)
  full <- adonis2(cbind(dat$PC1, dat$PC2) ~ clade + species + cluster + gait + mass,
                  data = dat, method = "euclidean", by = "terms") |> as.data.frame() |> strip_rownames()
  full$Model <- "Full"
  
  # Reduced model (cluster + gait + mass)
  red  <- adonis2(cbind(dat$PC1, dat$PC2) ~ cluster + mass + gait,
                  data = dat, method = "euclidean", by = "terms") |> as.data.frame() |> strip_rownames()
  red$Model <- "Reduced"
  
  bind_rows(full, red) |>
    filter(Term %in% c("clade","species","cluster","gait","mass")) |>
    transmute(
      Dataset = gait_label,
      Model,
      Term = recode(Term, clade="Clade", species="Species", cluster="Cluster", gait="Gait", mass="Mass"),
      Df, `Sum Sq` = round(SumOfSqs, 3),
      R2 = round(R2, 3),
      F  = round(F, 2),
      p  = fmt_p(`Pr(>F)`)
    )
}


permanova_tables <- function(dat, gait_label) {
  # Ensure factors
  dat <- dat |>
    dplyr::mutate(across(c(clade, species, cluster, gait), as.factor))
  
  # Response matrix (avoids cbind scoping)
  Y <- as.matrix(dat[, c("PC1","PC2")])
  
  # Full model
  full <- adonis2(Y ~ clade + species + cluster + gait + mass,
                  data = dat, method = "euclidean", by = "terms") |>
    as.data.frame()
  full$Term  <- rownames(full); rownames(full) <- NULL
  full$Model <- "Full"
  
  # Reduced model
  red <- adonis2(Y ~ cluster + mass + gait,
                 data = dat, method = "euclidean", by = "terms") |>
    as.data.frame()
  red$Term  <- rownames(red); rownames(red) <- NULL
  red$Model <- "Reduced"
  
  dplyr::bind_rows(full, red) |>
    dplyr::filter(Term %in% c("clade","species","cluster","gait","mass")) |>
    dplyr::transmute(
      Dataset = gait_label,
      Model,
      Term = dplyr::recode(Term, clade="Clade", species="Species",
                           cluster="Cluster", gait="Gait", mass="Mass"),
      Df,
      `Sum Sq` = round(SumOfSqs, 3),
      R2       = round(R2, 3),
      F        = round(F, 2),
      p        = ifelse(`Pr(>F)` < .001, "<0.001", sprintf("%.3f", `Pr(>F)`))
    )
}


axiswise_tables <- function(dat, gait_label) {
  m1 <- lmer(PC1 ~ cluster + mass + gait + (1|species), data = dat)
  m2 <- lmer(PC2 ~ cluster + mass + gait + (1|species), data = dat)
  
  # Type II (Wald) tests
  a1 <- car::Anova(m1, type = 2) |> as.data.frame() |> strip_rownames() |>
    filter(Term %in% c("cluster","gait","mass")) |>
    transmute(Dataset = gait_label, Axis = "PC1",
              Term = str_to_title(Term), `Chi²` = round(Chisq, 2), Df,
              p = fmt_p(`Pr(>Chisq)`))
  a2 <- car::Anova(m2, type = 2) |> as.data.frame() |> strip_rownames() |>
    filter(Term %in% c("cluster","gait","mass")) |>
    transmute(Dataset = gait_label, Axis = "PC2",
              Term = str_to_title(Term), `Chi²` = round(Chisq, 2), Df,
              p = fmt_p(`Pr(>Chisq)`))
  
  # Partial η² (Type III SS by default for LMMs; commonly reported)
  e1 <- effectsize::eta_squared(m1, partial = TRUE) |> as.data.frame() |>
    transmute(Term = str_to_title(Parameter), `η²p` = round(Eta2_partial, 2))
  e2 <- effectsize::eta_squared(m2, partial = TRUE) |> as.data.frame() |>
    transmute(Term = str_to_title(Parameter), `η²p` = round(Eta2_partial, 2))
  
  tab1 <- left_join(a1, e1, by = "Term")
  tab2 <- left_join(a2, e2, by = "Term")
  ax_main <- bind_rows(tab1, tab2)
  
  # R² (reported separately)
  r21 <- performance::r2_nakagawa(m1); r22 <- performance::r2_nakagawa(m2)
  ax_r2 <- tibble::tibble(
    Dataset = gait_label,
    Axis    = c("PC1","PC2"),
    `Marginal R²`    = round(c(r21$R2_marginal, r22$R2_marginal), 3),
    `Conditional R²` = round(c(r21$R2_conditional, r22$R2_conditional), 3)
  )
  
  list(main = ax_main, r2 = ax_r2)
}

## ---- load your data (edit paths as needed) ----
# Symmetrical (remove thylacine as you’ve been doing)
sym_dat <- read.csv("bigdat_Symmetrical_v2.csv") |>
  dplyr::filter(species != "Thylacinus cynocephalus") |>
  mutate(across(c(clade, species, cluster, gait), as.factor))

# Asymmetrical
asym_dat <- read.csv("bigdat_Asymmetrical_v2.csv") |>
  mutate(across(c(clade, species, cluster, gait), as.factor))

## ---- build tables ----
tab_perm <- bind_rows(
  permanova_tables(sym_dat,  "Symmetrical"),
  permanova_tables(asym_dat, "Asymmetrical")
)

ax_sym  <- axiswise_tables(sym_dat,  "Symmetrical")
ax_asym <- axiswise_tables(asym_dat, "Asymmetrical")

tab_axis   <- bind_rows(ax_sym$main, ax_asym$main)
tab_axis_r2<- bind_rows(ax_sym$r2,   ax_asym$r2)

## ---- Word tables (flextable + officer) ----
ft_perm <- flextable(tab_perm) |>
  autofit() |>
  flextable::theme_vanilla()

ft_axis <- flextable(tab_axis) |>
  autofit() |>
  flextable::theme_vanilla()

ft_r2 <- flextable(tab_axis_r2) |>
  autofit() |>
  flextable::theme_vanilla()

doc <- read_docx() |>
  body_add_par("Table S1. PERMANOVA results (by terms)", style = "heading 2") |>
  body_add_flextable(ft_perm) |>
  body_add_par("Notes: Full model includes clade + species + cluster + gait + mass; Reduced model omits clade and species to assess cluster independent of phylogenetic nesting.", style = "Normal") |>
  body_add_par(" ", style = "Normal") |>
  body_add_par("Table S2. Axis-wise mixed models (Type II Wald tests and partial η²)", style = "heading 2") |>
  body_add_flextable(ft_axis) |>
  body_add_par("Notes: Models are PC ~ cluster + gait + mass + (1|species). η²p from effectsize (Type III SS); tests are Type II Wald χ².", style = "Normal") |>
  body_add_par(" ", style = "Normal") |>
  body_add_par("Table S3. Axis-wise variance explained (Nakagawa’s R²)", style = "heading 2") |>
  body_add_flextable(ft_r2)

print(doc, target = "Supplementary_Tables.docx")
