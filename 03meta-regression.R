library(tidyverse)
library(meta)
library(metafor)
library(readxl)

# Function for performing meta-regression for a given effect size
perform_meta_regression <- function(df, effect, se_col, group_col, meta_covs, dis_levels) {
  results <- data.frame(
    contrast = character(),
    covariate = character(),
    estimate = numeric(),
    ci.lb = numeric(),
    ci.ub = numeric(),
    p.value = numeric(),
    tau = numeric(),
    I2 = numeric(),
    H2 = numeric(),
    R2 = numeric(),
    stringsAsFactors = FALSE
  )

  for (dis in dis_levels) {
    for (cov in meta_covs) {
      model <- tryCatch({
        rma(
          yi = df[[effect]],
          sei = df[[se_col]],
          data = df[df[[group_col]] == dis, ],
          method = "ML",
          mods = as.formula(paste0("~", cov)),
          test = "knha"
        )
      }, error = function(e) NULL)

      if (!is.null(model)) {
        results <- rbind(results, data.frame(
          contrast = dis,
          covariate = cov,
          estimate = coef(model)[2],
          ci.lb = model[["ci.lb"]][2],
          ci.ub = model[["ci.ub"]][2],
          p.value = coef(summary(model))[2, "pval"],
          tau = sqrt(model$tau2),
          I2 = model$I2,
          H2 = model$H2,
          R2 = model$R2,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(results)
}

###--- Meta-regression for Risk Difference (RD) ---###

df.meta2 <- df.meta |> 
  merge(read_xlsx("./data/meta_data.xlsx"), by = "cohort")

# Random effects meta-analysis (no overall effect, subgroup by contrast)
meta_rd2 <- metagen(
  TE = RD,
  seTE = se,
  studlab = cohort,
  data = df.meta2,
  sm = "SMD",
  random = TRUE,
  common = FALSE,
  overall = FALSE,
  subgroup = contrast
)

# Define unique contrasts and meta-covariates
meta.dis <- unique(df.meta2$contrast)
meta.cov <- names(df.meta2)[7:16]

# Perform meta-regression for RD
results_rd <- perform_meta_regression(
  df = df.meta2,
  effect = "RD",
  se_col = "se",
  group_col = "contrast",
  meta_covs = meta.cov,
  dis_levels = meta.dis
)

write.csv(results_rd, file = "./results/meta_regression_RD.csv", row.names = FALSE)

###--- Meta-regression for Odds Ratio (OR) ---###

df.meta2 <- df.or |> 
  merge(read_xlsx("./data/meta_data.xlsx"), by = "cohort") |> 
  mutate(se = (CI_Upper - CI_Lower) / (2 * qnorm(0.975)))

# Perform meta-regression for OR
results_or <- perform_meta_regression(
  df = df.meta2,
  effect = "OR",
  se_col = "se",
  group_col = "Variable",
  meta_covs = meta.cov,
  dis_levels = meta.dis
)

write.csv(results_or, file = "./results/meta_regression_OR.csv", row.names = FALSE)


