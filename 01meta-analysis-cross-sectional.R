if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, haven, psych, lme4, poLCA, forestplot, table1,
  gridExtra, ggpubr, autoReg, readxl, marginaleffects, rrtable,
  meta, metafor
)

# Load datasets
load("./data/charls.RData")
load("./data/elsa.RData")
load("./data/hrs.RData")
load("./data/mhas.RData")
load("./data/share.RData")

# Function to extract odds ratios from logistic regression models
get_or <- function(data) {
  fit.log <- glm(
    adl1 ~ age + sex + marr + edu + exercise + smoke + drink + income.c + bmi.c + detailed_category,
    family = "binomial",
    data = data
  )
  
  summary_fit <- summary(fit.log)
  last_13_coef <- tail(summary_fit$coefficients, 13)
  
  exp_coef <- exp(last_13_coef[, "Estimate"])
  exp_ci_lower <- exp(last_13_coef[, "Estimate"] - 1.96 * last_13_coef[, "Std. Error"])
  exp_ci_upper <- exp(last_13_coef[, "Estimate"] + 1.96 * last_13_coef[, "Std. Error"])
  p_values <- last_13_coef[, "Pr(>|z|)"]
  
  results <- data.frame(
    Variable = rownames(last_13_coef),
    OR = exp_coef,
    CI_Lower = exp_ci_lower,
    CI_Upper = exp_ci_upper,
    P_Value = p_values
  ) |>
    mutate(Variable = gsub("detailed_category", "", Variable))
  
  test <- data |>
    dplyr::select(detailed_category, adl1) |>
    table() |>
    as.data.frame() |>
    pivot_wider(names_from = adl1, values_from = Freq) |>
    rename("Variable" = detailed_category, c.n = `0`, e.n = `1`) |>
    mutate(n1 = rowSums(across(c(c.n, e.n))))
  
  test2 <- test[2:14, ] |>
    mutate(
      e.n2 = as.numeric(test[1, 3]),
      n2 = as.numeric(test[1, 4])
    )
  
  results <- left_join(results, test2, by = "Variable")
  return(results)
}

# Combine results across cohorts
df.or <- bind_rows(
  get_or(charls.cross) |> mutate(cohort = "CHARLS"),
  get_or(elsa.cross) |> mutate(cohort = "ELSA"),
  get_or(hrs.cross) |> mutate(cohort = "HRS"),
  get_or(mhas.cross) |> mutate(cohort = "MHAS"),
  get_or(share.cross) |> mutate(cohort = "SHARE")
) |>
  mutate(
    Variable = case_when(
      Variable == "No disease" ~ "No disorders",
      Variable == "Chronic lung disease" ~ "Lung disease",
      Variable == "Mental disease only" ~ "Psychological disorder",
      Variable == "Cognitive impairment only" ~ "Cognitive disorder",
      Variable == "Physical-Mental comorbidity" ~ "Physical–psychological multimorbidity",
      Variable == "Physical-Cognitive comorbidity" ~ "Physical–cognitive multimorbidity",
      Variable == "Mental-Cognitive comorbidity" ~ "Psychological–cognitive multimorbidity",
      Variable == "Physical-Mental-Cognitive comorbidity" ~ "Physical–psychological–cognitive multimorbidity",
      TRUE ~ Variable
    ),
    Variable = factor(
      Variable,
      levels = c(
        "No disorders", "Hypertension", "Diabetes", "Cancer", "Lung disease", "Heart disease",
        "Stroke", "Arthritis", "Psychological disorder", "Cognitive disorder",
        "Physical–psychological multimorbidity", "Physical–cognitive multimorbidity",
        "Psychological–cognitive multimorbidity", "Physical–psychological–cognitive multimorbidity"
      )
    )
  )

# Meta-analysis: part 1
df.or1 <- df.or |>
  dplyr::filter(Variable %in% c(
    "No disorders", "Hypertension", "Diabetes", "Cancer", 
    "Heart disease", "Psychological disorder", "Arthritis"
  ))

m.bin <- metabin(
  event.e = e.n, n.e = n1, event.c = e.n2, n.c = n2,
  studlab = cohort, data = df.or1, sm = "OR", method = "MH",
  MH.exact = TRUE, fixed = FALSE, random = TRUE,
  method.tau = "PM", hakn = TRUE, overall = FALSE,
  subgroup = Variable, overall.hetstat = FALSE, test.subgroup = FALSE,
  outclab = "Disability", label.e = "Chronic diseases",
  label.c = "No disease", title = ""
)

study_colors <- c(
  "HRS" = "#4E79A7", "CHARLS" = "#F28E2B", "ELSA" = "#E15759",
  "MHAS" = "#76B7B2", "SHARE" = "#59A14F"
)
color_vector <- study_colors[m.bin$studlab]

pdf("./results/20230321/or_forest1.pdf", width = 10, height = 13)
forest(
  m.bin, print.subgroup.name = FALSE, subgroup = TRUE,
  leftlabs = c("Study", "Effect Size", "Std. Error"),
  rightcols = c("effect", "ci", "w.random"), xlab = "Effect Size",
  weight.study = "random", comb.random = TRUE,
  col.square = color_vector, col.square.lines = color_vector,
  col.diamond.random = "black", col.study = color_vector,
  arrow.type = "closed", arrow.length = 0.1, cex = 1
)
dev.off()

# Meta-analysis: part 2 (remaining 7 conditions)
df.or2 <- df.or |>
  dplyr::filter(Variable %in% c(
    "No disorders", "Lung disease", "Stroke", "Cognitive disorder",
    "Physical–psychological multimorbidity", "Physical–cognitive multimorbidity",
    "Psychological–cognitive multimorbidity", "Physical–psychological–cognitive multimorbidity"
  ))

m.bin <- metabin(
  event.e = e.n, n.e = n1, event.c = e.n2, n.c = n2,
  studlab = cohort, data = df.or2, sm = "OR", method = "MH",
  MH.exact = TRUE, fixed = FALSE, random = TRUE,
  method.tau = "PM", hakn = TRUE, overall = FALSE,
  subgroup = Variable, overall.hetstat = FALSE, test.subgroup = FALSE,
  outclab = "Disability", label.e = "Chronic diseases",
  label.c = "No disease", title = ""
)

color_vector <- study_colors[m.bin$studlab]

pdf("./results/20230321/or_forest2.pdf", width = 10, height = 15)
forest(
  m.bin, print.subgroup.name = FALSE, subgroup = TRUE,
  leftlabs = c("Study", "Effect Size", "Std. Error"),
  rightcols = c("effect", "ci", "w.random"), xlab = "Effect Size",
  weight.study = "random", comb.random = TRUE,
  col.square = color_vector, col.square.lines = color_vector,
  col.diamond.random = "black", col.study = color_vector,
  arrow.type = "closed", arrow.length = 0.1, cex = 1
)
dev.off()

