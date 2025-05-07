if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, haven, psych, lme4, poLCA, forestplot, table1,
               gridExtra, ggpubr, autoReg, readxl, marginaleffects,
               meta, metafor)

source("./function/extract_RD.R")
source("./function/tidy_RD.R")

# Load datasets
datasets <- list(
  CHARLS = "./data/charls.RData",
  ELSA   = "./data/elsa.RData",
  HRS    = "./data/hrs.RData",
  MHAS   = "./data/mhas.RData",
  SHARE  = "./data/share.RData"
)
lapply(datasets, load)

# Filter data
filter_data <- function(data, max_p_na, max_na_adl, join_cov = NULL) {
  if (!is.null(join_cov)) {
    data <- data %>%
      left_join(join_cov %>% select(mergeid, p.na), by = "mergeid")
  }
  data %>% filter(p.na <= max_p_na, na.adl < max_na_adl)
}

elsa.long  <- filter_data(elsa.long, 0.2, 4)
hrs.long   <- filter_data(hrs.long, 0.2, 5)
mhas.long  <- filter_data(mhas.long, 0.2, 2)
share.long <- filter_data(share.long, 0.2, 4, cov.share)

# Model formulas
f1 <- adl ~ age + sex + marr + edu + exercise + smoke + drink + income.c + bmi.c + detailed_category
f2 <- adl ~ sex + marr + edu + exercise + smoke + drink + income.c + bmi.c + detailed_category
f3 <- adl ~ age + marr + edu + exercise + smoke + drink + income.c + bmi.c + detailed_category

# Run tidy_RD
rd_list <- list(
  CHARLS = tidy_RD(charls.long, f1, f2, f3, "detailed_category"),
  ELSA   = tidy_RD(elsa.long, f1, f2, f3, "detailed_category"),
  HRS    = tidy_RD(hrs.long, f1, f2, f3, "detailed_category"),
  MHAS   = tidy_RD(mhas.long, f1, f2, f3, "detailed_category"),
  SHARE  = tidy_RD(share.long, f1, f2, f3, "detailed_category")
)

# Combine all cohorts
df.rd <- bind_rows(lapply(names(rd_list), function(name) {
  rd_list[[name]] %>% mutate(cohort = name)
})) %>%
  mutate(
    contrast = str_remove(contrast, " - No disease"),
    contrast = case_when(
      contrast == "No disease" ~ "No disorders",
      contrast == "Chronic lung disease" ~ "Lung disease",
      contrast == "Mental disease only" ~ "Psychological disorder",
      contrast == "Cognitive impairment only" ~ "Cognitive disorder",
      contrast == "Physical-Mental comorbidity" ~ "Physical–psychological multimorbidity",
      contrast == "Physical-Cognitive comorbidity" ~ "Physical–cognitive multimorbidity",
      contrast == "Mental-Cognitive comorbidity" ~ "Psychological–cognitive multimorbidity",
      contrast == "Physical-Mental-Cognitive comorbidity" ~ "Physical–psychological–cognitive multimorbidity",
      TRUE ~ contrast
    ),
    contrast = factor(contrast, levels = c(
      "No disorders", "Hypertension", "Diabetes", "Cancer", "Lung disease",
      "Heart disease", "Stroke", "Arthritis", "Psychological disorder",
      "Cognitive disorder", "Physical–psychological multimorbidity",
      "Physical–cognitive multimorbidity", "Psychological–cognitive multimorbidity",
      "Physical–psychological–cognitive multimorbidity"
    ))
  )

# Select data for meta-analysis
selected_contrasts <- c(
  "No disorders", "Hypertension", "Heart disease", "Cancer",
  "Cognitive disorder", "Physical–cognitive multimorbidity",
  "Psychological–cognitive multimorbidity"
)

df.meta <- df.rd %>%
  filter(group == "Total", contrast %in% selected_contrasts) %>%
  select(contrast, RD, se, cohort)

# Meta-analysis using 'meta' package
meta_rd <- metagen(
  TE = RD, seTE = se, studlab = cohort, data = df.meta,
  sm = "SMD", random = TRUE, common = FALSE,
  overall = FALSE, overall.hetstat = FALSE,
  test.subgroup = FALSE, outclab = "Disability",
  label.e = "Chronic diseases", label.c = "No disease",
  subgroup = contrast
)

print(meta_rd)

# Color palette
study_colors <- c(
  "HRS" = "#4E79A7", "CHARLS" = "#F28E2B", "ELSA" = "#E15759",
  "MHAS" = "#76B7B2", "SHARE" = "#59A14F"
)
color_vector <- study_colors[meta_rd$studlab]

# Save forest plot
pdf("./results/20230321/rd_forest1.pdf", width = 10, height = 13)
forest(meta_rd,
       print.subgroup.name = FALSE,
       subgroup = TRUE,
       leftlabs = c("Study", "Effect Size", "Std. Error"),
       rightcols = c("effect", "ci", "w.random"),
       xlab = "Effect Size",
       weight.study = "random",
       comb.random = TRUE,
       col.square = color_vector,
       col.square.lines = color_vector,
       col.diamond.random = "black",
       col.study = color_vector,
       arrow.type = "closed",
       arrow.length = 0.1,
       smlab = "Risk differences",
       cex = 1)
dev.off()

