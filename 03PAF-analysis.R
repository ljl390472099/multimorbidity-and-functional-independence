
# Load required libraries
library(AF)
library(patchwork)
source("./function/extract_RD.R")
source("./function/tidy_RD.R")
source("./function/run_AF_analysis.R")

# Define function to compute Attributable Fraction (AF)
get_af <- function(data, f) {
  library(dplyr)
  library(tidyr)
  
  if (!"detailed_category" %in% colnames(data)) {
    stop("Missing 'detailed_category' column in input data")
  }
  
  data <- data %>%
    mutate(id = row_number()) %>%
    pivot_wider(
      names_from = detailed_category,
      values_from = detailed_category,
      values_fn = length,
      values_fill = 0
    ) %>%
    as.data.frame()
  
  # Rename variables (if present)
  rename_map <- c(
    "Chronic lung disease" = "lung2",
    "Heart disease" = "heart2",
    "Mental disease only" = "mental2",
    "Cognitive impairment only" = "cognition",
    "Physical-Mental comorbidity" = "physical_mental",
    "Physical-Cognitive comorbidity" = "phy_cog",
    "Mental-Cognitive comorbidity" = "mental_cog",
    "Physical-Mental-Cognitive comorbidity" = "phy_cog_mental"
  )
  
  set.seed(123)
  existing_renames <- rename_map[names(rename_map) %in% colnames(data)]
  data <- data %>% rename(!!!setNames(names(existing_renames), existing_renames))
  
  # Adjust formula to include only existing variables
  f_vars <- all.vars(f)
  existing_vars <- f_vars[f_vars %in% colnames(data)]
  f <- reformulate(existing_vars[-1], response = existing_vars[1])
  
  # Chronic disease list
  chronic_diseases <- intersect(c("mental2", "Diabetes", "physical_mental",
                                  "mental_cog", "Hypertension", "Arthritis",      
                                  "Stroke", "lung2", "heart2",          
                                  "Cancer", "phy_cog", "cognition",      
                                  "phy_cog_mental"), colnames(data))
  
  fit.log <- glm(f, family = "binomial", data = data)
  
  results_df <- data.frame(Disease = character(),
                           AF = numeric(),
                           Std.Error = numeric(),
                           z_value = numeric(),
                           Pr_z = numeric(),
                           stringsAsFactors = FALSE)
  
  for (disease in chronic_diseases) {
    AF_result <- AFglm(object = fit.log, data = data, exposure = disease)
    AF_summary <- summary(AF_result)
    disease_result <- data.frame(Disease = disease, AF_summary$AF)
    results_df <- rbind(results_df, disease_result)
  }
  return(results_df)
}

# Define model formulas
f1 <- adl ~ age + sex + marr + edu + exercise + smoke + drink + income.c + bmi.c +
  Hypertension + Diabetes + Cancer + lung2 + heart2 + Stroke + Arthritis + 
  mental2 + cognition + physical_mental + phy_cog + mental_cog + phy_cog_mental

f2 <- adl ~ sex + marr + edu + exercise + smoke + drink + income.c + bmi.c +
  Hypertension + Diabetes + Cancer + lung2 + heart2 + Stroke + Arthritis + 
  mental2 + cognition + physical_mental + phy_cog + mental_cog + phy_cog_mental

f3 <- adl ~ age + marr + edu + exercise + smoke + drink + income.c + bmi.c +
  Hypertension + Diabetes + Cancer + lung2 + heart2 + Stroke + Arthritis + 
  mental2 + cognition + physical_mental + phy_cog + mental_cog + phy_cog_mental

# Example individual AF computation
af.result <- AFglm(object = fit.log, data = data, exposure = "Hypertension")
summary(af.result)

# Combine AF results across subgroups and cohorts
df.af <- get_af(charls.long, f = f1) |> 
  mutate(group = "Total") |> 
  rbind(get_af(charls.long[charls.long$age.c == "<60",], f = f2) |> mutate(group = "<60")) |> 
  rbind(get_af(charls.long[charls.long$age.c == "60-74",], f = f2) |> mutate(group = "60-74")) |> 
  rbind(get_af(charls.long[charls.long$age.c == "75+",], f = f2) |> mutate(group = "75+")) |>
  rbind(get_af(charls.long[charls.long$sex == "Male",], f = f3) |> mutate(group = "Male")) |>
  rbind(get_af(charls.long[charls.long$sex == "Female",], f = f3) |> mutate(group = "Female")) |> 
  mutate(cohort = "CHARLS (China)") |> 
  rbind(
    get_af(elsa.long, f = f1) |> mutate(group = "Total") |> 
    rbind(get_af(elsa.long[elsa.long$age.c == "<60",], f = f2) |> mutate(group = "<60")) |> 
    rbind(get_af(elsa.long[elsa.long$age.c == "60-74",], f = f2) |> mutate(group = "60-74")) |> 
    rbind(get_af(elsa.long[elsa.long$age.c == "75+",], f = f2) |> mutate(group = "75+")) |>
    rbind(get_af(elsa.long[elsa.long$sex == "Male",], f = f3) |> mutate(group = "Male")) |>
    rbind(get_af(elsa.long[elsa.long$sex == "Female",], f = f3) |> mutate(group = "Female")) |> 
    mutate(cohort = "ELSA (UK)")
  ) |> 
  rbind(
    get_af(hrs.long, f = f1) |> mutate(group = "Total") |> 
    rbind(get_af(hrs.long[hrs.long$age.c == "<60",], f = f2) |> mutate(group = "<60")) |> 
    rbind(get_af(hrs.long[hrs.long$age.c == "60-74",], f = f2) |> mutate(group = "60-74")) |> 
    rbind(get_af(hrs.long[hrs.long$age.c == "75+",], f = f2) |> mutate(group = "75+")) |>
    rbind(get_af(hrs.long[hrs.long$sex == "Male",], f = f3) |> mutate(group = "Male")) |>
    rbind(get_af(hrs.long[hrs.long$sex == "Female",], f = f3) |> mutate(group = "Female")) |> 
    mutate(cohort = "HRS (USA)")
  ) |> 
  rbind(
    get_af(mhas.long, f = f1) |> mutate(group = "Total") |> 
    rbind(get_af(mhas.long[mhas.long$age.c == "<60",], f = f2) |> mutate(group = "<60")) |> 
    rbind(get_af(mhas.long[mhas.long$age.c == "60-74",], f = f2) |> mutate(group = "60-74")) |> 
    rbind(get_af(mhas.long[mhas.long$age.c == "75+",], f = f2) |> mutate(group = "75+")) |>
    rbind(get_af(mhas.long[mhas.long$sex == "Male",], f = f3) |> mutate(group = "Male")) |>
    rbind(get_af(mhas.long[mhas.long$sex == "Female",], f = f3) |> mutate(group = "Female")) |> 
    mutate(cohort = "MHAS (Mexico)")
  ) |> 
  rbind(
    get_af(share.long, f = f1) |> mutate(group = "Total") |> 
    rbind(get_af(share.long[share.long$age.c == "<60",], f = f2) |> mutate(group = "<60")) |> 
    rbind(get_af(share.long[share.long$age.c == "60-74",], f = f2) |> mutate(group = "60-74")) |> 
    rbind(get_af(share.long[share.long$age.c == "75+",], f = f2) |> mutate(group = "75+")) |>
    rbind(get_af(share.long[share.long$sex == "Male",], f = f3) |> mutate(group = "Male")) |>
    rbind(get_af(share.long[share.long$sex == "Female",], f = f3) |> mutate(group = "Female")) |> 
    mutate(cohort = "SHARE (Europe)")
  )

# Optional: disease variable name mapping for visualization
disease <- c("mental2", "Diabetes", "physical_mental", "mental_cog", "Hypertension", 
             "Arthritis", "Stroke", "lung2", "heart2", "Cancer", 
             "phy_cog", "cognition", "phy_cog_mental")

dis_name <- c("Mental disease only", "Diabetes", "Physical-Mental comorbidity", 
              "Mental-Cognitive comorbidity", "Hypertension", "Arthritis", 
              "Stroke", "Chronic lung disease", "Heart disease", "Cancer", 
              "Physical-Cognitive comorbidity", "Cognitive impairment only", 
              "Physical-Mental-Cognitive comorbidity")
