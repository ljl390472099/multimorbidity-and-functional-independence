# charls Biological Age Mediation Analysis

library(nhanesR)
library(BioAge)
library(tidyverse)
library(mediation)

#-------Adding biological age as mediator---------
library(charlsR)

bio_charls <- dex_biologicalAge(waves = 2011)

# Organizing results

tidy_med.charls <- function(data, med, binary = FALSE) {
  
  # Define variables
  disease <- c(
    "mental2", "Diabetes", "physical_mental", "mental_cog", "Hypertension", 
    "Arthritis", "Stroke", "lung2", "heart2", "Cancer", 
    "phy_cog", "cognition", "phy_cog_mental"
  )
  dis_name <- c(
    "Mental disease only", "Diabetes", "Physical-Mental comorbidity", 
    "Mental-Cognitive comorbidity", "Hypertension", "Arthritis", 
    "Stroke", "Chronic lung disease", "Heart disease", "Cancer", 
    "Physical-Cognitive comorbidity", "Cognitive impairment only", 
    "Physical-Mental-Cognitive comorbidity"
  )
  
  f_m <- as.formula(paste0(med, "~ sex + marr + edu + 
                 exercise + smoke + drink + income.c + bmi.c + detailed_category"))
  f_y <- as.formula(paste0("adl ~", med, "+ sex + marr + edu + 
                 exercise + smoke + drink + income.c + bmi.c + detailed_category"))
  
  results_summary <- data.frame(
    Treatment = character(),
    Mediator = character(),
    EffectType = character(),
    Estimate = numeric(),
    LowerCI = numeric(),
    UpperCI = numeric(),
    PValue = numeric(),
    stringsAsFactors = FALSE
  )
  
  set.seed(123)
  for (i in seq_along(disease)) {
    
    mediator <- med
    
    med_nhanes <- data |> 
      Left_Join(bio_charls, by = c("ID" = "id")) |> 
      mutate(ba = BiologicalAge - age) |>  # Focus on Biological Age
      dplyr::select(ID, age.c, sex, marr, edu, exercise, smoke, drink, income.c,
                    bmi.c, BiologicalAge, detailed_category, ba, adl) |> 
      filter(detailed_category == "No disease" | detailed_category == dis_name[i]) |> 
      mutate(detailed_category = ifelse(detailed_category == "No disease", 0, 1)) |> 
      mutate(detailed_category = factor(detailed_category)) |> 
      na.omit()
    
    # Using tryCatch to wrap each iteration's logic
    tryCatch({
      
      # Fit models
      if (binary) {
        model.m <- eval(bquote(glm(.(f_m), family = binomial, data = med_nhanes))) 
      } else { 
        model.m <- eval(bquote(lm(.(f_m), data = med_nhanes)))
      }
      
      model.y <- eval(bquote(glm(.(f_y), family = binomial, data = med_nhanes)))
      
      # Mediation analysis
      result <- mediate(model.m, model.y, 
                        treat = "detailed_category", 
                        mediator = mediator, 
                        sims = 1000, boot = TRUE)
      
      # Extract results
      summary_result <- summary(result)
      
      # Add required values to the results data frame
      results_summary <- rbind(
        results_summary,
        data.frame(
          Treatment = dis_name[i],  # Use corresponding disease name
          Mediator = mediator,
          EffectType = "ACME (Average)",
          Estimate = summary_result$d.avg,
          LowerCI = summary_result$d.avg.ci[1],
          UpperCI = summary_result$d.avg.ci[2],
          PValue = summary_result$d.avg.p
        ),
        data.frame(
          Treatment = dis_name[i],
          Mediator = mediator,
          EffectType = "ADE (Average)",
          Estimate = summary_result$z.avg,
          LowerCI = summary_result$z.avg.ci[1],
          UpperCI = summary_result$z.avg.ci[2],
          PValue = summary_result$z.avg.p
        ),
        data.frame(
          Treatment = dis_name[i],
          Mediator = mediator,
          EffectType = "Total Effect",
          Estimate = summary_result$tau.coef,
          LowerCI = summary_result$tau.ci[1],
          UpperCI = summary_result$tau.ci[2],
          PValue = summary_result$tau.p
        ),
        data.frame(
          Treatment = dis_name[i],
          Mediator = mediator,
          EffectType = "Prop. Mediated (Average)",
          Estimate = summary_result$n.avg,
          LowerCI = summary_result$n.avg.ci[1],
          UpperCI = summary_result$n.avg.ci[2],
          PValue = summary_result$n.avg.p
        )
      )
      
    }, error = function(e) {
      # If error occurs, print the error message and proceed to the next iteration
      message("Error in iteration ", i, ": ", e$message)
    })
  }
  
  return(results_summary)
}

