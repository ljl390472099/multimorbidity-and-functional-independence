

extract_OR <- function(data,formula,extract_n){

    fit.log <- glm(formula,
                   family = "binomial",
                   data = data)
    

 
  summary_fit <- summary(fit.log)
  
  last_13_coef <- tail(summary_fit$coefficients, extract_n)
  

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
  )
  return(results)
}
