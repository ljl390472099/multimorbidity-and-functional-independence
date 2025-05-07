library(marginaleffects)

extract_RD <- function(data,formula, start_var) {
  
  model <- glm(formula,
                 family = "binomial",
                 data = data)
  
  

  var_names <- attr(model$terms, "term.labels")
  

  start_idx <- which(var_names == start_var)
  if (length(start_idx) == 0) {
     stop("start_var not found in the model.")
   }
  

  selected_vars <- var_names[start_idx:length(var_names)]
  

  effects <- avg_comparisons(model)
  

  filtered_effects <- effects |> 
    dplyr::filter(term %in% selected_vars) |> 
    dplyr::mutate(p.value = 2 * pnorm(-abs(estimate / std.error))) |>
    dplyr::select(term, contrast,estimate, std.error, conf.low, conf.high, p.value)
  
  names(filtered_effects) <- c("disease","contrast","RD","se","ci_low","ci_high","P")
  return(filtered_effects)
}


