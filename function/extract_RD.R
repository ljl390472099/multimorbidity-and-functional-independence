library(marginaleffects)

# 自定义函数提取边际效应
extract_RD <- function(data,formula, start_var) {
  
  model <- glm(formula,
                 family = "binomial",
                 data = data)
  
  
  # 获取模型的所有变量名
  var_names <- attr(model$terms, "term.labels")
  
  # 找到 start_var 的起始位置
  start_idx <- which(var_names == start_var)
  if (length(start_idx) == 0) {
     stop("start_var not found in the model.")
   }
  
  # 筛选从 start_var 开始的变量
  selected_vars <- var_names[start_idx:length(var_names)]
  
  # 计算边际效应
  effects <- avg_comparisons(model)
  
  # 筛选出感兴趣的变量
  filtered_effects <- effects |> 
    dplyr::filter(term %in% selected_vars) |> 
    dplyr::mutate(p.value = 2 * pnorm(-abs(estimate / std.error))) |> # 根据 z 值计算 p 值
    dplyr::select(term, contrast,estimate, std.error, conf.low, conf.high, p.value)
  
  names(filtered_effects) <- c("disease","contrast","RD","se","ci_low","ci_high","P")
  return(filtered_effects)
}


