
# meta 元回归分析

library(tidyverse)
# 加载必要的包
library(meta)
library(metafor)
library(readxl)

df.meta2 <-  df.meta |> 
  merge(read_xlsx("./data/meta_data.xlsx"),by="cohort")


meta_rd2 <- metagen(
  TE = RD,  # 效应值
  seTE = se,  # 标准误差
  studlab = cohort,   # 研究名称
  data=df.meta2,
  sm = "SMD",                # 效应大小类型 (Standardized Mean Difference)
  random = TRUE,        # 使用随机效应模型
  common = FALSE,         # 不使用固定效应模型
  overall = F,
  subgroup = contrast
)


meta.dis <- df.meta2 |> dplyr::select(contrast) |> pull() |> unique()

meta.cov <- names(df.meta2)[7:16]

m.qual.rep <- rma(yi = RD, 
                  sei = se, 
                  data = df.meta2[df.meta2$contrast == meta.dis[1],], 
                  method = "ML", 
                  mods = ~ meta.cov[1], 
                  test = "knha")

# 创建一个空的数据框用于存储结果
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

# 循环遍历每个对比和每个协变量
for (dis in meta.dis) {
  for (cov in meta.cov) {
    # 进行meta回归
    model <- tryCatch({
      rma(yi = RD,
          sei = se,
          data = df.meta2[df.meta2$contrast == dis, ],
          method = "ML",
          mods = as.formula(paste0("~", cov)),
          test = "knha")
    }, error = function(e) NULL)  # 捕获错误并返回NULL
    
    # 如果模型成功运行，提取结果
    if (!is.null(model)) {
      results <- rbind(
        results,
        data.frame(
          contrast = dis,
          covariate = cov,
          estimate = coef(model)[2], # 提取协变量的估计值
          ci.lb = model[["ci.lb"]][2], # 提取置信区间下界
          ci.ub = model[["ci.ub"]][2], # 提取置信区间上界
          p.value = coef(summary(model))[2, "pval"], # 提取P值
          tau = sqrt(model$tau2), # 提取 tau（tau 的平方根）
          I2 = model$I2, # 提取 I^2
          H2 = model$H2, # 提取 H^2
          R2 = model$R2, # 提取 R^2
          stringsAsFactors = FALSE
        )
      )
    }
  }
}
# 查看汇总结果
print(results)

write.csv(results,file = "./results/meta_regression_RD.csv")


#-----OR的meta元回归---------

df.meta2 <-  df.or |> 
  merge(read_xlsx("./data/meta_data.xlsx"),by="cohort") |> 
  mutate(se = (CI_Upper - CI_Lower) / (2 * qnorm(0.975)))

# 创建一个空的数据框用于存储结果
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

# 循环遍历每个对比和每个协变量
for (dis in meta.dis) {
  for (cov in meta.cov) {
    # 进行meta回归
    model <- tryCatch({
      rma(yi = OR,
          sei = se,
          data = df.meta2[df.meta2$Variable == dis, ],
          method = "ML",
          mods = as.formula(paste0("~", cov)),
          test = "knha")
    }, error = function(e) NULL)  # 捕获错误并返回NULL
    
    # 如果模型成功运行，提取结果
    if (!is.null(model)) {
      results <- rbind(
        results,
        data.frame(
          contrast = dis,
          covariate = cov,
          estimate = coef(model)[2], # 提取协变量的估计值
          ci.lb = model[["ci.lb"]][2], # 提取置信区间下界
          ci.ub = model[["ci.ub"]][2], # 提取置信区间上界
          p.value = coef(summary(model))[2, "pval"], # 提取P值
          tau = sqrt(model$tau2), # 提取 tau（tau 的平方根）
          I2 = model$I2, # 提取 I^2
          H2 = model$H2, # 提取 H^2
          R2 = model$R2, # 提取 R^2
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

write.csv(results,file = "./results/meta_regression_OR.csv")


