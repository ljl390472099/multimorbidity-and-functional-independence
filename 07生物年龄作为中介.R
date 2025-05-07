# charls 生物年龄中介分析

library(nhanesR)
library(BioAge)
library(tidyverse)
library(mediation)

#-------加入生物年龄作为中介---------
library(charlsR)

bio_charls <- dex_biologicalAge(waves = 2011)


# 整理结果

tidy_med.charls <- function(data, med, binary = FALSE) {
  
  # 定义变量
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
      Left_Join(bio_charls,by= c("ID"="id")) |> 
      mutate(aa = BiologicalAge - age,
             aa.c = as.integer(aa > 0),
             aa.c2 = as.integer(aa > -5),
             aa.c3 = as.integer(aa > 5),
             aa.c4 = as.integer(aa > 10)) |> 
      dplyr::select(ID,age.c,sex,marr,edu,exercise,smoke,drink,income.c,
                    bmi.c,BiologicalAge,detailed_category,aa,aa.c,adl) |> 
      filter(detailed_category == "No disease" | detailed_category == dis_name[i]) |> 
      mutate(detailed_category = ifelse(detailed_category == "No disease",0,1)) |> 
      mutate(detailed_category = factor(detailed_category)) |> 
      na.omit()
    
    
    # 使用 tryCatch 包裹每次迭代的逻辑
    tryCatch({
      
      # 拟合模型
      if (binary) {
        model.m <- eval(bquote(glm(.(f_m), family = binomial, data = med_nhanes))) 
      } else { 
        model.m <- eval(bquote(lm(.(f_m), data = med_nhanes)))
      }
      
      model.y <- eval(bquote(glm(.(f_y), family = binomial, data = med_nhanes)))
      
      # 中介分析
      result <- mediate(model.m, model.y, 
                        treat = "detailed_category", 
                        mediator = mediator, 
                        sims = 1000, boot = TRUE)
      
      # 提取结果
      summary_result <- summary(result)
      
      # 将需要的值添加到结果数据框
      results_summary <- rbind(
        results_summary,
        data.frame(
          Treatment = dis_name[i],  # 使用对应的疾病名称
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
      # 如果出现错误，打印错误信息并继续下一个循环
      message("Error in iteration ", i, ": ", e$message)
    })
  }
  
  return(results_summary)
}

# 添加样本量的函数
tidy_med.charls <- function(data, med, binary = FALSE) {
  
  # 定义变量
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
    SampleSize = integer(), # 添加样本量列
    stringsAsFactors = FALSE
  )
  
  set.seed(123)
  for (i in seq_along(disease)) {
    
    mediator <- med
    
    med_nhanes <- data |> 
      Left_Join(bio_charls, by = c("ID" = "id")) |> 
      mutate(aa = BiologicalAge - age,
             aa.c = as.integer(aa > 0),
             aa.c2 = as.integer(aa > -5),
             aa.c3 = as.integer(aa > 5),
             aa.c4 = as.integer(aa > 10)) |> 
      dplyr::select(ID, age.c, sex, marr, edu, exercise, smoke, drink, income.c,
                    bmi.c, BiologicalAge, detailed_category, aa, aa.c, adl) |> 
      filter(detailed_category == "No disease" | detailed_category == dis_name[i]) |> 
      mutate(detailed_category = ifelse(detailed_category == "No disease", 0, 1)) |> 
      mutate(detailed_category = factor(detailed_category)) |> 
      na.omit()
    
    sample_size <- nrow(med_nhanes) # 计算样本量
    
    # 使用 tryCatch 包裹每次迭代的逻辑
    tryCatch({
      
      # 拟合模型
      if (binary) {
        model.m <- eval(bquote(glm(.(f_m), family = binomial, data = med_nhanes))) 
      } else { 
        model.m <- eval(bquote(lm(.(f_m), data = med_nhanes)))
      }
      
      model.y <- eval(bquote(glm(.(f_y), family = binomial, data = med_nhanes)))
      
      # 中介分析
      result <- mediate(model.m, model.y, 
                        treat = "detailed_category", 
                        mediator = mediator, 
                        sims = 1000, boot = TRUE)
      
      # 提取结果
      summary_result <- summary(result)
      
      # 将需要的值添加到结果数据框
      results_summary <- rbind(
        results_summary,
        data.frame(
          Treatment = dis_name[i],  # 使用对应的疾病名称
          Mediator = mediator,
          EffectType = "ACME (Average)",
          Estimate = summary_result$d.avg,
          LowerCI = summary_result$d.avg.ci[1],
          UpperCI = summary_result$d.avg.ci[2],
          PValue = summary_result$d.avg.p,
          SampleSize = sample_size # 添加样本量
        ),
        data.frame(
          Treatment = dis_name[i],
          Mediator = mediator,
          EffectType = "ADE (Average)",
          Estimate = summary_result$z.avg,
          LowerCI = summary_result$z.avg.ci[1],
          UpperCI = summary_result$z.avg.ci[2],
          PValue = summary_result$z.avg.p,
          SampleSize = sample_size # 添加样本量
        ),
        data.frame(
          Treatment = dis_name[i],
          Mediator = mediator,
          EffectType = "Total Effect",
          Estimate = summary_result$tau.coef,
          LowerCI = summary_result$tau.ci[1],
          UpperCI = summary_result$tau.ci[2],
          PValue = summary_result$tau.p,
          SampleSize = sample_size # 添加样本量
        ),
        data.frame(
          Treatment = dis_name[i],
          Mediator = mediator,
          EffectType = "Prop. Mediated (Average)",
          Estimate = summary_result$n.avg,
          LowerCI = summary_result$n.avg.ci[1],
          UpperCI = summary_result$n.avg.ci[2],
          PValue = summary_result$n.avg.p,
          SampleSize = sample_size # 添加样本量
        )
      )
      
    }, error = function(e) {
      # 如果出现错误，打印错误信息并继续下一个循环
      message("Error in iteration ", i, ": ", e$message)
    })
  }
  
  return(results_summary)
}

tidy_med.charls <- function(data, med, binary = FALSE) {
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
  
  results_summary <- data.frame(
    Treatment = character(),
    Mediator = character(),
    EffectType = character(),
    Estimate = numeric(),
    LowerCI = numeric(),
    UpperCI = numeric(),
    PValue = numeric(),
    SampleSize = integer(),
    stringsAsFactors = FALSE
  )
  
  set.seed(123)
  for (i in seq_along(disease)) {
    tryCatch({
      med_nhanes <- data |> 
        filter(detailed_category == "No disease" | detailed_category == dis_name[i]) |> 
        mutate(detailed_category = ifelse(detailed_category == "No disease", 0, 1)) |> 
        mutate(detailed_category = factor(detailed_category)) |> 
        na.omit()
      
      sample_size <- nrow(med_nhanes)  # 计算样本量
      
      if (sample_size < 10) {
        message("Skipping disease: ", dis_name[i], " due to small sample size")
        next
      }
      
      f_m <- as.formula(paste0(med, "~ sex + marr + edu + 
                               exercise + smoke + drink + income.c + bmi.c + detailed_category"))
      f_y <- as.formula(paste0("adl ~", med, "+ sex + marr + edu + 
                               exercise + smoke + drink + income.c + bmi.c + detailed_category"))
      
      # 模型拟合
      if (binary) {
        model.m <- eval(bquote(glm(.(f_m), family = binomial, data = med_nhanes))) 
      } else { 
        model.m <- eval(bquote(lm(.(f_m), data = med_nhanes)))
      }
      model.y <- eval(bquote(glm(.(f_y), family = binomial, data = med_nhanes)))
      
      # 中介分析
      result <- mediate(model.m, model.y, 
                        treat = "detailed_category", 
                        mediator = med, 
                        sims = 1000, boot = TRUE)
      
      summary_result <- summary(result)
      
      # 保存结果
      results_summary <- rbind(
        results_summary,
        data.frame(
          Treatment = dis_name[i],
          Mediator = med,
          EffectType = "ACME (Average)",
          Estimate = summary_result$d.avg,
          LowerCI = summary_result$d.avg.ci[1],
          UpperCI = summary_result$d.avg.ci[2],
          PValue = summary_result$d.avg.p,
          SampleSize = sample_size
        ),
        data.frame(
          Treatment = dis_name[i],
          Mediator = med,
          EffectType = "ADE (Average)",
          Estimate = summary_result$z.avg,
          LowerCI = summary_result$z.avg.ci[1],
          UpperCI = summary_result$z.avg.ci[2],
          PValue = summary_result$z.avg.p,
          SampleSize = sample_size
        ),
        data.frame(
          Treatment = dis_name[i],
          Mediator = med,
          EffectType = "Total Effect",
          Estimate = summary_result$tau.coef,
          LowerCI = summary_result$tau.ci[1],
          UpperCI = summary_result$tau.ci[2],
          PValue = summary_result$tau.p,
          SampleSize = sample_size
        ),
        data.frame(
          Treatment = dis_name[i],
          Mediator = med,
          EffectType = "Prop. Mediated (Average)",
          Estimate = summary_result$n.avg,
          LowerCI = summary_result$n.avg.ci[1],
          UpperCI = summary_result$n.avg.ci[2],
          PValue = summary_result$n.avg.p,
          SampleSize = sample_size
        )
      )
      
    }, error = function(e) {
      message("Error in iteration ", i, ": ", e$message)
    })
  }
  
  return(results_summary)
}



med.bio.charls <- tidy_med.charls(charls.long,med = "aa.c",binary = T) |> 
  rbind(tidy_med.charls(charls.long,med = "aa.c2",binary = T)) |>
  rbind(tidy_med.charls(charls.long,med = "aa.c3",binary = T)) |> 
  rbind(tidy_med.charls(charls.long,med = "aa.c4",binary = T)) |> 
  rbind(tidy_med.charls(charls.long,med = "aa")) |> 
  rbind(tidy_med.charls(charls.long,med = "BiologicalAge")) 

med.bio.charls <- tidy_med.charls(charls.long,med = "BiologicalAge")

write.csv(med.bio.charls,"./results/df.med.bio_charls.csv")

tidy_med.charls.sen <- function(data, med, binary = FALSE) {
  
  # 定义变量
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
  
  results_list <- list()
  
  set.seed(123)
  for (i in seq_along(disease)) {
    
    mediator <- med
    
    med_nhanes <- data |> 
      left_join(bio_charls, by = c("ID" = "id")) |> 
      mutate(aa = BiologicalAge - age,
             aa.c = as.integer(aa > 0),
             aa.c2 = as.integer(aa > -5),
             aa.c3 = as.integer(aa > 5),
             aa.c4 = as.integer(aa > 10)) |> 
      dplyr::select(ID, age.c, sex, marr, edu, exercise, smoke, drink, income.c,
                    bmi.c, BiologicalAge, detailed_category, aa, aa.c, adl) |> 
      filter(detailed_category == "No disease" | detailed_category == dis_name[i]) |> 
      mutate(detailed_category = ifelse(detailed_category == "No disease", 0, 1)) |> 
      mutate(detailed_category = factor(detailed_category)) |> 
      na.omit()
    
    tryCatch({
      # 拟合模型
      if (binary) {
        model.m <- eval(bquote(glm(.(f_m), family = binomial, data = med_nhanes))) 
      } else { 
        model.m <- eval(bquote(lm(.(f_m), data = med_nhanes)))
      }
      
      model.y <- eval(bquote(glm(.(f_y), family = binomial("probit"), data = med_nhanes)))
      
      # 中介分析
      result <- mediate(model.m, model.y, 
                        treat = "detailed_category", 
                        mediator = mediator, 
                        sims = 1000, boot = TRUE)
      
      # 保存结果到列表
      results_list[[dis_name[i]]] <- result
      
    }, error = function(e) {
      # 如果出现错误，打印错误信息并继续下一个循环
      message("Error in iteration ", i, ": ", e$message)
    })
  }
  
  return(results_list)
}


med.result.charls <- tidy_med.charls.sen(charls.long,med = "BiologicalAge")


set.seed(123)
sens.out <- medsens(med.result.charls[["Hypertension"]], 
                    rho.by = 0.1, 
                    effect.type = "indirect", 
                    sims = 100) 
pdf("./results/network/med_Hypertension.pdf", width = 10,height = 10)
par(mfrow = c(2,2))
sens.out |> 
  plot(sens.par = "rho", 
       main = "Hypertension", 
       ylim = c(-0.2, 0.2))

sens.out |> 
  plot(
    sens.par = "R2", 
    r.type = "total", 
    sign.prod = "positive")
dev.off()

generate_med_sensitivity_plot <- function(disease_name, med_result, output_dir = "./results/network/") {
  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 设置随机种子以保证结果可重复
  set.seed(123)
  
  # 生成敏感性分析结果
  sens.out <- medsens(med_result[[disease_name]], 
                      rho.by = 0.1, 
                      effect.type = "indirect", 
                      sims = 100)
  
  # 定义输出文件路径
  output_file <- file.path(output_dir, paste0("med_", disease_name, ".pdf"))
  
  # 保存结果为 PDF
  pdf(output_file, width = 10, height = 10)
  par(mfrow = c(2, 2))
  
  # 绘制第一个子图
  sens.out |> 
    plot(sens.par = "rho", 
         main = "", 
         ylim = c(-0.2, 0.2))
 # mtext("A", side = 3, line = 0.5, adj = -0.2, font = 2)
  
  # 绘制第二个子图
  sens.out |> 
    plot(sens.par = "R2", 
         r.type = "total", 
         sign.prod = "positive")
 # mtext("B", side = 3, line = 0.5, adj = -0.2, font = 2)
  
  dev.off()
  
  cat("Sensitivity plot saved to:", output_file, "\n")
}

dis.med.name <- names(med.result.charls)

for (i in dis.med.name){
  generate_med_sensitivity_plot(disease_name = i,
                                med_result = med.result.charls)
}


