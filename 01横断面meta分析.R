

if(!require("pacman"))install.packages("pacman")
pacman::p_load(tidyverse, haven, psych,lme4,poLCA,forestplot,table1,
               gridExtra,ggpubr,autoReg,readxl,marginaleffects)

library(rrtable)
load("./data/charls.RData")
load("./data/elsa.RData")
load("./data/hrs.RData")
load("./data/mhas.RData")
load("./data/share.RData")

# 加载必要的包
library(meta)
library(metafor)


#-----提取横断面OR----------

get_or <- function(data){
  fit.log <- glm(adl1 ~ age + sex + marr + edu + exercise + smoke + drink + income.c + bmi.c +
                   detailed_category,
                 family = "binomial",
                 data = data)
  # 提取summary对象
  summary_fit <- summary(fit.log)
  # 提取最后13个变量的系数
  last_13_coef <- tail(summary_fit$coefficients, 13)
  
  # 计算OR值和95%置信区间
  exp_coef <- exp(last_13_coef[, "Estimate"])
  exp_ci_lower <- exp(last_13_coef[, "Estimate"] - 1.96 * last_13_coef[, "Std. Error"])
  exp_ci_upper <- exp(last_13_coef[, "Estimate"] + 1.96 * last_13_coef[, "Std. Error"])
  p_values <- last_13_coef[, "Pr(>|z|)"]
  
  # 将结果合并成数据框
  results <- data.frame(
    Variable = rownames(last_13_coef),
    OR = exp_coef,
    CI_Lower = exp_ci_lower,
    CI_Upper = exp_ci_upper,
    P_Value = p_values
  ) |> 
    mutate(Variable = gsub("detailed_category","",Variable))
  
  test <- data |> 
    dplyr::select(detailed_category,adl1) |> 
    table() |> 
    as.data.frame() |> 
    pivot_wider(
      names_from = adl1,  # 列名来源
      values_from = Freq     # 数据来源
    ) |> 
    rename("Variable" = detailed_category,c.n = `0`,e.n = `1`) |> 
    mutate(n1 = rowSums(across(c(c.n, e.n))))
  
  test2 <- test[2:14,] |> 
    mutate(e.n2 = as.numeric(test[1,3]),
           n2 = as.numeric(test[1,4]))
  
  
  results <- left_join(results,test2,by=c("Variable"="Variable"))
  return(results)
}

df.or <- get_or(charls.cross) |> 
  mutate(cohort = "CHARLS") |> 
  rbind(get_or(elsa.cross) |> 
          mutate(cohort = "ELSA")) |> 
  rbind(get_or(hrs.cross) |> 
          mutate(cohort = "HRS")) |> 
  rbind(get_or(mhas.cross) |> 
          mutate(cohort = "MHAS")) |> 
  rbind(get_or(share.cross) |> 
          mutate(cohort = "SHARE")) |> 
  mutate(Variable = case_when(Variable == "No disease" ~ "No disorders",
                                       Variable == "Chronic lung disease" ~ "Lung disease",
                                       Variable == "Mental disease only" ~ "Psychological disorder",
                                       Variable == "Cognitive impairment only" ~ "Cognitive disorder",
                                       Variable == "Physical-Mental comorbidity" ~ "Physical–psychological multimorbidity",
                                       Variable == "Physical-Cognitive comorbidity" ~ "Physical–cognitive multimorbidity",
                                       Variable == "Mental-Cognitive comorbidity" ~ "Psychological–cognitive multimorbidity",
                                       Variable == "Physical-Mental-Cognitive comorbidity" ~ "Physical–psychological–cognitive multimorbidity",
                                       TRUE ~ Variable),
         Variable = factor(Variable, 
                                    levels = c(
                                      # Reference group
                                      "No disorders",
                                      # Single physical diseases
                                      "Hypertension",
                                      "Diabetes",
                                      "Cancer",
                                      "Lung disease",
                                      "Heart disease",
                                      "Stroke",
                                      "Arthritis",
                                      # Non-physical diseases
                                      "Psychological disorder",
                                      "Cognitive disorder",
                                      # Comorbidities
                                      "Physical–psychological multimorbidity",
                                      "Physical–cognitive multimorbidity",
                                      "Psychological–cognitive multimorbidity",
                                      "Physical–psychological–cognitive multimorbidity"
                                    )))
  

#-----横断面meta分析-------
df.or1 <- df.or |> 
  dplyr::filter(Variable %in% c("No disorders","Hypertension",
                                "Diabetes",
                                "Cancer",
                                "Heart disease",
                                "Psychological disorder",
                                "Arthritis"))

m.bin <- metabin(event.e = e.n, # 实验组
                 n.e = n1,
                 event.c = e.n2, # 对照组
                 n.c = n2,
                 studlab = cohort,
                 data = df.or1,
                 sm = "OR",
                 method = "MH",
                 MH.exact = TRUE,
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "PM",
                 hakn = TRUE,
                 overall = F,
                 subgroup = Variable,
                 overall.hetstat =F,
                 test.subgroup = F,
                 outclab = "Disability",
                 label.e = "Chronic diseases",
                 label.c = "No disease",
                 title = "")


forest(m.bin,
       #layout = "BMJ",
       subgroup = TRUE,
       print.subgroup.name = F,
       overall.hetstat = F,
       pch.random = 15,                          # 设置每个研究点的形状为方块（OR值）
       cex.random = 1,                         # 设置点的大小
       lwd.random = 2,                           # 设置线条宽度
       lty.random = 2,                    # 设置线条类型
       addpoly = T, 
       # smlab = "",
       # 显示汇总效应（菱形）
       weight.study = "random",
       type.subgroup.random = "square")


study_colors <- c(
  "HRS" = "#4E79A7",
  "CHARLS" = "#F28E2B",
  "ELSA" = "#E15759",
  "MHAS" = "#76B7B2",
  "SHARE" = "#59A14F"
)

# 为每个研究的颜色创建对应向量
color_vector <- study_colors[m.bin$studlab]

pdf("./results/20230321/or_forest1.pdf", width = 10,height = 13)
# 绘制森林图
forest(m.bin, 
       print.subgroup.name = FALSE,                  # 是否打印子组名称
       subgroup = TRUE,                              # 启用子组分析
       leftlabs = c("Study", "Effect Size", "Std. Error"),  # 左侧标签
       rightcols = c("effect", "ci", "w.random"),   # 右侧列
       xlab = "Effect Size",                        # x轴标签
       weight.study = "random",                     # 权重类型
       comb.random = TRUE,                          # 使用随机效应模型
       col.square = color_vector,                   # 方点颜色（研究对应颜色）
       col.square.lines = color_vector,             # 方点两端线段颜色
       col.diamond.random = "black",                 # 随机效应菱形颜色
       #lty.random = "dashed",                       # 随机效应汇总线虚线
       col.study = color_vector,
       arrow.type = "closed",
       arrow.length = 0.1,
       #xlim = c(-0.5,0.3),
       cex = 1)                                     # 字体大小调整

dev.off()


#-----后面6个分析-----


df.or1 <- df.or |> 
  dplyr::filter(Variable %in% c("No disorders", 
                                "Lung disease",
                                "Stroke",
                                "Cognitive disorder",
                                # Comorbidities
                                "Physical–psychological multimorbidity",
                                "Physical–cognitive multimorbidity",
                                "Psychological–cognitive multimorbidity",
                                "Physical–psychological–cognitive multimorbidity"))

m.bin <- metabin(event.e = e.n, # 实验组
                 n.e = n1,
                 event.c = e.n2, # 对照组
                 n.c = n2,
                 studlab = cohort,
                 data = df.or1,
                 sm = "OR",
                 method = "MH",
                 MH.exact = TRUE,
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "PM",
                 hakn = TRUE,
                 overall = F,
                 subgroup = Variable,
                 overall.hetstat =F,
                 test.subgroup = F,
                 outclab = "Disability",
                 label.e = "Chronic diseases",
                 label.c = "No disease",
                 title = "")


# 为每个研究的颜色创建对应向量
color_vector <- study_colors[m.bin$studlab]

pdf("./results/20230321/or_forest2.pdf", width = 10,height = 15)
# 绘制森林图
forest(m.bin, 
       print.subgroup.name = FALSE,                  # 是否打印子组名称
       subgroup = TRUE,                              # 启用子组分析
       leftlabs = c("Study", "Effect Size", "Std. Error"),  # 左侧标签
       rightcols = c("effect", "ci", "w.random"),   # 右侧列
       xlab = "Effect Size",                        # x轴标签
       weight.study = "random",                     # 权重类型
       comb.random = TRUE,                          # 使用随机效应模型
       col.square = color_vector,                   # 方点颜色（研究对应颜色）
       col.square.lines = color_vector,             # 方点两端线段颜色
       col.diamond.random = "black",                 # 随机效应菱形颜色
       #lty.random = "dashed",                       # 随机效应汇总线虚线
       col.study = color_vector,
       arrow.type = "closed",
       arrow.length = 0.1,
       #xlim = c(-0.5,0.3),
       cex = 1)                                     # 字体大小调整

dev.off()

