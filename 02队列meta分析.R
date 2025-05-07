# RD meta分析

if(!require("pacman"))install.packages("pacman")
pacman::p_load(tidyverse, haven, psych,lme4,poLCA,forestplot,table1,
               gridExtra,ggpubr,autoReg,readxl,marginaleffects)

source("./function/extract_RD.R")
source("./function/tidy_RD.R")

load("./data/charls.RData")
load("./data/elsa.RData")
load("./data/hrs.RData")
load("./data/mhas.RData")
load("./data/share.RData")

#----整理RD-----

elsa.long <- elsa.long |> 
  filter(p.na <= 0.2) |> 
  filter(na.adl < 4) 

hrs.long <- hrs.long |> 
  filter(p.na <= 0.2) |> 
  filter(na.adl < 5) 

mhas.long <- mhas.long |> 
  filter(p.na <= 0.2) |> 
  filter(na.adl < 2)

share.long <- share.long |> 
  left_join(cov.share |> 
              dplyr::select(mergeid,p.na),by="mergeid") |> 
  filter(p.na <= 0.2) |> 
  filter(na.adl < 4)

f1 <- adl ~ age + sex + marr + edu + exercise + smoke + drink + income.c+ bmi.c+
  detailed_category

f2 <- adl ~  sex + marr + edu + exercise + smoke + drink + income.c+ bmi.c+
  detailed_category

f3 <- adl ~  age  + marr + edu + exercise + smoke + drink + income.c+ bmi.c+
  detailed_category

rd.charls.long <- tidy_RD(data = charls.long,f1,f2,f3,disease_n = "detailed_category")
rd.elsa.long <- tidy_RD(data = elsa.long,f1,f2,f3,disease_n = "detailed_category")
rd.hrs.long <- tidy_RD(data = hrs.long,f1,f2,f3,disease_n = "detailed_category")
rd.mhas.long <- tidy_RD(data = mhas.long,f1,f2,f3,disease_n = "detailed_category")
rd.share.long <- tidy_RD(data = share.long,f1,f2,f3,disease_n = "detailed_category")

df.rd <- rd.charls.long |> 
  mutate(cohort = "CHARLS") |> 
  rbind(rd.elsa.long |> 
          mutate(cohort = "ELSA")) |> 
  rbind(rd.hrs.long |> 
          mutate(cohort = "HRS")) |> 
  rbind(rd.mhas.long |> 
          mutate(cohort = "MHAS")) |> 
  rbind(rd.share.long |> 
          mutate(cohort = "SHARE")) |> 
  mutate(contrast = str_remove(contrast, " - No disease")) |> 
  mutate(contrast = case_when(contrast == "No disease" ~ "No disorders",
                              contrast == "Chronic lung disease" ~ "Lung disease",
                              contrast == "Mental disease only" ~ "Psychological disorder",
                              contrast == "Cognitive impairment only" ~ "Cognitive disorder",
                              contrast == "Physical-Mental comorbidity" ~ "Physical–psychological multimorbidity",
                              contrast == "Physical-Cognitive comorbidity" ~ "Physical–cognitive multimorbidity",
                              contrast == "Mental-Cognitive comorbidity" ~ "Psychological–cognitive multimorbidity",
                              contrast == "Physical-Mental-Cognitive comorbidity" ~ "Physical–psychological–cognitive multimorbidity",
                              TRUE ~ contrast),
         contrast = factor(contrast, 
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


#-----meta分析-------

# 加载必要的包
library(meta)
library(metafor)

df.meta <- df.rd |> 
  filter(group == "Total") |> 
  dplyr::select(contrast,RD,se,cohort) 

df.meta1 <- df.meta |> 
  dplyr::filter(contrast %in% c("No disorders","Hypertension",
                                "Heart disease",
                                "Cancer",
                                "Cognitive disorder",
                                "Physical–cognitive multimorbidity",
                                "Psychological–cognitive multimorbidity"))
# 使用 meta 包
meta_rd <- metagen(
  TE = RD,  # 效应值
  seTE = se,  # 标准误差
  studlab = cohort,   # 研究名称
  data=df.meta1,
  sm = "SMD",                # 效应大小类型 (Standardized Mean Difference)
  random = TRUE,        # 使用随机效应模型
  common = FALSE,         # 不使用固定效应模型
  overall = F,
  overall.hetstat =F,
  test.subgroup = F,
  outclab = "Disability",
  label.e = "Chronic diseases",
  label.c = "No disease",
  title = "",
  subgroup = contrast
)

# 显示结果
print(meta_rd)

# 绘制森林图
forest(meta_rd, 
       #layout = "BMJ",
       print.subgroup.name = F,
       subgroup = TRUE,
       leftlabs = c("Study", "Effect Size", "Std. Error"),
       rightcols = c("effect", "ci", "w.random"),
       xlab = "Effect Size",
       weight.study = "random",
       comb.random = TRUE,
       smlab = "Risk differences",
       col.diamond = "blue",                 # 汇总值（菱形）的颜色
       col.study = "black",                  # 单独研究点的颜色
       col.square = "black",                 # 非汇总点（方点）颜色
       col.square.lines = "black",           # 方点的线段颜色
       lty.summary = 5                # 设置 x=1 的线为虚线
       )


metafor::forest(meta_rd, 
       lty.random = 10, 
       print.subgroup.name = FALSE,                 # 是否打印子组名称
       subgroup = TRUE,                             # 启用子组分析
       leftlabs = c("Study", "Effect Size", "Std. Error"),  # 左侧标签
       rightcols = c("effect", "ci", "w.random"),  # 右侧列
       xlab = "Effect Size",                       # x轴标签
       weight.study = "random",                    # 权重类型
       comb.random = TRUE,                         # 使用随机效应模型
       col.diamond = "black",                       # 汇总值（菱形）的颜色
       col.square = rainbow(length(meta_rd$TE)),   # 每个研究一个不同颜色
       col.square.lines = rainbow(length(meta_rd$TE)), # 方点线段颜色
       cex = 1,                                    # 字体大小调整
       col.study = rainbow(length(meta_rd$studlab)), # 每个研究颜色
       pch = 15 ,                                   # 固定为方点（或其他形状）
       psize = 1,                                   # 固定方点大小
       lwd.square = 0.2
)


study_colors <- c(
  "HRS" = "#4E79A7",
  "CHARLS" = "#F28E2B",
  "ELSA" = "#E15759",
  "MHAS" = "#76B7B2",
  "SHARE" = "#59A14F"
)
df.meta <- df.rd |> 
  filter(group == "Total") |> 
  dplyr::select(contrast,RD,se,cohort) 

df.meta1 <- df.meta |> 
  dplyr::filter(contrast %in% c("No disorders","Hypertension",
                                "Heart disease",
                                "Cancer",
                                "Cognitive disorder",
                                "Physical–cognitive multimorbidity",
                                "Psychological–cognitive multimorbidity"))
# 为每个研究的颜色创建对应向量
color_vector <- study_colors[meta_rd$studlab]

pdf("./results/20230321/rd_forest1.pdf", width = 10,height = 13)
# 绘制森林图
forest(meta_rd, 
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
       smlab = "Risk differences",
       #xlim = c(-0.5,0.3),
       cex = 1)                                     # 字体大小调整

dev.off()


#-----第二段森林图绘制-----

df.meta1 <- df.meta |> 
  dplyr::filter(contrast %in% c("No disorders",
                                "Diabetes",
                                "Lung disease",
                                "Stroke",
                                "Arthritis",
                                "Psychological disorder",
                               
                                # Comorbidities
                                "Physical–psychological multimorbidity",
                                
                                "Physical–psychological–cognitive multimorbidity"))
# 使用 meta 包
meta_rd <- metagen(
  TE = RD,  # 效应值
  seTE = se,  # 标准误差
  studlab = cohort,   # 研究名称
  data=df.meta1,
  sm = "SMD",                # 效应大小类型 (Standardized Mean Difference)
  random = TRUE,        # 使用随机效应模型
  common = FALSE,         # 不使用固定效应模型
  overall = F,
  overall.hetstat =F,
  test.subgroup = F,
  outclab = "Disability",
  label.e = "Chronic diseases",
  label.c = "No disease",
  title = "",
  subgroup = contrast
)

# 为每个研究的颜色创建对应向量
color_vector <- study_colors[meta_rd$studlab]

pdf("./results/20230321/rd_forest2.pdf", width = 10,height = 15)
# 绘制森林图
forest(meta_rd, 
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
       smlab = "Risk differences",
       #xlim = c(-0.5,0.3),
       cex = 1)                                     # 字体大小调整

dev.off()
