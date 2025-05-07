# 网络分析


library(networktools)
library(bootnet)
library(mgm)
library(dplyr)
library(magrittr)
library(NetworkComparisonTest)
library(psych)
library(qgraph)

library(IsingSampler)
library(IsingFit)


##回归掉协变量
library(tidyverse)

df_net <- elsa.long |>
  dplyr::select(hb , diabetes , cancer , lung , heart_new , stroke , psy.dep, arth,mmse.c,
                adl, age, sex,marr, edu, exercise,smoke,drink,income.c,bmi.c) |>
  mutate(across(everything(), as.numeric)) |> 
  na.omit()
df_net.charls <- charls.long |>
  dplyr::select(hb , diabetes , cancer , lung , heart , stroke , psy.dep, arth,memary.mmse,
                adl, age, sex,marr, edu, exercise,smoke,drink,income.c,bmi.c) |>
  mutate(across(everything(), as.numeric)) |> 
  na.omit()

df_net.hrs <- hrs.long |>
  dplyr::select(hp, diabetes , cancer , lung , heart_new , stroke , psy.dep, arth,memary.mmse,
                adl, age, sex,marr, edu, exercise,smoke,drink,income.c,bmi.c) |>
  mutate(across(everything(), as.numeric)) |> 
  na.omit()

df_net.mhas <- mhas.long |>
  dplyr::select(hb, diabetes , cancer , lung , hrtatt , stroke , dep, arth,mmse.c,
                adl, age, sex,marr, edu, exercise,smoke,drink,income.c,bmi.c) |>
  # mutate(across(everything(), as.numeric)) |> 
  mutate(across(c(hb, diabetes , cancer , lung , hrtatt , stroke , dep, arth,mmse.c,
                  adl), as.numeric)) |>
  na.omit()

df_net.share <- share.long |>
  dplyr::select(hb, diabetes , cancer , lung , heart , stroke , mental, arth,cognitive,
                adl, age, sex,marr, edu, exercise,smoke,drink,income.c,bmi.c) |>
  # mutate(across(everything(), as.numeric)) |> 
  mutate(across(c(hb, diabetes , cancer , lung , heart , stroke , mental, arth,cognitive,
                  adl), as.numeric)) |> 
  na.omit()

generate_network <- function(data,set.layout = F){
  set.seed(123)
  df_net <- data
  names(df_net)<-c("D1","D2","D3","D4","D5","D6","D7","D8","D9","ADL",
                   "age", "sex","marr", "edu", "exercise","smoke","drink","income.c","bmi.c")
  
  depression_mat <- as.data.frame(matrix(nrow=nrow(df_net),ncol=10))
  colnames(depression_mat)<-colnames(df_net)[1:10]#这个函数的意义是将每一个变量进行回归，然后将残差存储在depression_mat中
  for (i in 1:(ncol(depression_mat))){
    formula<-paste(colnames(df_net)[i],'~ age + sex + marr + edu + exercise + smoke + drink + income.c+ bmi.c')#这里的age和gender是协变量，可以根据实际情况进行修改
    lm1 <- glm(formula = formula , family = "binomial",data = df_net)#这里的family = "gaussian"是假定误差服从正态分布
    y_res<-lm1$residuals#这里的lm1$residuals是获取残差
    depression_mat[,i]<-y_res#将残差存储在depression_mat中
  }
  
  network <- as.matrix(depression_mat)
  p<-ncol(network)
  fit_obj <- mgm(data = network,
                 type = rep('g', p),
                 level = rep(1, p),
                 lambdaSel = 'CV',
                 ruleReg = 'OR',
                 pbar = TRUE)
  
  pred_obj <- predict(object = fit_obj,
                      data = network,
                      errorCon = 'R2')
  pred_obj$error#可预测性计算结果
  
  #估计网络
  Network <- estimateNetwork(data=network,
                             default = "EBICglasso",
                             tuning = 0.5,
                             corMethod="cor",#相关
                             corArgs=list(method="spearman",
                                          use="pairwise.complete.obs"))
  if (set.layout == F) {
    p.network <- plot(Network,
                      layout = "spring",
                      groups = groups,
                      label.cex = 1,#调整节点中标签的字体大小
                      label.color = 'black',
                      # 调整节点中标签的字体颜色
                      negDashed = T,
                      #负相关的边为虚线
                      legend=T,
                      nodeNames = items,
                      legend.cex = 0.5,#图例的字体大小
                      legend.mode = 'style1',#图例的风格，默认是'style1'
                      pie = pred_obj$error[,2])
  }else{
    p.network <- plot(Network,
                      layout = p.network.elsa[[2]]$layout,
                      groups = groups,
                      label.cex = 1,#调整节点中标签的字体大小
                      label.color = 'black',
                      # 调整节点中标签的字体颜色
                      negDashed = T,
                      #负相关的边为虚线
                      legend=T,
                      nodeNames = items,
                      legend.cex = 0.5,#图例的字体大小
                      legend.mode = 'style1',#图例的风格，默认是'style1'
                      pie = pred_obj$error[,2])
  }
  
  return(list(Network, p.network))
}


items <-list("Hypertension","Diabetes","Cancer","Lung disease","Heart disease",
             "Stroke","Psychological disorder","Arthritis","Cognitive disorder","Disability")
groups<-list("Chronic Diseases"=c(1:9),"Outcome" = 10)


p.network.elsa <- generate_network(df_net) 

p.network.charls <- generate_network(df_net.charls,set.layout = T) 

p.network.hrs <- generate_network(df_net.hrs,set.layout = T) 

p.network.mhas <- generate_network(df_net.mhas,set.layout = T) 
p.network.share <- generate_network(df_net.share,set.layout = T) 

pdf("./results/network/1charls.pdf", width = 8,height = 5)
plot(p.network.charls[[2]])
dev.off()

pdf("./results/network/2elsa.pdf", width = 8,height = 5)
plot(p.network.elsa[[2]])
dev.off()

pdf("./results/network/3hrs.pdf", width = 8,height = 5)
plot(p.network.hrs[[2]])
dev.off()

pdf("./results/network/4mhas.pdf", width = 8,height = 5)
plot(p.network.mhas[[2]])
dev.off()

pdf("./results/network/5share.pdf", width = 8,height = 5)
plot(p.network.share[[2]])
dev.off()

pdf("./results/network/1_1charls.pdf", width = 8,height = 5)
centralityPlot(p.network.charls[[1]],
               labels = items,
               include = c("Strength","Closeness","Betweenness"))
dev.off()

pdf("./results/network/2_1elsa.pdf", width = 8,height = 5)
centralityPlot(p.network.elsa[[1]],
               labels = items,
               include = c("Strength","Closeness","Betweenness"))
dev.off()

pdf("./results/network/3_1hrs.pdf", width = 8,height = 5)
centralityPlot(p.network.hrs[[1]],
               labels = items,
               include = c("Strength","Closeness","Betweenness"))
dev.off()

pdf("./results/network/4_1mhas.pdf", width = 8,height = 5)
centralityPlot(p.network.mhas[[1]],
               labels = items,
               include = c("Strength","Closeness","Betweenness"))
dev.off()

pdf("./results/network/5_1share.pdf", width = 8,height = 5)
centralityPlot(p.network.share[[1]],
               labels = items,
               include = c("Strength","Closeness","Betweenness"))
dev.off()

summary(p.network.share[[1]])


#------敏感性分析-------

network_sensitive <- function(data,simn){
  set.seed(123)
  Result1 <- bootnet(data, 
                     statistics = c("Strength","Closeness","Betweenness","edge"), 
                     nBoots = 100, 
                     nCores = 6)
  
  p1 <- plot(Result1, order = "sample") + #检验节点或边线之间的差异
    theme(
      axis.text.x = element_text(size = 16,colour = "black"),  # 仅修改X轴刻度标签字体大小并加粗
      axis.text.y = element_text(size = 16,colour = "black"),  # 仅修改Y轴刻度标签字体大小并加粗
      axis.title = element_text(size = 16, face = "bold"),   # X和Y轴标题
      legend.text = element_text(size = 18),   # 图例文本
      strip.text = element_text(size = 16, face = "bold"),    # 分面标签文本
      text = element_text(size = 16)           # 其他文本（如标签）
    )
  
  p2 <- plot(Result1, "strength", plot = "difference") + #分析节点在强度上的差异检验
    theme(
      axis.text.x = element_text(size = 16,colour = "black"),  # 仅修改X轴刻度标签字体大小并加粗
      axis.text.y = element_text(size = 16,colour = "black"),  # 仅修改Y轴刻度标签字体大小并加粗
      axis.title = element_text(size = 16, face = "bold"),   # X和Y轴标题
      legend.text = element_text(size = 18),   # 图例文本
      strip.text = element_text(size = 16, face = "bold"),    # 分面标签文本
      text = element_text(size = 16)           # 其他文本（如标签）
    )
  
  p3 <- plot(Result1, "closeness", plot = "difference") + #分析节点在紧密性上的差异检验
    theme(
      axis.text.x = element_text(size = 16,colour = "black"),  # 仅修改X轴刻度标签字体大小并加粗
      axis.text.y = element_text(size = 16,colour = "black"),  # 仅修改Y轴刻度标签字体大小并加粗
      axis.title = element_text(size = 16, face = "bold"),   # X和Y轴标题
      legend.text = element_text(size = 18),   # 图例文本
      strip.text = element_text(size = 16, face = "bold"),    # 分面标签文本
      text = element_text(size = 16)           # 其他文本（如标签）
    )
  
  p4 <- plot(Result1, "betweenness", plot = "difference") + #分析节点在中介性上的差异检验
    theme(
      axis.text.x = element_text(size = 16,colour = "black"),  # 仅修改X轴刻度标签字体大小并加粗
      axis.text.y = element_text(size = 16,colour = "black"),  # 仅修改Y轴刻度标签字体大小并加粗
      axis.title = element_text(size = 16, face = "bold"),   # X和Y轴标题
      legend.text = element_text(size = 18),   # 图例文本
      strip.text = element_text(size = 16, face = "bold"),    # 分面标签文本
      text = element_text(size = 16)           # 其他文本（如标签）
    )
  
  p5 <- plot(Result1,"edge", plot = "difference", onlyNonZero = TRUE, order = "sample") + 
    theme(
      axis.text.x = element_text(size = 16,colour = "black"),  # 仅修改X轴刻度标签字体大小并加粗
      axis.text.y = element_text(size = 16,colour = "black"),  # 仅修改Y轴刻度标签字体大小并加粗
      axis.title = element_text(size = 16, face = "bold"),   # X和Y轴标题
      legend.text = element_text(size = 18),   # 图例文本
      strip.text = element_text(size = 16, face = "bold"),    # 分面标签文本
      text = element_text(size = 16)           # 其他文本（如标签）
    )
  
  
  Result2 <- bootnet(data, 
                     statistics = c("Strength","Closeness","Betweenness"),
                     nBoots = simn, nCores = 6, type = "case")
  
  p6 <-plot(Result2, statistics=c("Strength","Closeness","Betweenness")) + #节点在中心性指标的稳定性
    theme(
      axis.text.x = element_text(size = 16,colour = "black"),  # 仅修改X轴刻度标签字体大小并加粗
      axis.text.y = element_text(size = 16,colour = "black"),  # 仅修改Y轴刻度标签字体大小并加粗
      axis.title = element_text(size = 16, face = "bold"),   # X和Y轴标题
      legend.text = element_text(size = 18),   # 图例文本
      strip.text = element_text(size = 16, face = "bold"),    # 分面标签文本
      text = element_text(size = 16)           # 其他文本（如标签）
    )
  
  return(list(p1,p2,p3,p4,p5,p6))
  
}


p.sen.charls <- network_sensitive(p.network.charls[[1]],simn = 1000)
p.sen.elsa <- network_sensitive(p.network.elsa[[1]],simn = 1000)
p.sen.hrs <- network_sensitive(p.network.hrs[[1]],simn = 1000)
p.sen.mhas <- network_sensitive(p.network.mhas[[1]],simn = 1000)
p.sen.share <- network_sensitive(p.network.share[[1]],simn = 1000)

library(patchwork)
library(tidyverse)
pdf("./results/network/sen_charls2.pdf", width = 36,height = 24)
(p.sen.charls[[1]] +  p.sen.charls[[6]] + p.sen.charls[[2]]) /
  (p.sen.charls[[3]] + p.sen.charls[[4]] + p.sen.charls[[5]]) +
  plot_annotation(tag_levels = "A")  # 添加标签 A, B, C 等
dev.off()

pdf("./results/network/sen_elsa.pdf", width = 24,height = 16)
(p.sen.elsa[[1]] + p.sen.elsa[[6]] + p.sen.elsa[[2]]) /
  (p.sen.elsa[[3]] + p.sen.elsa[[4]] + p.sen.elsa[[5]]) +
  plot_annotation(tag_levels = "A")  # 添加标签 A, B, C 等
dev.off()

pdf("./results/network/sen_hrs.pdf", width = 24,height = 16)
(p.sen.hrs[[1]] + p.sen.hrs[[6]] + p.sen.hrs[[2]]) /
  (p.sen.hrs[[3]] + p.sen.hrs[[4]] + p.sen.hrs[[5]]) +
  plot_annotation(tag_levels = "A")  # 添加标签 A, B, C 等
dev.off()

pdf("./results/network/sen_mhas.pdf", width = 24,height = 16)
(p.sen.mhas[[1]] + p.sen.mhas[[6]] + p.sen.mhas[[2]]) /
  (p.sen.mhas[[3]] + p.sen.mhas[[4]] + p.sen.mhas[[5]]) +
  plot_annotation(tag_levels = "A")  # 添加标签 A, B, C 等
dev.off()

pdf("./results/network/sen_share.pdf", width = 24,height = 16)
(p.sen.share[[1]] + p.sen.share[[6]] + p.sen.share[[2]]) /
  (p.sen.share[[3]] + p.sen.share[[4]] + p.sen.share[[5]]) +
  plot_annotation(tag_levels = "A")  # 添加标签 A, B, C 等
dev.off()



save(p.sen.charls,
     p.sen.elsa,
     p.sen.hrs,
     p.sen.mhas,
     p.sen.share,
     file = "./data/network_sensitive_results.RData")

load("./data/network_sensitive_results.RData")