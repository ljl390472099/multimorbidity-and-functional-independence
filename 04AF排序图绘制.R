
#-----计算AF-----
library(AF)
library(patchwork)
source("./function/extract_RD.R")
source("./function/tidy_RD.R")
source("./function/run_AF_analysis.R")

get_af <- function(data,f){
  library(dplyr)
  library(tidyr)
  
  # 检查输入数据
  if (!"detailed_category" %in% colnames(data)) {
    stop("数据中缺少列 'detailed_category'")
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
  
  # 动态重命名变量，仅重命名存在的变量
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
  # 检查并执行变量重命名
  existing_renames <- rename_map[names(rename_map) %in% colnames(data)]
  data <- data %>% rename(!!!setNames(names(existing_renames), existing_renames))
  
  # 检查并移除不存在的变量
  f_vars <- all.vars(f)
  existing_vars <- f_vars[f_vars %in% colnames(data)]
  f <- reformulate(existing_vars[-1], response = existing_vars[1])
  # 定义慢性疾病变量列表
  chronic_diseases <- intersect(c("mental2", "Diabetes", "physical_mental",
                                  "mental_cog", "Hypertension", "Arthritis",      
                                  "Stroke", "lung2", "heart2",          
                                  "Cancer", "phy_cog", "cognition",      
                                  "phy_cog_mental"), colnames(data))
  fit.log <- glm(f,
                 family = "binomial",
                 data = data)
  
  # 创建一个空数据框用于存储结果
  results_df <- data.frame(Disease = character(),
                           AF = numeric(),
                           Std.Error = numeric(),
                           z_value = numeric(),
                           Pr_z = numeric(),
                           stringsAsFactors = FALSE)
  
  # 遍历每个慢性疾病变量
  for (disease in chronic_diseases) {
    # 执行 AF.cs 计算
    AF_result <- AFglm(object = fit.log, 
                       data = data, 
                       exposure = disease)
    
    # 获取摘要结果
    AF_summary <- summary(AF_result)
    
    # 提取所需指标
    disease_result <- data.frame(Disease = disease,
                                 AF_summary$AF)
    
    # 将结果行追加到数据框
    results_df <- rbind(results_df, disease_result)
  } 
  return(results_df)
}




f1 <- adl ~ age + sex + marr + edu + exercise + smoke + drink + income.c+ bmi.c+
  Hypertension+ Diabetes+Cancer+lung2+heart2+Stroke+Arthritis+mental2+cognition+
  physical_mental+phy_cog+mental_cog+phy_cog_mental

f2 <- adl ~   sex + marr + edu + exercise + smoke + drink + income.c+ bmi.c+
  Hypertension+ Diabetes+Cancer+lung2+heart2+Stroke+Arthritis+mental2+cognition+
  physical_mental+phy_cog+mental_cog+phy_cog_mental

f3 <- adl ~ age + marr + edu + exercise + smoke + drink + income.c+ bmi.c+
  Hypertension+ Diabetes+Cancer+lung2+heart2+Stroke+Arthritis+mental2+cognition+
  physical_mental+phy_cog+mental_cog+phy_cog_mental

af.result <- AFglm(object = fit.log, data = data, exposure = "Hypertension")
summary(af.result)

df.af <- get_af(charls.long,f=f1) |> 
  mutate(group = "Total") |> 
  rbind(get_af(charls.long[charls.long$age.c == "<60",],f=f2) |> 
          mutate(group = "<60")) |> 
  rbind(get_af(charls.long[charls.long$age.c == "60-74",],f=f2) |> 
          mutate(group = "60-74")) |> 
  rbind(get_af(charls.long[charls.long$age.c == "75+",],f=f2) |> 
          mutate(group = "75+")) |>
  rbind(get_af(charls.long[charls.long$sex == "Male",],f=f3) |> 
          mutate(group = "Male")) |>
  rbind(get_af(charls.long[charls.long$sex == "Female",],f=f3) |> 
          mutate(group = "Female")) |> 
  mutate(cohort = "CHARLS (China)") |> 
  rbind(get_af(elsa.long,f=f1) |> 
          mutate(group = "Total") |> 
          rbind(get_af(elsa.long[elsa.long$age.c == "<60",],f=f2) |> 
                  mutate(group = "<60")) |> 
          rbind(get_af(elsa.long[elsa.long$age.c == "60-74",],f=f2) |> 
                  mutate(group = "60-74")) |> 
          rbind(get_af(elsa.long[elsa.long$age.c == "75+",],f=f2) |> 
                  mutate(group = "75+")) |>
          rbind(get_af(elsa.long[elsa.long$sex == "Male",],f=f3) |> 
                  mutate(group = "Male")) |>
          rbind(get_af(elsa.long[elsa.long$sex == "Female",],f=f3) |> 
                  mutate(group = "Female")) |> 
          mutate(cohort = "ELSA (UK)")) |>
  rbind(get_af(hrs.long,f=f1) |> 
          mutate(group = "Total") |> 
          rbind(get_af(hrs.long[hrs.long$age.c == "<60",],f=f2) |> 
                  mutate(group = "<60")) |> 
          rbind(get_af(hrs.long[hrs.long$age.c == "60-74",],f=f2) |> 
                  mutate(group = "60-74")) |> 
          rbind(get_af(hrs.long[hrs.long$age.c == "75+",],f=f2) |> 
                  mutate(group = "75+")) |>
          rbind(get_af(hrs.long[hrs.long$sex == "Male",],f=f3) |> 
                  mutate(group = "Male")) |>
          rbind(get_af(hrs.long[hrs.long$sex == "Female",],f=f3) |> 
                  mutate(group = "Female")) |> 
          mutate(cohort = "HRS (USA)")) |> 
  rbind(get_af(mhas.long,f=f1) |> 
          mutate(group = "Total") |> 
          rbind(get_af(mhas.long[mhas.long$age.c == "<60",],f=f2) |> 
                  mutate(group = "<60")) |> 
          rbind(get_af(mhas.long[mhas.long$age.c == "60-74",],f=f2) |> 
                  mutate(group = "60-74")) |> 
          rbind(get_af(mhas.long[mhas.long$age.c == "75+",],f=f2) |> 
                  mutate(group = "75+")) |>
          rbind(get_af(mhas.long[mhas.long$sex == "Male",],f=f3) |> 
                  mutate(group = "Male")) |>
          rbind(get_af(mhas.long[mhas.long$sex == "Female",],f=f3) |> 
                  mutate(group = "Female")) |> 
          mutate(cohort = "MHAS(Mexico)")) |> 
  rbind(get_af(share.long,f=f1) |> 
          mutate(group = "Total") |> 
          rbind(get_af(share.long[share.long$age.c == "<60",],f=f2) |> 
                  mutate(group = "<60")) |> 
          rbind(get_af(share.long[share.long$age.c == "60-74",],f=f2) |> 
                  mutate(group = "60-74")) |> 
          rbind(get_af(share.long[share.long$age.c == "75+",],f=f2) |> 
                  mutate(group = "75+")) |>
          rbind(get_af(share.long[share.long$sex == "Male",],f=f3) |> 
                  mutate(group = "Male")) |>
          rbind(get_af(share.long[share.long$sex == "Female",],f=f3) |> 
                  mutate(group = "Female")) |> 
          mutate(cohort = "SHARE(Europe)"))

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

dis_name2 <- dis_name |>
  as.data.frame() |> 
  mutate(dis_name = case_when(dis_name == "No disease" ~ "No disorders",
                              dis_name == "Chronic lung disease" ~ "Lung disease",
                              dis_name == "Mental disease only" ~ "Psychological disorder",
                              dis_name == "Cognitive impairment only" ~ "Cognitive disorder",
                              dis_name == "Physical-Mental comorbidity" ~ "Physical–psychological multimorbidity",
                              dis_name == "Physical-Cognitive comorbidity" ~ "Physical–cognitive multimorbidity",
                              dis_name == "Mental-Cognitive comorbidity" ~ "Psychological–cognitive multimorbidity",
                              dis_name == "Physical-Mental-Cognitive comorbidity" ~ "Physical–psychological–cognitive multimorbidity",
                              TRUE ~ dis_name)) |> 
  pull()
  
# 创建映射关系
dis_match <- data.frame(
  disease = disease,
  dis_name = dis_name) 
names(df.af) <-c("disease","af","se","z","p","group","cohort")
df.af <- left_join(df.af,dis_match,by="disease")|> 
  mutate(dis_name = case_when(dis_name == "No disease" ~ "No disorders",
                              dis_name == "Chronic lung disease" ~ "Lung disease",
                              dis_name == "Mental disease only" ~ "Psychological disorder",
                              dis_name == "Cognitive impairment only" ~ "Cognitive disorder",
                              dis_name == "Physical-Mental comorbidity" ~ "Physical–psychological multimorbidity",
                              dis_name == "Physical-Cognitive comorbidity" ~ "Physical–cognitive multimorbidity",
                              dis_name == "Mental-Cognitive comorbidity" ~ "Psychological–cognitive multimorbidity",
                              dis_name == "Physical-Mental-Cognitive comorbidity" ~ "Physical–psychological–cognitive multimorbidity",
                              TRUE ~ dis_name))


names(df.af) <-c("disease","af","se","z","p","group","cohort","dis_name")

df.af <- df.af |> 
  mutate(lci = af-1.96*se,
         uci = af+1.96*se)

write.csv(df.af,"./results/df_af.csv",row.names = F)

# 绘制分组af图
ggplot(df.af , aes(x = reorder(dis_name, af), y = af*100, fill = cohort)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ cohort + group , nrow = 5) +
  coord_flip() +
  labs(x = "Disease", y = "Attributable Fraction (%)", title = "Attributable Fraction by Disease") +
  theme_minimal()  


#-----绘制AF排序图:按年龄分组 -----
library(tidyverse)
library(ggtext)
library(magrittr)
library(patchwork)

plot_af_order <- function(cohort1){
  pastel_nature_colors <- data.frame(dis_name=dis_name2,
                                     colour = c(
                                       "#8FDBF3","#A2C4F1","#CEBAF0","#C6C3E1",
                                       "#EEC2E5","#FFCFD1","#FBE3CD","#F9F9CA",
                                       "#FAA09C","#F9B29C","#F9C89B","#FBDF9D","#E9E4AF"))
  
  df <- df.af |> 
    group_by(group,cohort) %>%
    arrange(desc(af), .by_group = TRUE,) %>%
    mutate(id = row_number()) %>%
    ungroup() |> 
    mutate(group = ifelse(group == "45-59", "<60", group)) |> 
    mutate(rei_2 = paste(id," ",dis_name),
           g = case_when(group == "<60" ~ 1,
                         group == "60-74" ~ 2,
                         group == "75+" ~ 3,
                         group == "Total" ~ 4)) |> 
    mutate(id=14-id) |> 
    filter(cohort == cohort1 & group %in% c("<60","60-74","75+","Total")) |> 
    left_join(pastel_nature_colors,by=c("dis_name"))
  
 
  
  colour1 <- df |> filter(g==1) |> dplyr::pull(colour)
  colour2 <- df |> filter(g==2) |> dplyr::pull(colour)
  colour3 <- df |> filter(g==3) |> dplyr::pull(colour)
  colour4 <- df |> filter(g==4) |> dplyr::pull(colour)
  
  
  p <- df %>% 
    filter(cohort == cohort1 & g %in% c("<60","60-74","75+","Total")) |> 
    ggplot(aes(x=group, y=id))+
    # 绘制第一列边框
    geom_rect(data=df %>% filter(g==1),
              aes(xmin=g, xmax=g-2,ymin=id-0.3,ymax=id+0.3),
              fill= colour1,color="black")+
    # 添加第一列数据文本
    geom_text(data=df %>% filter(g==1),
              aes(x=g-1,label=rei_2,y=id),size=3.5,color="black")+
    # 绘制第2列数据框
    geom_rect(data=df |> filter(g==1),aes(xmin=1.02,xmax=1.42,ymin=id-0.3,ymax=id+0.3),
              fill=colour1,color="black")+
    geom_text(data=df |> filter(g==1),aes(x=1.2,label=round(af*100,digits=2),
                                          y=id),size=3.5,color="black")+
    scale_y_continuous(limits=c(0,14))+
    # 绘制第3列边框
    geom_rect(data=df %>% filter(g==2),
              aes(xmin=g,xmax=g+2,ymin=id-0.3,ymax=id+0.3),
              fill=colour2,color="black")+
    # 添加第3列数据文本
    geom_text(data=df %>% filter(g==2),
              aes(x=g+1,label=rei_2,y=id),size=3.5,color="black")+
    # 添加连接线
    geom_line(data=df |> filter(g ==1|g==2),aes(x=ifelse(g==1,g+0.42,g),
                                                y=id,group=dis_name),color="#5785C1")+
    # 绘制第4列数据
    geom_rect(data=df %>% filter(g==2),aes(xmin=4.02, xmax=4.4,ymin=id-0.3,ymax=id+0.3),
              fill=colour2,color="black")+
    geom_text(data=df |> filter(g==2),aes(x=4.2,label=round(af*100,digits=2),
                                          y=id),size=3.5,color="black")+
    # 绘制第5列边框
    geom_rect(data=df %>% filter(g==3),
              aes(xmin=g+2,xmax=g+4,ymin=id-0.3,ymax=id+0.3),
              fill=colour3,color="black")+
    # 添加第5列数据文本
    geom_text(data=df %>% filter(g==3),
              aes(x=g+3,label=rei_2,y=id),size=3.5,color="black")+
    # 添加连接线
    geom_line(data=df |> filter(g ==3|g==2),aes(x=ifelse(g==2,4.4,g+2),
                                                y=id,group=dis_name),color="#5785C1")+
    # 绘制第6列数据
    geom_rect(data= df %>% filter(g==3),aes(xmin=7.02, xmax=7.4,ymin=id-0.3,ymax=id+0.3),
              fill=colour3,color="black")+
    geom_text(data=df |> filter(g==3),aes(x=7.2,label=round(af*100,digits=2),
                                          y=id),size=3.5,color="black")+
    # 绘制第7列
    geom_rect(data=df %>% filter(g==4),
              aes(xmin=g+4,xmax=g+6,ymin=id-0.3,ymax=id+0.3),
              fill=colour4,color="black")+
    # 添加第7列数据文本
    geom_text(data=df %>% filter(g==4),
              aes(x=g+5,label=rei_2,y=id),size=3.5,color="black")+
    # 添加连接线
    geom_line(data=df |> filter(g ==3|g==4),aes(x=ifelse(g==3,7.4,g+4),
                                                y=id,group=dis_name),color="#5785C1")+
    # 绘制第8列数据
    geom_rect(data=df |> filter(g==4),aes(xmin=10.02, xmax=10.4,ymin=id-0.3,ymax=id+0.3),
              fill=colour4,color="black")+
    geom_text(data=df |> filter(g==4),aes(x=10.2,label=round(af*100,digits=2),
                                          y=id),size=3.5,color="black")+
    scale_fill_manual(values = colour1) +
    
    # 添加注释文本内容
    annotate(geom="text",y=14,x=0, label="Age < 60",size=4.3, fontface="bold") +
    annotate(geom="text",y=14,x=3, label="Age 60-74",size=4.3, fontface="bold") +
    annotate(geom="text",y=14,x=6, label="Age 75+",size=4.3, fontface="bold") +
    annotate(geom="text",y=14,x=9, label="Total",size=4.3, fontface="bold") +
    theme_void()
  
  return(p)
  
}

p.share <- plot_af_order(cohort1 = "SHARE(Europe)")
p.elsa <- plot_af_order(cohort1 = "ELSA (UK)")
p.charls <- plot_af_order(cohort1 = "CHARLS (China)")
p.hrs <- plot_af_order(cohort1 = "HRS (USA)")
p.mhas <- plot_af_order(cohort1 = "MHAS(Mexico)")

pdf("./results/network/af_order.pdf", width = 20,height = 25)
p.charls / p.elsa / p.hrs / p.mhas / p.share + 
  plot_annotation(tag_levels="A")# 添加标签“A”，“B”
dev.off()

#-----绘制AF排序图：按性别分--------

plot_af_sex <- function(cohort1){
  pastel_nature_colors <- data.frame(dis_name=dis_name2,
                                     colour = c(
                                       "#8FDBF3","#A2C4F1","#CEBAF0","#C6C3E1",
                                       "#EEC2E5","#FFCFD1","#FBE3CD","#F9F9CA",
                                       "#FAA09C","#F9B29C","#F9C89B","#FBDF9D","#E9E4AF"))
  
  df <- df.af |> 
    group_by(group,cohort) %>%
    arrange(desc(af), .by_group = TRUE,) %>%
    mutate(id = row_number()) %>%
    ungroup() |> 
    mutate(group = ifelse(group == "45-59", "<60", group)) |> 
    mutate(rei_2 = paste(id," ",dis_name),
           g = case_when(group == "Female" ~ 1,
                         group == "Male" ~ 2,
                         group == "Total" ~ 3)) |> 
    mutate(id=14-id) |> 
    filter(cohort == cohort1 & group %in% c("Female","Male","Total")) |> 
    left_join(pastel_nature_colors,by=c("dis_name"))
  
  
  
  colour1 <- df |> filter(g==1) |> dplyr::pull(colour)
  colour2 <- df |> filter(g==2) |> dplyr::pull(colour)
  colour3 <- df |> filter(g==3) |> dplyr::pull(colour)
 
  
  p <- df %>% 
    filter(cohort == cohort1 & g %in% c("Female","Male","Total")) |> 
    ggplot(aes(x=group, y=id))+
    # 绘制第一列边框
    geom_rect(data=df %>% filter(g==1),
              aes(xmin=g, xmax=g-2,ymin=id-0.3,ymax=id+0.3),
              fill= colour1,color="black")+
    # 添加第一列数据文本
    geom_text(data=df %>% filter(g==1),
              aes(x=g-1,label=rei_2,y=id),size=3.5,color="black")+
    # 绘制第2列数据框
    geom_rect(data=df |> filter(g==1),aes(xmin=1.02,xmax=1.42,ymin=id-0.3,ymax=id+0.3),
              fill=colour1,color="black")+
    geom_text(data=df |> filter(g==1),aes(x=1.2,label=round(af*100,digits=2),
                                          y=id),size=3.5,color="black")+
    scale_y_continuous(limits=c(0,14))+
    # 绘制第3列边框
    geom_rect(data=df %>% filter(g==2),
              aes(xmin=g,xmax=g+2,ymin=id-0.3,ymax=id+0.3),
              fill=colour2,color="black")+
    # 添加第3列数据文本
    geom_text(data=df %>% filter(g==2),
              aes(x=g+1,label=rei_2,y=id),size=3.5,color="black")+
    # 添加连接线
    geom_line(data=df |> filter(g ==1|g==2),aes(x=ifelse(g==1,g+0.42,g),
                                                y=id,group=dis_name),color="#5785C1")+
    # 绘制第4列数据
    geom_rect(data=df %>% filter(g==2),aes(xmin=4.02, xmax=4.4,ymin=id-0.3,ymax=id+0.3),
              fill=colour2,color="black")+
    geom_text(data=df |> filter(g==2),aes(x=4.2,label=round(af*100,digits=2),
                                          y=id),size=3.5,color="black")+
    # 绘制第5列边框
    geom_rect(data=df %>% filter(g==3),
              aes(xmin=g+2,xmax=g+4,ymin=id-0.3,ymax=id+0.3),
              fill=colour3,color="black")+
    # 添加第5列数据文本
    geom_text(data=df %>% filter(g==3),
              aes(x=g+3,label=rei_2,y=id),size=3.5,color="black")+
    # 添加连接线
    geom_line(data=df |> filter(g ==3|g==2),aes(x=ifelse(g==2,4.4,g+2),
                                                y=id,group=dis_name),color="#5785C1")+
    # 绘制第6列数据
    geom_rect(data= df %>% filter(g==3),aes(xmin=7.02, xmax=7.4,ymin=id-0.3,ymax=id+0.3),
              fill=colour3,color="black")+
    geom_text(data=df |> filter(g==3),aes(x=7.2,label=round(af*100,digits=2),
                                          y=id),size=3.5,color="black")+
    scale_fill_manual(values = colour1) +
    
    # 添加注释文本内容
    annotate(geom="text",y=14,x=0, label="Female",size=4.3, fontface="bold") +
    annotate(geom="text",y=14,x=3, label="Male",size=4.3, fontface="bold") +
    annotate(geom="text",y=14,x=6, label="Total",size=4.3, fontface="bold") +
   # annotate(geom="text",y=14,x=9, label="Total",size=4.3, fontface="bold") +
    theme_void()
  
  return(p)
  
}

p.share.sex <- plot_af_sex(cohort1 = "SHARE(Europe)")
p.elsa.sex <- plot_af_sex(cohort1 = "ELSA (UK)")
p.charls.sex <- plot_af_sex(cohort1 = "CHARLS (China)")
p.hrs.sex <- plot_af_sex(cohort1 = "HRS (USA)")
p.mhas.sex <- plot_af_sex(cohort1 = "MHAS(Mexico)")

pdf("./results/network/af_order_sex.pdf", width = 15,height = 25)
p.charls.sex / p.elsa.sex / p.hrs.sex / p.mhas.sex / p.share.sex + 
  plot_annotation(tag_levels="A")# 添加标签“A”，“B”
dev.off()

