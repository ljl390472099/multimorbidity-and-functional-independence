# 03 RD热图绘制
library(patchwork)
#------绘制排序热图------

loc_rank <- df.rd |> 
  dplyr::select(contrast,rank,group,cohort) |> 
  mutate(group = paste0(cohort," - ",group)) |> 
  dplyr::select(contrast,rank,group) |> 
  mutate(contrast = str_remove(contrast, " - No disease") )

loc_rank1 <- df.rd |> 
  dplyr::select(contrast,rank,group,cohort) |> 
  mutate(group = paste0(cohort," - ",group)) |> 
  dplyr::select(contrast,rank,group) |> 
  spread(key = group,value = rank) |> 
  mutate(contrast = str_remove(contrast, " - No disease")) 
loc_rank1$cause_rank <- rowSums(loc_rank1[,2:ncol(loc_rank1)],na.rm = T)
loc_rank1 <- loc_rank1[order(loc_rank1$cause_rank),]

loc_rank$rank <- factor(loc_rank$rank,levels = sort(unique(loc_rank$rank)))
loc_rank$contrast <- factor(loc_rank$contrast,levels = rev(loc_rank1$contrast))

rank_color <- c("#d73027",#1                
                "#fb8c59",#2                
                "#fddf8f","#fddf8f",#3-4                
                "#e0f3f8","#e0f3f8","#e0f3f8","#e0f3f8",#5-8                
                "#91bfdb","#91bfdb","#91bfdb","#91bfdb","#91bfdb",#9-13                
                "#4575b4","#4575b4","#4575b4","#4575b4","#4575b4","#4575b4","#4575b4","#4575b4","#4575b4","#4575b4")#13+
# 颜色与排名对应
names(rank_color) <- sort(unique(loc_rank$rank))

loc_rank <- loc_rank |> 
  mutate(contrast = case_when(contrast == "No disease" ~ "No disorders",
                              contrast == "Chronic lung disease" ~ "Lung disease",
                              contrast == "Mental disease only" ~ "Psychological disorder",
                              contrast == "Cognitive impairment only" ~ "Cognitive disorder",
                              contrast == "Physical-Mental comorbidity" ~ "Physical–psychological multimorbidity",
                              contrast == "Physical-Cognitive comorbidity" ~ "Physical–cognitive multimorbidity",
                              contrast == "Mental-Cognitive comorbidity" ~ "Psychological–cognitive multimorbidity",
                              contrast == "Physical-Mental-Cognitive comorbidity" ~ "Physical–psychological–cognitive multimorbidity",
                              TRUE ~ contrast))



p1 <- ggplot(loc_rank,aes(x=group,y=contrast))+  
  geom_tile(aes(fill=rank),colour="white",size=0.3)+  # scale_fill_brewer(palette ="RdBu")+  
  geom_text(aes(label=rank),family="mono",colour="black")+  
  scale_x_discrete(limits = unique(loc_rank$group))+  
  scale_y_discrete(limits = rev(unique(loc_rank1$contrast)))+  
  scale_fill_manual(values = rank_color)+  xlab("")+ylab("")+  
  theme(    
    axis.text.x=element_text(angle = 45,colour="black",size=10,
                             hjust = 1,vjust = 1,family = "Arial"),    
    axis.text.y=element_text(colour="black",size=10,family = "Arial"),    
    panel.background = element_blank(),    
    panel.grid = element_blank()  )+  
  guides(fill = FALSE)


#-----按队列绘图------
df.rd2 <- df.rd |> 
  mutate(contrast = str_remove(contrast, " - No disease")) |> 
  mutate(contrast = case_when(contrast == "No disease" ~ "No disorders",
                              contrast == "Chronic lung disease" ~ "Lung disease",
                              contrast == "Mental disease only" ~ "Psychological disorder",
                              contrast == "Cognitive impairment only" ~ "Cognitive disorder",
                              contrast == "Physical-Mental comorbidity" ~ "Physical–psychological multimorbidity",
                              contrast == "Physical-Cognitive comorbidity" ~ "Physical–cognitive multimorbidity",
                              contrast == "Mental-Cognitive comorbidity" ~ "Psychological–cognitive multimorbidity",
                              contrast == "Physical-Mental-Cognitive comorbidity" ~ "Physical–psychological–cognitive multimorbidity",
                              TRUE ~ contrast)) |> 
  filter(group != "Total")

p2 <- ggplot(df.rd2, aes(x = reorder(contrast, RD), y = RD * 100, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  facet_wrap(~ cohort + group, nrow = 5, strip.position = "top") +
  coord_flip() +
  labs(
    x = "Disease",
    y = "RD (%)"
  ) +
  scale_fill_brewer(palette = "Set2") + # 使用更柔和的配色
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(), # 去除网格线
    panel.border = element_blank(), # 去除面板边框
    strip.background = element_blank(), # 去除分面标题背景
    strip.text = element_text(size = 12, face = "bold"), # 分面标题样式
    axis.text = element_text(color = "black", size = 10), # 坐标轴文字样式
    axis.title = element_text(size = 12, face = "bold"), # 坐标轴标题样式
    legend.position = "right", # 调整图例位置
    legend.title = element_blank(), # 去除图例标题
    legend.text = element_text(size = 10)
  )

p1 / p2 + 
  plot_layout(heights = c(1, 2))  + # 调整 p1 和 p2 的相对高度
  plot_annotation(tag_levels = "A") 

#-----绘制气泡图------# 绘制按年龄段分的气泡图
library(tidyverse)
library(ggh4x)
library(scales)
# remotes::install_github("EvaMaeRey/ggcirclepack")
library(ggcirclepack)


library(scales)  # 用于调整透明度

library(ggpackets)  # 用于 geom_circlepack
plot_rd <- function(cohort_name){
  df.rd |> 
    mutate(contrast = str_remove(contrast, " - No disease")) |> 
    mutate(contrast = case_when(
      contrast == "Psychological disorder" ~ "Psychological disorder",
      contrast == "Diabetes" ~ "Diabetes",
      contrast == "Physical–psychological multimorbidity" ~ "Physical–psychological multimorbidity",
      contrast == "Psychological–cognitive multimorbidity" ~ "Psychological–cognitive multimorbidity",
      contrast == "Hypertension" ~ "Hypertension",
      contrast == "Arthritis" ~ "Arthritis",
      contrast == "Stroke" ~ "Stroke",
      contrast == "Lung disease" ~ "Lung disease",
      contrast == "Heart disease" ~ "Heart disease",
      contrast == "Cancer" ~ "Cancer",
      contrast == "Physical–cognitive multimorbidity" ~ "Physical–cognitive multimorbidity",
      contrast == "Cognitive disorder" ~ "Cognitive disorder",
      contrast == "Physical–psychological–cognitive multimorbidity" ~ "Physical–psychological–cognitive multimorbidity",
      TRUE ~ contrast)) |> 
    filter(cohort == cohort_name) |>
    filter(group != "Total"  & RD > 0) |>  # 筛选所需数据
    group_by(group) |> 
    arrange(desc(RD), .by_group = TRUE) |>  # 根据每个分面 (group) 的 RD 值降序排列
    mutate(contrast = factor(contrast, levels = unique(contrast))) |>  # 转换为有序因子
    ggplot(aes(id = contrast, area = RD, fill = contrast)) +  # 用面积表示 value
    geom_circlepack() +
    geom_circlepack_text(aes(label = round(RD, 2)), size = 3) +  # 添加数值标签
    facet_nested(~ group) +
    scale_size_continuous(range = c(0, 1)) +
    scale_fill_manual(
      values = c(
        "Psychological disorder" = alpha("#88CCEE", 0.5),  # 蓝色，80% 不透明度
        "Diabetes" = alpha("#CC6677", 0.5),  # 紫红色，80% 不透明度
        "Physical–psychological multimorbidity" = alpha("#DDCC77", 0.5),  # 浅黄色，80% 不透明度
        "Psychological–cognitive multimorbidity" = alpha("#117733", 0.5),  # 深蓝绿，80% 不透明度
        "Hypertension" = alpha("#332288", 0.5),  # 深蓝色，80% 不透明度
        "Arthritis" = alpha("#AA4499", 0.5),  # 紫色，80% 不透明度
        "Stroke" = alpha("#44AA99", 0.5),  # 蓝绿，80% 不透明度
        "Lung disease" = alpha("#999933", 0.5),  # 橄榄色，80% 不透明度
        "Heart disease" = alpha("#882255", 0.5),  # 深紫红，80% 不透明度
        "Cancer" = alpha("#661100", 0.5),  # 深棕红，80% 不透明度
        "Physical–cognitive multimorbidity" = alpha("#6699CC", 0.5),  # 浅蓝灰，80% 不透明度
        "Cognitive disorder" = alpha("#888888", 0.5),  # 灰色，80% 不透明度
        "Physical–psychological–cognitive multimorbidity" = alpha("#E69F00", 0.5)  # 橙色，80% 不透明度
      )) +
    theme(
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm"),
      strip.background = element_rect(fill = "grey90", color = "white"),
      strip.text = element_text(color = "black", face = "bold"),
      panel.spacing.x = unit(0, "cm"),
      panel.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    guides(
      fill = guide_legend(title = NULL), 
      size = guide_legend(title = NULL)
    )  # 去除图例标题
  
}

plot_rd("CHARLS") / plot_rd("ELSA") / plot_rd("HRS") / plot_rd("MHAS") / plot_rd("SHARE")+ 
  
  plot_annotation(tag_levels = "A") 


#-----环状棒棒糖图绘制------
library(tidyverse)
library(circlize)

pastel_nature_colors <- data.frame(dis_name=sort(dis_name2),
                                   color = c(
                                     "#8FDBF3","#A2C4F1","#CEBAF0","#C6C3E1",
                                     "#EEC2E5","#FFCFD1","#fbd37c","#f5c75b",
                                     "#FAA09C","#F9B29C","#F9C89B","#FBDF9D","#E9E4AF"),
                                   label = paste0("D", 1:length(dis_name2)))
df.rd3 <-  df.rd |> 
  mutate(contrast = str_remove(contrast, " - No disease")) |> 
  mutate(contrast = case_when(
    contrast == "Psychological disorder" ~ "Psychological disorder",
    contrast == "Diabetes" ~ "Diabetes",
    contrast == "Physical–psychological multimorbidity" ~ "Physical–psychological multimorbidity",
    contrast == "Psychological–cognitive multimorbidity" ~ "Psychological–cognitive multimorbidity",
    contrast == "Hypertension" ~ "Hypertension",
    contrast == "Arthritis" ~ "Arthritis",
    contrast == "Stroke" ~ "Stroke",
    contrast == "Lung disease" ~ "Lung disease",
    contrast == "Heart disease" ~ "Heart disease",
    contrast == "Cancer" ~ "Cancer",
    contrast == "Physical–cognitive multimorbidity" ~ "Physical–cognitive multimorbidity",
    contrast == "Cognitive disorder" ~ "Cognitive disorder",
    contrast == "Physical–psychological–cognitive multimorbidity" ~ "Physical–psychological–cognitive multimorbidity",
    TRUE ~ contrast)) |> 
  filter(cohort == "CHARLS") |>
  filter(group != "Total"  & RD > 0) |> 
  group_by(group) %>%
  arrange(desc(RD), .by_group = TRUE) %>%
  mutate(x_pos = (seq_along(RD) - 0.5) / n()) %>%
  mutate(y_pos = RD / max(RD)) %>%
  ungroup() %>%
  left_join(
    pastel_nature_colors,
    by = c("contrast" = "dis_name")
  )

# 初始化
# circlizecircos.clear()
circos.par(start.degree =90,
           gap.after = c(rep(2, n_distinct(df.rd3$group) - 1),20),
           track.margin = c(0, 0.01),
           cell.padding = c(0, 0, 0,0))
circos.initialize(factors = df.rd3$group, 
                  xlim = c(0, 1))

# 提取每个组的前5个基因，用于添加标签名
label_data <- df.rd3 %>%
  group_by(group) %>%
  # slice_head(n = 5) %>%
  mutate(start = x_pos, end = x_pos) %>%
  ungroup() %>%
  select(group, start, end, label, color)

# 添加基因名称标签
circos.genomicLabels(label_data,
                   labels.column = 4,
                     facing = "reverse.clockwise",
                     side = "outside",
                     cex = 0.7, 
                    col = "black",
                     connection_height = convert_height(4,"mm"),
                     line_lwd = 0.5)
# 添加扇形分类标签
circos.track( factors = df.rd3$group, 
              ylim = c(0, 1),
              track.height = 0.08,
              bg.col = c("#f8bf99","#fda453","#f47720","#87bdc0","#5fadac"),
              panel.fun = function(x, y) {
                circos.text(
                  CELL_META$xcenter,
                  CELL_META$ylim[2]-0.7,
                  CELL_META$sector.index,
                  facing = "bending.inside",
                  cex = 0.8,
                  adj = c(0.5,0)
                  ) })
#绘制棒棒图
circos.trackPlotRegion(
  factors = df.rd3$group, 
  ylim = c(0, 1), 
  track.height = 0.37, 
  panel.fun = function(x, y) {
    sector_data <- df.rd3[df.rd3$group == CELL_META$sector.index, ]
    
    for (i in 1:nrow(sector_data)) {
      x_pos <- sector_data$x_pos[i]
      y_pos <- sector_data$y_pos[i]
      
      # 添加线条
      circos.lines(
        c(x_pos, x_pos), 
        c(0, y_pos), 
        col = sector_data$color[i], 
        lwd = 3
      )
      
      # 添加点
      circos.points(
        x_pos, 
        y_pos, 
        pch = 16, 
        col = sector_data$color[i], 
        cex = 1.5
      )
    }
    
    # 在第一个扇区添加y轴刻度线
    if (CELL_META$sector.index == "<60") {
      circos.yaxis(
        side = "left", 
        at = seq(0, 1, by = 0.2), 
        labels.cex = 0.8,
        tick.length = convert_length(2, "mm")
      )
      
      # 添加y轴标题
      circos.text(
        CELL_META$cell.xlim[1] - convert_x(10, "mm"), 
        0.5, 
        "RD", 
        facing = "reverse.clockwise", 
        niceFacing = TRUE, 
        adj = c(0.5, 0.5), 
        cex = 0.8
      )
    }
  }
)

#------编写环形棒棒图函数--------
pastel_nature_colors <- data.frame(dis_name=sort(dis_name2),
                                   color = c(
                                     "#8FDBF3","#A2C4F1","#CEBAF0","#C6C3E1",
                                     "#EEC2E5","#FFCFD1","#fbd37c","#f5c75b",
                                     "#FAA09C","#F9B29C","#F9C89B","#FBDF9D","#E9E4AF"),
                                   label = paste0("D", 1:length(dis_name2)))

create_circos_plot <- function(cohort_value) {
  # 处理数据
  df.rd3 <- df.rd %>%
    mutate(contrast = str_remove(contrast, " - No disease")) %>%
    mutate(contrast = case_when(
      contrast == "Psychological disorder" ~ "Psychological disorder",
      contrast == "Diabetes" ~ "Diabetes",
      contrast == "Physical–psychological multimorbidity" ~ "Physical–psychological multimorbidity",
      contrast == "Psychological–cognitive multimorbidity" ~ "Psychological–cognitive multimorbidity",
      contrast == "Hypertension" ~ "Hypertension",
      contrast == "Arthritis" ~ "Arthritis",
      contrast == "Stroke" ~ "Stroke",
      contrast == "Lung disease" ~ "Lung disease",
      contrast == "Heart disease" ~ "Heart disease",
      contrast == "Cancer" ~ "Cancer",
      contrast == "Physical–cognitive multimorbidity" ~ "Physical–cognitive multimorbidity",
      contrast == "Cognitive disorder" ~ "Cognitive disorder",
      contrast == "Physical–psychological–cognitive multimorbidity" ~ "Physical–psychological–cognitive multimorbidity",
      TRUE ~ contrast)) %>%
    filter(cohort == cohort_value) %>%
    filter(group != "Total" & RD > 0) %>%
    group_by(group) %>%
    arrange(desc(RD), .by_group = TRUE) %>%
    mutate(x_pos = (seq_along(RD) - 0.5) / n()) %>%
    mutate(y_pos = RD / max(RD)) %>%
    ungroup() %>%
    left_join(pastel_nature_colors, by = c("contrast" = "dis_name"))
  
  # 初始化
  circos.par(start.degree = 90,
             gap.after = c(rep(2, n_distinct(df.rd3$group) - 1), 20),
             track.margin = c(0, 0.03),
             cell.padding = c(0, 0, 0, 0))
  circos.initialize(factors = df.rd3$group, xlim = c(0, 1))
  
  # 提取每个组的标签数据
  label_data <- df.rd3 %>%
    group_by(group) %>%
    mutate(start = x_pos, end = x_pos) %>%
    ungroup() %>%
    select(group, start, end, label, color)
  
  # 添加基因名称标签
  circos.genomicLabels(label_data,
                       labels.column = 4,
                       facing = "reverse.clockwise",
                       side = "outside",
                       cex = 0.7,
                       col = "black",
                       connection_height = convert_height(4, "mm"),
                       line_lwd = 0.5)
  
  # 添加扇形分类标签
  circos.track(factors = df.rd3$group,
               ylim = c(0, 1),
               track.height = 0.08,
               bg.col = c("#f8bf99", "#fda453", "#f47720", "#87bdc0", "#5fadac"),
               panel.fun = function(x, y) {
                 circos.text(
                   CELL_META$xcenter,
                   CELL_META$ylim[2] - 0.7,
                   CELL_META$sector.index,
                   facing = "bending.inside",
                   cex = 0.8,
                   adj = c(0.5, 0)
                 )
               })
  
  # 绘制棒棒图
  circos.trackPlotRegion(
    factors = df.rd3$group,
    ylim = c(0, 1),
    track.height = 0.37,
    panel.fun = function(x, y) {
      sector_data <- df.rd3[df.rd3$group == CELL_META$sector.index, ]
      
      for (i in 1:nrow(sector_data)) {
        x_pos <- sector_data$x_pos[i]
        y_pos <- sector_data$y_pos[i]
        
        # 添加线条
        circos.lines(
          c(x_pos, x_pos),
          c(0, y_pos),
          col = sector_data$color[i],
          lwd = 3
        )
        
        # 添加点
        circos.points(
          x_pos,
          y_pos,
          pch = 16,
          col = sector_data$color[i],
          cex = 1.5
        )
      }
      
      # 在第一个扇区添加y轴刻度线
      if (CELL_META$sector.index == "<60") {
        circos.yaxis(
          side = "left",
          at = seq(0, 1, by = 0.2),
          labels.cex = 0.8,
          tick.length = convert_length(2, "mm")
        )
        
        # 添加y轴标题
        circos.text(
          CELL_META$cell.xlim[1] - convert_x(10, "mm"),
          0.5,
          "RD",
          facing = "reverse.clockwise",
          niceFacing = TRUE,
          adj = c(0.5, 0.5),
          cex = 0.8
        )
      }
    }
  )
}


pdf("./results/20230321/rd_charls.pdf", width = 6, height = 6)  # 创建 PDF 输出
create_circos_plot(cohort_value = "CHARLS")
dev.off()  # 关闭设备


pdf("./results/20230321/rd_elsa.pdf", width = 6, height = 6)  # 创建 PDF 输出
create_circos_plot(cohort_value = "ELSA")
dev.off()  # 关闭设备

pdf("./results/20230321/rd_hrs.pdf", width = 6, height = 6)  # 创建 PDF 输出
create_circos_plot(cohort_value = "HRS")
dev.off()  # 关闭设备

pdf("./results/20230321/rd_mhas.pdf", width = 6, height = 6)  # 创建 PDF 输出
create_circos_plot(cohort_value = "MHAS")
dev.off()  # 关闭设备

pdf("./results/20230321/rd_share.pdf", width = 6, height = 6)  # 创建 PDF 输出
create_circos_plot(cohort_value = "SHARE")
dev.off()  # 关闭设备




#-----绘制环状热图-------

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
# Load R packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(eulerr)))
suppressMessages(suppressWarnings(library(circlize)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(ComplexHeatmap)))
library(tidyverse)

pastel_nature_colors <- data.frame(dis_name=sort(dis_name2),
                                   label = paste0("D", 1:length(dis_name2)))

df.long <- loc_rank1 %>%
  pivot_longer(
    cols = -c(contrast, cause_rank),  # 排除 contrast 和 cause_rank 列
    names_to = c("group", "age_sex"),
    names_sep = " - ",  # 使用 " - " 作为分隔符
    values_to = "value"
  ) |> 
  pivot_wider(
    names_from = age_sex,  # 使用 age_sex 列的值来创建新的列
    values_from = value     # 使用 value 列的值填充新列
  ) |> 
  left_join(
    pastel_nature_colors,
    by = c("contrast" = "dis_name")
  )

# 转换数据为矩阵
dat_matrix <- df.long %>%
  select(`<60`, `60-74`, `75+`, `Female`, `Male`,`Total`) %>%
  as.matrix()
dat_matrix_rev <- dat_matrix[, ncol(dat_matrix):1]

# 设置行名
rownames(dat_matrix) <- df.long$label

# 定义颜色
group_colors <- c("CHARLS" = alpha("#DDCC77", 0.8), "ELSA" = alpha("#FF8C00", 0.8), 
                  "HRS" = alpha("#699ECA", 0.8), "MHAS" = alpha("#F898CB", 0.8), "SHARE" = alpha("#4DAF4A", 0.8))

green_pink <- colorRamp2(
  breaks = c(min(dat_matrix,na.rm = T), mean(dat_matrix,na.rm = T), max(dat_matrix,na.rm = T)),
  colors = c("#fbb01a", "#F7F7F7","#0070b8" )
)


# 设置分组的展示顺序
df.long$group <- factor(df.long$group, levels = c("CHARLS","ELSA","HRS","MHAS","SHARE"))


# Draw Plot
pdf("./results/20230321/rd_order.pdf", width = 10, height = 10)
# pdf(glue('{wkdir}/{no}.pdf'), width = 12, height = 12)

# 初始化 circlize
circos.clear()
circos.par(
  start.degree = 60,
  canvas.xlim = c(-0.9, 0.9),
  canvas.ylim = c(-0.9, 0.9),
  gap.after = c(rep(5, length(levels(df.long$group)) - 1), 30),
  track.margin = c(0, 0.03),
  cell.padding = c(0, 0, 0, 0)
)

# 绘制热图
circos.heatmap(
  dat_matrix,
  split = df.long$group,
  cluster = FALSE,
  bg.border = "black",
  bg.lwd = 1,
  cell.border = "white",
  cell.lwd = 0.5,
  rownames.side = "outside",
  rownames.cex = 1.2,
  col = green_pink,
  track.height = 0.25  # 增加track.height
)

# 添加数值标签
for(sector.index in get.all.sector.index()) {
  idx = which(df.long$group == sector.index)
  values = dat_matrix_rev[idx, , drop = FALSE]
  xlim = get.cell.meta.data("xlim", sector.index)
  ylim = get.cell.meta.data("ylim", sector.index)
  n_row = nrow(values)  # 这是每个扇区的行数 (14个)
  n_col = ncol(values)  # 这是列数 (5个)
  
  # 计算每个单元格的高度和宽度
  x_points = seq(xlim[1], xlim[2], length.out = n_row + 1)
  y_points = seq(ylim[1], ylim[2], length.out = n_col + 1)
  
  for(i in 1:n_row) {  # 遍历14个行
    x = (x_points[i] + x_points[i + 1]) / 2    # 计算x轴中心点
    for(j in 1:n_col) {  # 遍历5个列
      y = (y_points[j] + y_points[j + 1]) / 2  # 计算y轴中心点
      circos.text(
        x, y, sprintf("%.0f", values[i, j]), 
        sector.index = sector.index, 
        facing = "bending.inside", 
        cex = 0.8
      )
    }
  }
}

# 添加列名：在最后一个扇区的右侧
circos.track(
  track.index = get.current.track.index(),
  bg.border = NA,
  panel.fun = function(x, y){
    # 判断当前的扇区是否是环形图中的最后一个扇区
    if(CELL_META$sector.numeric.index == length(levels(df.long$group))) {
      # 计算每个标签的y坐标，确保均匀分布
      cn <- colnames(dat_matrix_rev)
      n <- length(cn)
      cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
      y_coords <- seq(CELL_META$cell.ylim[1] + cell_height / 2, CELL_META$cell.ylim[2] - cell_height / 2, length.out = n)
      #添加线段
      for (i in 1:n) {
        circos.lines(
          c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")), # x坐标，1mm偏移量
          c(y_coords[i], y_coords[i]),
          col = "black",
          lwd = 2
        )
      } 
      
      #添加文本
      circos.text(
        rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"), # x坐标，1.5mm偏移量
        y_coords, # y坐标
        cn,
        cex = 1.2,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }
)

# 添加扇形分组标签
circos.track(
  ylim = c(0, 1), 
  track.height = 0.065, 
  bg.col = adjustcolor(group_colors[levels(df.long$group)], alpha.f = 0.3),
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[2] - 0.75,
      CELL_META$sector.index, 
      facing = "bending.inside", 
      cex = 1.5, 
      adj = c(0.5, 0)
    )
  }
)

# Save to file
dev.off()

