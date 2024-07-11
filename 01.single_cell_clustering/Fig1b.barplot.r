#log10(number of cells/nuclei)

cols = c('blood'='#234eaf',
'brain'='#11abbc',
'dorsal_fin'='#72dfb5',
'gill'='#56a22f',
'heart'='#0b6a59',
'intestine'='#614224',
'kidney'='#6e6f62',
'liver'='#91909d',
'muscle'='#9ebddb',
'notochord'='#f495a6',
'skin'='#a95dd0',
'spinal_cord'='#8d204a',
'supraneural_body'='#ee5443',
'tailfin'='#dcb150')

library(ggplot2)
data$class <- factor(data$class, levels = rev(unique(data$class)))
cols<-c("#af2934","#ffe327","#2f4e87","#b0b9b8","#f0eedf",
        "#aed4e9","#f4a69a","#4593c3","#f18e0c",
        "#262a35","#c5942e","#a2a7ab","#3ba889")
ggplot(data, aes(x = count, y = class, fill = class)) +
    geom_col(width = 0.9) +
    labs(x = "Cell Number", y = "Tissue")+
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 22),  # 调整xy轴标题大小为22
          axis.text = element_text(size = 16),   # 调整xy轴刻度标签大小为16
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white")) +
    scale_fill_manual(values = cols)