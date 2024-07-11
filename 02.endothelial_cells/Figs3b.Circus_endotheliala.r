"f:/0-work/1.项目/4.全组织图谱/9.作图/13.endothelial/GO-BP"
#读取并查看数据集
df<-read.csv("23620.csv")
View(df)

# 计算百分比
df$fraction = df$n / sum(df$n)
# 计算累计百分比(每个矩形的顶部)
df$ymax = cumsum(df$fraction)
# 计算每个矩形的底部
df$ymin = c(0, head(df$ymax, n=-1))
#查看数据集
df

#加载r包
library(ggplot2)
#可视化绘图
ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=tissue)) +
  geom_rect() +
  coord_polar(theta="y") + 
  xlim(c(2, 4))

# 计算标签位置
df$labelPosition <- (df$ymax + df$ymin) / 2
# 计算标签内容
df$label <- paste0(df$Treat, "\n value: ", df$bio)
#查看数据集
df

# 可视化绘图
p1 <- ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = tissue)) +
    geom_rect() +
    geom_label(x = 3.5, aes(y = labelPosition, label = tissue), size = 6) +
    scale_fill_manual(values = tissue_colour) +  # 使用 scale_fill_manual() 函数
    coord_polar(theta = "y") +
    xlim(c(0.5, 4)) +
    theme_void() +
    theme(legend.position = "none")
# 可视化绘图
ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=tissue)) +
  geom_rect() +
  geom_text( x=2, aes(y=labelPosition, label=tissue, color=tissue), size=6) +
  scale_fill_manual(values = c('blood'='#ee5443',
'muscle'='#11abbc',
'liver'='#dcb150',
'heart'='#72dfb5',
'tailfin'='#56a22f',
'intestine'='#6e6f62',
'dorsal_fin'='#8d204a',
'kidney'='#234eaf',
'spinal_cord'='#91909d',
'gill'='#0b6a59',
'skin'='#a95dd0',
'notochord'='#614224',
'supraneural_body'='#f495a6',
'brain'='#9ebddb',
                   'other'='#dddfe6')) +
  scale_color_manual(values = c('black','black','black','black','black','black','black','black','black','black','black','black')) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")


tissue_colour <- c('blood'='#ee5443',
'muscle'='#11abbc',
'liver'='#dcb150',
'heart'='#72dfb5',
'tailfin'='#56a22f',
'intestine'='#6e6f62',
'dorsal_fin'='#8d204a',
'kidney'='#234eaf',
'spinal_cord'='#91909d',
'gill'='#0b6a59',
'skin'='#a95dd0',
'notochord'='#614224',
'supraneural_body'='#f495a6',
'brain'='#9ebddb',
                   'other'='#dddfe6'
)

pdf('Circus_endothelial.pdf',width=10,height = 10)
p1
dev.off()