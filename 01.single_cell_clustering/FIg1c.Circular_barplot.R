

library(tidyverse)
library(viridis)
 library(Seurat)
library(tidyverse)
library(ggrepel)
data <- read.csv('melted_dat.csv')
head(data)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(data$observation))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(group, individual)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label
data$value <- as.numeric(data$value)
label_data <- data %>% group_by(id, individual) %>% summarize(tot=sum(value))
##label_data <- data %>% group_by(id, individual) %>% summarize(tot = sum(value), .groups = 'drop')
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
 
# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]



celltype_to_color <- c('blood'='#ee5443',
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
)

allcolor = c()
p <- ggplot(data) +      
    
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.5) +
    scale_fill_viridis(discrete=TRUE) +
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 150, xend = start, yend = 150), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    
    # Add text showing the value of each 100/75/50/25 lines
    ggplot2::annotate("text", x = rep(max(data$id),5), y = c(0, 50, 100, 150, 200), label = c("0", "50", "100", "150", "200") , color="grey", size=6 , angle=0, fontface="bold", hjust=1) +
    
    ylim(-150,max(label_data$tot, na.rm=T)) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    
    # Add labels on top of each bar
    geom_text(data=label_data, aes(x=id, y=tot+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
    
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
    geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1, 0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

p1 <- ggplot(data) +      
    geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=1) +
    scale_fill_manual(values = celltype_to_color) +
    geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0", "25", "50", "75", "100") , size=6 , angle=0, fontface="bold", hjust=1) +
    theme(legend.position = "none", axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
    ylim(-150,max(label_data$tot, na.rm=T)) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    geom_text(data=label_data, aes(x=id, y=tot+10, label=individual, hjust=hjust), colour = "black", fontface="bold", alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
    #geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1, 0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
    geom_text(data=base_data, aes(x = (start+end)/2, y = -18, label=group), hjust=0.5, colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
p1


p2 <- ggplot(data) +      
    geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=1) +
    scale_fill_manual(values = celltype_to_color) +
    geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0", "25", "50", "75", "100") , size=6 , angle=0, fontface="bold", hjust=1) +
    theme(legend.position = "none", axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
    ylim(-150,max(label_data$tot, na.rm=T)) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    geom_text(data=label_data, aes(x=id, y=tot+10, label=id, hjust=hjust), colour = "black", fontface="bold", alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
    geom_point(data=data, aes(x=as.factor(id), y=-20, fill=as.factor(id)), shape=21, colour="black", size=5, position=position_dodge(width=0)) +
    scale_fill_manual(values = allcolor)
p2


allcolor = c('1'='#34362d','2'='#ff90c9','3'='#eec3ff','4'='#885578','5'='#63ffac','6'='#00a6aa',
'7'='#b903aa','8'='#372101','9'='#FFFFFF','10'='#FFFFFF','11'='#d16100','12'='#6f0062',
'13'='#c2ffed','14'='#a30059','15'='#ffb500','16'='#300018','17'='#00489c','18'='#7a87a1',
'19'='#a1c299','20'='#00846f','21'='#008941','22'='#809693','23'='#a3c8c9','24'='#1ce6ff',
'25'='#000035','26'='#c2ff99','27'='#b4a8bd','28'='#b77b68','29'='#ff913f','30'='#788d66',
'31'='#00c2a0','32'='#fad09f','33'='#452c2c','34'='#ff4a46','35'='#FFFFFF','36'='#FFFFFF',
'37'='#013349','38'='#ff2f80','39'='#ba0900','40'='#61615a','41'='#7a4900','42'='#0aa6d8',
'43'='#6a3a4c','44'='#004d43','45'='#5a0007','46'='#ffff00','47'='#ddefff','48'='#636375',
'49'='#ff8a9a','50'='#c0b9b2','51'='#6b7900','52'='#938a81','53'='#0086ed','54'='#4a3b53',
'55'='#001e09','56'='#FFFFFF','57'='#FFFFFF','58'='#997d87','59'='#bec459','60'='#1b4400',
'61'='#d157a0','62'='#ffdbe5','63'='#886f4c','64'='#8fb0ff','65'='#7b4f4b','66'='#456648',
'67'='#4fc601','68'='#cc0744','69'='#ffaa92','70'='#0000a6','71'='#3b5dff','72'='#006fa6',
'73'='#b79762',
'74'='#456d75','75'='#ff34ff','76'='#0cbd66','77'='#FFFFFF','78'='#FFFFFF')


#12.11
p3 <- ggplot(data) +      
    geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=1) +
    scale_fill_manual(values = celltype_to_color) +
    geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0", "25", "50", "75", "100") , size=6 , angle=0, fontface="bold", hjust=1) +
    theme(legend.position = "none", axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
    ylim(-150,max(label_data$tot, na.rm=T)) +
    theme_minimal() +
    theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    geom_text(data=label_data, aes(x=id, y=tot+10, label=label, hjust=hjust), colour = "black", fontface="bold", alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
    geom_point(data=data, aes(x=as.factor(id), y=-20, fill=as.factor(id)), shape=21, colour="black", size=5, position=position_dodge(width=0)) +
    scale_fill_manual(values = allcolor)
p3