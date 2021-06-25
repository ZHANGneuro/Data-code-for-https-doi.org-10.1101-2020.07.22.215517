



############# MRI
rm(list=ls(all=TRUE)) 
library(ggplot2) #ggplot() for plotting
library(reshape2) #melt(), dcast() for data reformatting
library(plyr) #ddply() for data reformatting
library(Cairo) #better aliasing of output images
library(gridExtra) # for grid.arrange
library(grid) 
library(qdapRegex)
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)


bar_matrix <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/position_mri/mri_fix_sac_visual_angle_back_2s.txt", stringsAsFactors=FALSE)
no_sub = length(unique(bar_matrix[,1]))

## anova  fix 3 4  saccade 5, 6
datapool<- data.frame(
  data = bar_matrix[,5],
  ego = rep(c("1","2","3"), no_sub)
)
fit <- aov(datapool$data ~ datapool$ego, data=datapool)
print(summary(fit))

datapool<- data.frame(
  data = bar_matrix[,6],
  ego = rep(c("1","2","3"), no_sub)
)
fit <- aov(datapool$data ~ datapool$ego, data=datapool)
print(summary(fit))


pool_left <- bar_matrix[which(bar_matrix[,2]==1),]
pool_right <- bar_matrix[which(bar_matrix[,2]==2),]
pool_back <- bar_matrix[which(bar_matrix[,2]==3),]
font_size = 35
scaleFUN <- function(x) sprintf("%.2f", x)
bar_width=0.7


table <- data.frame(
  mean = c(rep(mean(pool_left[,3]), no_sub), rep(mean(pool_left[,4]), no_sub),
          rep(mean(pool_right[,3]), no_sub), rep(mean(pool_right[,4]), no_sub),
          rep(mean(pool_back[,3]), no_sub), rep(mean(pool_back[,4]), no_sub)),
  se = c(rep(sd(pool_left[,3]/sqrt(no_sub)), no_sub), rep(sd(pool_left[,4]/sqrt(no_sub)), no_sub),
        rep(sd(pool_right[,3]/sqrt(no_sub)), no_sub), rep(sd(pool_right[,4]/sqrt(no_sub)), no_sub),
        rep(sd(pool_back[,3]/sqrt(no_sub)), no_sub), rep(sd(pool_back[,4]/sqrt(no_sub)), no_sub)),
  individual = c(pool_left[,3], pool_left[,4],
                 pool_right[,3], pool_right[,4],
                 pool_back[,3], pool_back[,4]),
  ego = c(rep('Left', no_sub*2),rep('Right', no_sub*2),rep('Back', no_sub*2)),
  coor = c(rep('Horizontal',no_sub), rep('Vertical',no_sub), rep('Horizontal',no_sub), rep('Vertical',no_sub), rep('Horizontal',no_sub), rep('Vertical',no_sub))
)
table$ego <- factor(table$ego, levels= c('Back','Left','Right'))

library(ggplot2)
figure <- ggplot(data=table, aes(x=coor, y=mean, fill=ego)) +  
  geom_bar(position=position_dodge(), stat="identity", width=bar_width) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=0.5, position=position_dodge(bar_width)) +
  scale_fill_manual("legend", values = c("Back" = "dodgerblue4", "Left" = "gold4", "Right" = "goldenrod1")) +
  #geom_point(aes(x=table$ego,y=table$individual),colour="red",size=6,shape=1) +
  theme_classic() +  
  #labs(x="",y="\n")+  #Performance
  #labs(x="",y="Performance\n")+  #
  labs(x="",y="Va.of fixation \n") + #Saccade
  scale_y_continuous(labels=scaleFUN, limits = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, by = 0.5)) +
  #ggtitle(paste("","\n", sep = "")) +
  theme(axis.title.y = element_text(size=font_size))+
  theme(axis.title.x = element_text(size=font_size))+
  #theme(legend.text=element_text(size=font_size), legend.direction = "vertical", legend.position=c(.7,.75))+
  theme(legend.text=element_text(size=font_size))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size=font_size, angle = 20, vjust = 0.5)) + 
  theme(axis.text.y = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
  theme(legend.title=element_blank())+
  guides(fill=FALSE)
ggsave(paste("/Users/bo/Desktop/mri_va_fixation_2s.pdf"), #meg_beh.pdf mri_beh.pdf mri_va_fixation_2s.pdf
       plot = figure, # or give ggplot object name as in myPlot,
       width = 6.5, height = 5,
       units = "in", # other options c("in", "cm", "mm"),
       dpi = 300)



table <- data.frame(
  mean = c(rep(mean(pool_left[,5]), no_sub), rep(mean(pool_left[,6]), no_sub),
           rep(mean(pool_right[,5]), no_sub), rep(mean(pool_right[,6]), no_sub),
           rep(mean(pool_back[,5]), no_sub), rep(mean(pool_back[,6]), no_sub)),
  se = c(rep(sd(pool_left[,5]/sqrt(no_sub)), no_sub), rep(sd(pool_left[,6]/sqrt(no_sub)), no_sub),
         rep(sd(pool_right[,5]/sqrt(no_sub)), no_sub), rep(sd(pool_right[,6]/sqrt(no_sub)), no_sub),
         rep(sd(pool_back[,5]/sqrt(no_sub)), no_sub), rep(sd(pool_back[,6]/sqrt(no_sub)), no_sub)),
  individual = c(pool_left[,5], pool_left[,6],
                 pool_right[,5], pool_right[,6],
                 pool_back[,5], pool_back[,6]),
  ego = c(rep('Left', no_sub*2),rep('Right', no_sub*2),rep('Back', no_sub*2)),
  coor = c(rep('Horizontal',no_sub), rep('Vertical',no_sub), rep('Horizontal',no_sub), rep('Vertical',no_sub), rep('Horizontal',no_sub), rep('Vertical',no_sub))
)
table$ego <- factor(table$ego, levels= c('Back','Left','Right'))

scaleFUN <- function(x) sprintf("%.1f", x)
library(ggplot2)
figure <- ggplot(data=table, aes(x=coor, y=mean, fill=ego)) +  
  geom_bar(position=position_dodge(), stat="identity", width=bar_width) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=0.5, position=position_dodge(bar_width)) +
  #geom_point(aes(x=table$ego,y=table$individual),colour="red",size=6,shape=1) +
  scale_fill_manual("legend", values = c("Back" = "dodgerblue4", "Left" = "gold4", "Right" = "goldenrod1")) +
  theme_classic() +  
  #labs(x="",y="\n")+  #Performance
  #labs(x="",y="Performance\n")+  #
  labs(x="",y="\n") + #Va.of saccade
  scale_y_continuous(labels=scaleFUN, limits = c(-1.5, 2), breaks = seq(-1.5, 2, by = 1)) +
  theme(axis.title.y = element_text(size=font_size))+
  theme(axis.title.x = element_text(size=font_size))+
  #theme(legend.text=element_text(size=font_size), legend.direction = "vertical", legend.position=c(.7,.75))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size=font_size, angle = 20, vjust = 0.5)) + 
  theme(axis.text.y = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
  theme(legend.title=element_blank())+
  guides(fill=FALSE)
ggsave(paste("/Users/bo/Desktop/mri_va_saccade_2s.pdf"), #meg_beh.pdf mri_beh.pdf
       plot = figure, # or give ggplot object name as in myPlot,
       width = 6.5, height = 5,
       units = "in", # other options c("in", "cm", "mm"),
       dpi = 300)













