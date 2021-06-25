

#fixation
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

eye_data_back <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/mri_eyedata_back_4s.txt", stringsAsFactors=FALSE)
eye_data_left <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/mri_eyedata_left_4s.txt", stringsAsFactors=FALSE)
eye_data_right <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/mri_eyedata_right_4s.txt", stringsAsFactors=FALSE)
eye_data_back <- eye_data_back/4
eye_data_left <- eye_data_left/4
eye_data_right <- eye_data_right/4

eye_data_back <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/mri_eyedata_back_2s.txt", stringsAsFactors=FALSE)
eye_data_left <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/mri_eyedata_left_2s.txt", stringsAsFactors=FALSE)
eye_data_right <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/mri_eyedata_right_2s.txt", stringsAsFactors=FALSE)
eye_data_back <- eye_data_back/2
eye_data_left <- eye_data_left/2
eye_data_right <- eye_data_right/2

mean(eye_data_back[,2])
mean(eye_data_left[,2])
mean(eye_data_right[,2])

eye_data_back = eye_data_back[-which(eye_data_back[,1]==0),]
eye_data_left = eye_data_left[-which(eye_data_left[,1]==0),]
eye_data_right = eye_data_right[-which(eye_data_right[,1]==0),]

no_sub = length(eye_data_back[,1])


## anova
datapool<- data.frame(
  mean = c(eye_data_back[,2], eye_data_left[,2], eye_data_right[,2]),
  x_label = c(rep('back',length(eye_data_back[,1])),rep('left',length(eye_data_back[,1])), rep('right',length(eye_data_back[,1])))
)
fit <- aov(datapool$mean ~ datapool$x_label, data=datapool)
print(summary(fit))


fix_or_sacca = 2
table <- data.frame(
  mean = c(rep(mean(eye_data_back[,fix_or_sacca]), no_sub), rep(mean(eye_data_left[,fix_or_sacca]), no_sub), rep(mean(eye_data_right[,fix_or_sacca]),no_sub)),
  se = c(rep(sd(eye_data_back[,fix_or_sacca])/sqrt(no_sub), no_sub), rep(sd(eye_data_left[,fix_or_sacca])/sqrt(no_sub), no_sub), rep(sd(eye_data_right[,fix_or_sacca])/sqrt(no_sub), no_sub)),
  sd = c(rep(sd(eye_data_back[,fix_or_sacca]), no_sub), rep(sd(eye_data_left[,fix_or_sacca]), no_sub), rep(sd(eye_data_right[,fix_or_sacca]), no_sub)),
  x_label = c(rep('Back', no_sub),rep('Left', no_sub),rep('Right', no_sub)),
  individual = c(eye_data_back[,fix_or_sacca], eye_data_left[,fix_or_sacca], eye_data_right[,fix_or_sacca]),
  color = c(rep("dodgerblue4",no_sub), rep("gold4",no_sub), rep("goldenrod1",no_sub))
)
table$x_label <- factor(table$x_label, levels= c('Back','Left','Right'))

font_size = 33
scaleFUN <- function(x) sprintf("%.1f", x)
library(ggplot2)
figure <- ggplot(data=table, aes(x=x_label, y=mean)) +  
  geom_bar(position=position_dodge(), stat="identity", fill=table$color, width=.3) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), colour="black", width=.15) +
  #geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
  scale_fill_manual("legend", values = c("Back" = "dodgerblue4", "Left" = "gold4", "Right" = "goldenrod1")) +
  theme_classic() +    
  #labs(x="",y="\n")+  #Performance
  #labs(x="",y="Performance\n")+  #
  labs(x="",y="\n") + 
  scale_y_continuous(labels=scaleFUN, limits = c(0, 4), breaks = seq(0, 4, by = 2)) +
  #ggtitle(paste("","\n", sep = "")) +
  theme(axis.title.y = element_text(size=font_size))+
  theme(axis.title.x = element_text(size=font_size))+
  theme(legend.text=element_text(size=font_size))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
  theme(axis.text.y = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
  theme(legend.title=element_blank())+
  guides(fill = guide_legend(override.aes = list(shape = NA)))
ggsave(paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/mri_saccade_4s.pdf"), #meg_beh.pdf mri_beh.pdf
       plot = figure, # or give ggplot object name as in myPlot,
       width = 5.5, height = 4.5,
       units = "in", # other options c("in", "cm", "mm"),
       dpi = 300)


