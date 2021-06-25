


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


# four col, 1) tp fix  2) tp sac 3) tnp fix 4)tnp sac

eye_data_back <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/meg_eyedata_back_both1&2s.txt", stringsAsFactors=FALSE)
eye_data_left <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/meg_eyedata_left_both1&2s.txt", stringsAsFactors=FALSE)
eye_data_right <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/meg_eyedata_right_both1&2s.txt", stringsAsFactors=FALSE)
eye_data_control <- read.table("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/meg_eyedata_control_both1&2s.txt", stringsAsFactors=FALSE)

no_sub = length(eye_data_back[,1])

fix_or_sacca_t = 2
fix_or_sacca_tn = 4
data_1 <- eye_data_back[,fix_or_sacca_t] + eye_data_back[,fix_or_sacca_tn]
data_2 <- eye_data_left[,fix_or_sacca_t] + eye_data_left[,fix_or_sacca_tn]
data_3 <- eye_data_right[,fix_or_sacca_t] + eye_data_right[,fix_or_sacca_tn]

data_1 <- eye_data_back[,fix_or_sacca_t] 
data_2 <- eye_data_left[,fix_or_sacca_t] 
data_3 <- eye_data_right[,fix_or_sacca_t] 

data_1 <- data_1/2
data_2 <- data_2/2
data_3 <- data_3/2

table <- data.frame(
  mean = c(rep(mean(data_1), no_sub), rep(mean(data_2), no_sub), rep(mean(data_3),no_sub)),
  se = c(rep(sd(data_1)/sqrt(no_sub), no_sub), rep(sd(data_2)/sqrt(no_sub), no_sub), rep(sd(data_3)/sqrt(no_sub), no_sub)),
  sd = c(rep(sd(data_1), no_sub), rep(sd(data_2), no_sub), rep(sd(data_3), no_sub)),
  x_label = c(rep('Back', no_sub),rep('Left', no_sub),rep('Right', no_sub)),
  individual = c(data_1, data_2, data_3),
  color = c(rep("dodgerblue4",no_sub), rep("gold4",no_sub), rep("goldenrod1",no_sub))
)
table$x_label <- factor(table$x_label, levels= c('Back','Left','Right'))

## anova
datapool<- data.frame(
  data = c(data_1, data_2, data_3),
  cond = c(rep('back',10),rep('left',10), rep('right',10))
)
fit <- aov(datapool$data ~ datapool$cond, data=datapool)
print(summary(fit))


font_size = 33
scaleFUN <- function(x) sprintf("%.1f", x)
library(ggplot2)
figure <- ggplot(data=table, aes(x=x_label, y=mean)) +  
  geom_bar(position=position_dodge(), stat="identity", fill=table$color, width=.3) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), colour="black", width=.15) +
  #geom_point(aes(x=x_label,y=individual),colour="red",size=6,shape=1) +
  scale_fill_manual("legend", values = c("Back" = "dodgerblue4", "Left" = "gold4", "Right" = "goldenrod1")) +
  theme_classic() +    
  #labs(x="",y="Frequency of saccade [Hz]\n") + #fixation saccade
  labs(x="",y="\n") +
  scale_y_continuous(labels=scaleFUN, limits = c(0, 4), breaks = seq(0, 4, by = 2)) +
  theme(axis.title.y = element_text(size=font_size))+
  theme(axis.title.x = element_text(size=font_size))+
  theme(legend.text=element_text(size=font_size))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
  theme(axis.text.y = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
  theme(legend.title=element_blank())+
  guides(fill = guide_legend(override.aes = list(shape = NA)))
ggsave(paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/meg_saccade_1s.pdf"), #meg_beh.pdf mri_beh.pdf fixation
       plot = figure, # or give ggplot object name as in myPlot,
       width = 5.5, height = 4.5,
       units = "in", # other options c("in", "cm", "mm"),
       dpi = 300)


