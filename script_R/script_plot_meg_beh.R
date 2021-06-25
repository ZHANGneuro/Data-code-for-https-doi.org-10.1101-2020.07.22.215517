library(ggplot2) #ggplot() for plotting

#rm(list=ls(all=TRUE)) 
# smt
rate_table_meg = c()
for (sub_num in c(2:13)) {
  if(sub_num==2){
    path_s3 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s3_rawdata.txt", sep = "")
    path_s4 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s4_rawdata.txt", sep = "")
    dataset_s3 <-  read.table(path_s3,stringsAsFactors=FALSE)
    dataset_s4 <-  read.table(path_s4,stringsAsFactors=FALSE)
    t_dataset <- rbind(dataset_s3,dataset_s4)
    t_dataset = add_correct_marker(t_dataset)
  }else if (sub_num==3){
    path_s1 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s1_rawdata.txt", sep = "")
    path_s2 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s2_rawdata.txt", sep = "")
    dataset_s1 <-  read.table(path_s1,stringsAsFactors=FALSE)
    dataset_s2 <-  read.table(path_s2,stringsAsFactors=FALSE)
    t_dataset <- rbind(dataset_s1,dataset_s2)
    t_dataset = add_correct_marker(t_dataset)
  }else{
    path_s1 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s1_rawdata.txt", sep = "")
    path_s2 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s2_rawdata.txt", sep = "")
    path_s3 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s3_rawdata.txt", sep = "")
    path_s4 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s4_rawdata.txt", sep = "")
    dataset_s1 <-  read.table(path_s1,stringsAsFactors=FALSE)
    dataset_s2 <-  read.table(path_s2,stringsAsFactors=FALSE)
    dataset_s3 <-  read.table(path_s3,stringsAsFactors=FALSE)
    dataset_s4 <-  read.table(path_s4,stringsAsFactors=FALSE)
    t_dataset <- rbind(dataset_s1,dataset_s2,dataset_s3,dataset_s4)
    t_dataset = add_correct_marker(t_dataset)
  }
  returned = correct_rate_ego_each_subject_meg(t_dataset,14)
  rate_table_meg <- rbind(rate_table_meg, returned[,1])
}

# MRI
path = "/Users/bo/Documents/data_yuji_lab/data_fmri/raw_behav_date_for_tempUse/"
files <- as.character(list.files(path=path, pattern = "formal_rawdata_t"))
rate_table_fmri = c()
for(i in 8:26){
  figure_name <- paste("Sub_", i, sep = "")
  dataset <-  read.table(paste(path, "sub_", as.character(i), "_formal_rawdata_t.txt", sep = ""),stringsAsFactors=FALSE)
  dataset_main <- dataset[which(is.na(dataset[,6])),]
  dataset_main = add_correct_marker_mri(dataset_main,10, 13, FALSE)
  returned = correct_rate_ego_each_subject(dataset_main,14)
  rate_table_fmri = rbind(rate_table_fmri, returned[,1])
}

meg_list <- data.frame(
  value = c(rate_table_meg[,1], rate_table_meg[,2], rate_table_meg[,3]),
  cate = c(rep('1', 12), rep('2',12), rep('3',12))
)
summary(aov(meg_list$value ~ meg_list$cate))

fmri_list <- data.frame(
  value = c(rate_table_fmri[,1], rate_table_fmri[,2], rate_table_fmri[,3]),
  cate = c(rep('1', 19), rep('2',19), rep('3',19))
)
summary(aov(meg_list$value ~ meg_list$cate))




#write.table(rate_table_fmri[,3], paste("/Users/boo/Desktop/MEG_data_script/correlation/performance_mri_back.txt", sep = ""), sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)



rate_table = rate_table_meg  #rate_table_fmri rate_table_meg

num_sub = length(rate_table[,1])
table <- data.frame(
  individual = c(rate_table[,1], rate_table[,2], rate_table[,3]),
  mean = c(rep(mean(rate_table[,1]),num_sub), rep(mean(rate_table[,2]),num_sub), rep(mean(rate_table[,3]),num_sub)),
  #se = c(rep(sd(rate_table[,1]),num_sub), rep(sd(rate_table[,2]),num_sub), rep(sd(rate_table[,3]),num_sub)),
  sd = c(rep(sd(rate_table[,1]),num_sub), rep(sd(rate_table[,2]),num_sub), rep(sd(rate_table[,3]),num_sub)),
  x_label = c(rep('Left',num_sub), rep('Right',num_sub), rep('Back',num_sub)),
  #color = c(rep(rgb(89/255,88/255,89/255),num_sub), rep(rgb(89/255,88/255,89/255),num_sub), rep(rgb(89/255,88/255,89/255),num_sub))
  color = c(rep("gold4",num_sub), rep("goldenrod1",num_sub),rep("dodgerblue4",num_sub))
  )
table$x_label = factor(table$x_label, levels=c('Back', 'Left', 'Right'))

font_size = 33
library(ggplot2)
figure <- ggplot(data=table, aes(x=x_label, y=mean)) +  
  geom_bar(position=position_dodge(), stat="identity", fill=table$color, width=.3) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), colour="black", width=.15) +
  scale_fill_manual("legend", values = c("Back" = "dodgerblue4", "Left" = "gold4", "Right" = "goldenrod1")) +
  #geom_hline(yintercept=0.25, linetype="dashed", color = "orange") +
  #geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
  theme_classic() +    
  #labs(x="",y="\n")+  #Performance
  labs(x="",y="\n")+  #Performance
  ggtitle(paste("","", sep = "")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, by = 0.5)) +
  theme(axis.title.y = element_text(size=font_size)) +
  theme(axis.title.x = element_text(size=font_size)) +
  theme(legend.text=element_text(size=font_size)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
  theme(axis.text.y = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
  theme(legend.title=element_blank())+
  guides(fill = guide_legend(override.aes = list(shape = NA)))
ggsave(paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyeData_script_result/meg_beh.pdf"), #meg_beh.pdf mri_beh.pdf
       plot = figure, # or give ggplot object name as in myPlot,
       width = 5.5, height = 4.5,
       units = "in", # other options c("in", "cm", "mm"),
       dpi = 300)



# fmri spatial task
table <- data.frame(
  individual = c(rate_table[,1],rate_table[,2],rate_table[,3],rate_table[,4]),
  mean = c(rep(mean(rate_table[,1]),10), rep(mean(rate_table[,2]),10),
           rep(mean(rate_table[,3]),10), rep(mean(rate_table[,4]),10)),
  se = c(rep(std(rate_table[,1])/sqrt(10),10), rep(std(rate_table[,2])/sqrt(10),10),
         rep(std(rate_table[,3])/sqrt(10),10), rep(std(rate_table[,4])/sqrt(10),10)),
  x_label = c(rep('dc',10), rep('sc',10),
              rep('cfp',10), rep('ctp',10)))

table$x_label = factor(table$x_label, levels=c('dc', 'sc', 'ctp', 'cfp'))

font_size = 30
library(ggplot2)
figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
  geom_bar(position=position_dodge(), stat="identity", width=.15) +
  geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.15) +
  geom_hline(yintercept=0.25, linetype="dashed", color = "orange")+
  geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
  theme_classic() +    
  labs(x="",y="Performance\n")+ 
  ggtitle(paste("","\n\n", sep = "")) +
  theme(axis.title.y = element_text(size=font_size))+
  theme(axis.title.x = element_text(size=font_size))+
  theme(legend.text=element_text(size=font_size))+
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + 
  theme(axis.text.x = element_text(size=font_size, angle = 90, vjust = 0.5)) + 
  theme(axis.text.y = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
  theme(legend.title=element_blank())+
  guides(fill = guide_legend(override.aes = list(shape = NA)))
figure















# false alarm rate day 1
rate_table_md= c()
for (sub_num in c(1:12)) {
  path_temp <- paste("/Users/boo/Desktop/MEG_data_script/day1_data/sub_", sub_num, "_day1_beh_data.txt", sep = "")
  dataset_hnd <-  read.table(path_temp,stringsAsFactors=FALSE)
  temp_rate <- training_false_alarm_rate_each_subject(dataset_hnd)
  rate_table_md <- rbind(rate_table_md, temp_rate$y)
}

# false alarm rate day 2
rate_table_fmri = c()
for (sub_num in c(4:13)) {
  path_s1 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s1_rawdata.txt", sep = "")
  path_s2 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s2_rawdata.txt", sep = "")
  path_s3 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s3_rawdata.txt", sep = "")
  path_s4 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s4_rawdata.txt", sep = "")
  
  dataset_s1 <-  read.table(path_s1,stringsAsFactors=FALSE)
  dataset_s2 <-  read.table(path_s2,stringsAsFactors=FALSE)
  dataset_s3 <-  read.table(path_s3,stringsAsFactors=FALSE)
  dataset_s4 <-  read.table(path_s4,stringsAsFactors=FALSE)
  dataset <- rbind(dataset_s1,dataset_s2,dataset_s3,dataset_s4)
  dataset = add_correct_marker(dataset)
  
  dataset_hnd <- dataset[which(dataset[,2]=="hnd"),]
  temp_rate <- meg_false_alarm_rate_each_subject(dataset_hnd)
  rate_table_fmri <- rbind(rate_table_fmri, temp_rate$y)
}

figure <- plot_combined_false_alarm_rate(rate_table_md, rate_table_fmri, 20)




# rt
pool_front = c()
pool_left = c()
pool_right = c()
pool_back = c()

for (sub_num in c(4:13)) {
  path_s1 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s1_rawdata.txt", sep = "")
  path_s2 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s2_rawdata.txt", sep = "")
  path_s3 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s3_rawdata.txt", sep = "")
  path_s4 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s4_rawdata.txt", sep = "")
  rt_path_s1 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s1_timing.txt", sep = "")
  rt_path_s2 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s2_timing.txt", sep = "")
  rt_path_s3 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s3_timing.txt", sep = "")
  rt_path_s4 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s4_timing.txt", sep = "")
  
  dataset_s1 <-  read.table(path_s1,stringsAsFactors=FALSE)
  dataset_s2 <-  read.table(path_s2,stringsAsFactors=FALSE)
  dataset_s3 <-  read.table(path_s3,stringsAsFactors=FALSE)
  dataset_s4 <-  read.table(path_s4,stringsAsFactors=FALSE)
  timing_s1 <-  read.table(rt_path_s1,stringsAsFactors=FALSE)
  timing_s2 <-  read.table(rt_path_s2,stringsAsFactors=FALSE)
  timing_s3 <-  read.table(rt_path_s3,stringsAsFactors=FALSE)
  timing_s4 <-  read.table(rt_path_s4,stringsAsFactors=FALSE)
  
  dataset <- rbind(dataset_s1,dataset_s2,dataset_s3,dataset_s4)
  dataset = add_correct_marker(dataset)
  timing <- rbind(timing_s1,timing_s2,timing_s3,timing_s4)
  if (sub_num==4){
    timing[,10] <- timing[,9] - timing[,8]
  } else {
    cue_onset <- c()
    press_onset <- c()
    for (ith_row in 1:length(timing[,1])) {
      split_str_cue <- strsplit(timing[ith_row, 8],"_")
      split_str_press <- strsplit(timing[ith_row, 9],"_")
      cue_onset <- c(cue_onset, as.numeric(split_str_cue[[1]][2]))
      press_onset <- c(press_onset, as.numeric(split_str_press[[1]][2]))
    }
    timing[,10] <- press_onset - cue_onset
  }
  
  temp_front <- timing[which(dataset[,11]=='front'), 10]
  temp_left <- timing[which(dataset[,11]=='left'), 10]
  temp_right <- timing[which(dataset[,11]=='right'), 10]
  temp_back <- timing[which(dataset[,11]=='back'), 10]
  
  pool_front = c(pool_front, mean(temp_front, na.rm=TRUE))
  pool_left = c(pool_left,  mean(temp_left, na.rm=TRUE))
  pool_right = c(pool_right,  mean(temp_right, na.rm=TRUE))
  pool_back = c(pool_back,  mean(temp_back, na.rm=TRUE))
}

mean(pool_front)
mean(pool_left)
mean(pool_right)
mean(pool_back)













# testing
sub_num = 13
path_s1 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s1_rawdata.txt", sep = "")
path_s2 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s2_rawdata.txt", sep = "")
path_s3 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s3_rawdata.txt", sep = "")
path_s4 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s4_rawdata.txt", sep = "")

dataset_s1 <-  read.table(path_s1,stringsAsFactors=FALSE)
dataset_s2 <-  read.table(path_s2,stringsAsFactors=FALSE)
dataset_s3 <-  read.table(path_s3,stringsAsFactors=FALSE)
dataset_s4 <-  read.table(path_s4,stringsAsFactors=FALSE)
dataset <- rbind(dataset_s1,dataset_s2,dataset_s3,dataset_s4)

length(which(dataset_s1[,11]=='left'))
length(which(dataset_s1[,11]=='right'))
length(which(dataset_s1[,11]=='back'))

length(which(dataset_s2[,11]=='left'))
length(which(dataset_s2[,11]=='right'))
length(which(dataset_s2[,11]=='back'))

length(which(dataset_s3[,11]=='left'))
length(which(dataset_s3[,11]=='right'))
length(which(dataset_s3[,11]=='back'))

length(which(dataset_s4[,11]=='left'))
length(which(dataset_s4[,11]=='right'))
length(which(dataset_s4[,11]=='back'))






# smt
rate_table = c()
for (sub_num in c(4:13)) {
  path_s1 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s1_rawdata.txt", sep = "")
  path_s2 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s2_rawdata.txt", sep = "")
  path_s3 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s3_rawdata.txt", sep = "")
  path_s4 <- paste("/Users/bo/Documents/data_yuji_lab/data_MEG/eyedata_meg/sub_", sub_num, "_s4_rawdata.txt", sep = "")
  
  dataset_s1 <-  read.table(path_s1,stringsAsFactors=FALSE)
  dataset_s2 <-  read.table(path_s2,stringsAsFactors=FALSE)
  dataset_s3 <-  read.table(path_s3,stringsAsFactors=FALSE)
  dataset_s4 <-  read.table(path_s4,stringsAsFactors=FALSE)
  
  write.table(dataset_s1, paste("/Users/boo/Desktop/temp/sub_", sub_num, "_s1_rawdata.txt", sep = ""), sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
  write.table(dataset_s2, paste("/Users/boo/Desktop/temp/sub_", sub_num, "_s2_rawdata.txt", sep = ""), sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
  write.table(dataset_s3, paste("/Users/boo/Desktop/temp/sub_", sub_num, "_s3_rawdata.txt", sep = ""), sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
  write.table(dataset_s4, paste("/Users/boo/Desktop/temp/sub_", sub_num, "_s4_rawdata.txt", sep = ""), sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
}





