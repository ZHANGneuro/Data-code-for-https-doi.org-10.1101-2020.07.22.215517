# loading packages
library(ggplot2)
rm(list=ls(all=TRUE)) 
se <- function(x){sd(x,na.rm=TRUE)/sqrt(length(x))} 
library(lme4)
library(lmerTest)

add_correct_marker <- function(dataset){
  
  counter = c()
  for(i in 1:length(dataset[,1])){
    
    if(dataset[i,2] == "sc"){
      if (grepl(substr(dataset[i,12], 1, 1), dataset[i,13])){
        counter = c(counter, 1)
      }else {
        counter = c(counter, 0)
      }
    }
    if(dataset[i,2] == "dc"){
      if(dataset[i,11]=="left"){
        if(grepl(substr(dataset[i,12], 4, 4), dataset[i,13])){
          counter = c(counter, 1)
        }else {
          counter = c(counter, 0)
        }
      }
      if(dataset[i,11]=="right"){
        if(grepl(substr(dataset[i,12], 2, 2), dataset[i,13])){
          counter = c(counter, 1)
        }else {
          counter = c(counter, 0)
        }
      }
      if(dataset[i,11]=="back"){
        if(grepl(substr(dataset[i,12], 3, 3), dataset[i,13])){
          counter = c(counter, 1)
        }else {
          counter = c(counter, 0)
        }
      }
    }
    if(dataset[i,2] == "fpc"){
      if(grepl(dataset[i,12], dataset[i,13])){
        counter = c(counter, 1)
      }else {
        counter = c(counter, 0)
      }
    }
    if(dataset[i,2] == "tpc"){
      if(grepl(dataset[i,12], dataset[i,13])){
        counter = c(counter, 1)
      }else {
        counter = c(counter, 0)
      }
    }
    if(dataset[i,2] == "hnd"){
      if(dataset[i,12] == "Correct"){
        counter = c(counter, 1)
      }else {
        counter = c(counter, 0)
      }
    }
  }
  
  dataset <- cbind(dataset,counter)
  return(dataset)
}


add_correct_marker_mri <- function(dataset, ith_correct, ith_choice, withQ2){
  if(withQ2){
    for(i in 1:length(dataset[,1])){
      if(dataset[i,ith_choice]=="Q2_left"){
        dataset[i,ith_choice]="left"
      }
      if(dataset[i,ith_choice]=="Q2_right"){
        dataset[i,ith_choice]="right"
      }
      if(dataset[i,ith_choice]=="Q2_back"){
        dataset[i,ith_choice]="back"
      }
    }
  }
  counter = c()
  for (i in 1:length(dataset[, 1])) {
    if (is.na(dataset[i, ith_correct])){
      if(dataset[i, 6]=="Incorrect"){
        counter = c(counter, 0)
      }
      if(dataset[i, 6]=="Correct"){
        counter = c(counter, 1)
      }
    } else {
      if (dataset[i, ith_correct] == dataset[i, ith_choice]) {
        counter = c(counter, 1)
      } else{
        counter = c(counter, 0)
      }
    }
  }
  dataset <- cbind(dataset,counter)
  return(dataset)
}


map_transformation <- function(dataset){
  if(dataset[1,1]==2){
    for(i in 1:length(dataset[,1])){
      if(dataset[i,2]==3){
        dataset[i,2]=4
      }
    }
  }else if(dataset[1,1]==3){
    for(i in 1:length(dataset[,1])){
      if(dataset[i,2]==2){
        dataset[i,2]=3
      }else if(dataset[i,2]==3){
        dataset[i,2]=5
      }
    }
  }else if(dataset[1,1]==4){
    for(i in 1:length(dataset[,1])){
      if(dataset[i,2]==1){
        dataset[i,2]=2
      }else if(dataset[i,2]==2){
        dataset[i,2]=4
      }else if(dataset[i,2]==3){
        dataset[i,2]=6
      }
    }
  }else if(dataset[1,1]==5){
    for(i in 1:length(dataset[,1])){
      if(dataset[i,2]==1){
        dataset[i,2]=3
      }else if(dataset[i,2]==2){
        dataset[i,2]=5
      }else if(dataset[i,2]==3){
        dataset[i,2]=6
      }
    }
  }else if(dataset[1,1]==6){
    for(i in 1:length(dataset[,1])){
      if(dataset[i,2]==1){
        dataset[i,2]=4
      }else if(dataset[i,2]==2){
        dataset[i,2]=5
      }else if(dataset[i,2]==3){
        dataset[i,2]=6
      }
    }
  }
  return(dataset)
}

correct_rate_map_each_subject <- function(t_dataset,ith_col){
  which_map <- unique(t_dataset[,2])
  which_map <- which_map[order(which_map)]
  map1_rate = length(which(t_dataset[,2]==which_map[1] & t_dataset[,ith_col]==1))/ length(which(t_dataset[,2]==which_map[1]))
  map2_rate = length(which(t_dataset[,2]==which_map[2] & t_dataset[,ith_col]==1))/ length(which(t_dataset[,2]==which_map[2]))
  map3_rate = length(which(t_dataset[,2]==which_map[3] & t_dataset[,ith_col]==1))/ length(which(t_dataset[,2]==which_map[3]))
  table <- data.frame(
    y = c(map1_rate, map2_rate, map3_rate),
    x = c(which_map[1],
          which_map[2],
          which_map[3]),
    trial_num = c(length(which(t_dataset[,2]==which_map[1] & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,2]==which_map[2] & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,2]==which_map[3] & t_dataset[,ith_col]==1))),
    total_trial_num = c(length(which(t_dataset[,2]==which_map[1])),
                        length(which(t_dataset[,2]==which_map[2])),
                        length(which(t_dataset[,2]==which_map[3]))))
  # table <- data.frame(
  #   y = c(map1_rate, map2_rate, map3_rate),
  #   x = c(which_map[1],
  #         which_map[2],
  #         which_map[3]),
  #   trial_num = c(length(which(t_dataset[,2]==which_map[1] & t_dataset[,ith_col]==1)),
  #                 length(which(t_dataset[,2]==which_map[2] & t_dataset[,ith_col]==1)),
  #                 length(which(t_dataset[,2]==which_map[3] & t_dataset[,ith_col]==1))),
  #   total_trial_num = c(length(which(t_dataset[,2]==which_map[1])),
  #                       length(which(t_dataset[,2]==which_map[2])),
  #                       length(which(t_dataset[,2]==which_map[3]))))
  table <- data.frame(
    y = c(map1_rate, map2_rate,map3_rate),
    x = c(which_map[1],
          which_map[2],
          which_map[3]),
    trial_num = c(length(which(t_dataset[,2]==which_map[1] & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,2]==which_map[2] & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,2]==which_map[3] & t_dataset[,ith_col]==1))),
    total_trial_num = c(length(which(t_dataset[,2]==which_map[1])),
                        length(which(t_dataset[,2]==which_map[2])),
                        length(which(t_dataset[,2]==which_map[3]))))
  return(table)
}




direction_transformation <- function(input_table){
  for(i in 1:length(input_table[,1])){
    if(input_table[i,2]==1 | input_table[i,2]==2){
      if(input_table[i,3]==1){
        input_table[i,3]=2
      } else if(input_table[i,3]==2){
        input_table[i,3]=3
      } else if(input_table[i,3]==3){
        input_table[i,3]=4
      } else if(input_table[i,3]==4){
        input_table[i,3]=5
      } else if(input_table[i,3]==5){
        input_table[i,3]=6
      } else if(input_table[i,3]==6){
        input_table[i,3]=7
      } else if(input_table[i,3]==7){
        input_table[i,3]=1
      }
    }
    if(input_table[i,2]==3 | input_table[i,2]==4){
      if(input_table[i,3]==1){
        input_table[i,3]=4
      } else if(input_table[i,3]==2){
        input_table[i,3]=5
      } else if(input_table[i,3]==3){
        input_table[i,3]=6
      } else if(input_table[i,3]==4){
        input_table[i,3]=7
      } else if(input_table[i,3]==5){
        input_table[i,3]=1
      } else if(input_table[i,3]==6){
        input_table[i,3]=2
      } else if(input_table[i,3]==7){
        input_table[i,3]=3
      }
    }
    if(input_table[i,2]==5 | input_table[i,2]==6){
      if(input_table[i,3]==1){
        input_table[i,3]=6
      } else if(input_table[i,3]==2){
        input_table[i,3]=7
      } else if(input_table[i,3]==3){
        input_table[i,3]=1
      } else if(input_table[i,3]==4){
        input_table[i,3]=2
      } else if(input_table[i,3]==5){
        input_table[i,3]=3
      } else if(input_table[i,3]==6){
        input_table[i,3]=4
      } else if(input_table[i,3]==7){
        input_table[i,3]=5
      }
    }
  }
  return(input_table)
}

correct_rate_direction_each_subject <- function(t_dataset, ith_col,num_of_direction){
  which_direction <- unique(t_dataset[,3])
  which_direction <- which_direction[order(which_direction)]
  direction1_rate = length(which(t_dataset[,3]==which_direction[1] & t_dataset[,ith_col]==1))/ length(which(t_dataset[,3]==which_direction[1]))
  direction2_rate = length(which(t_dataset[,3]==which_direction[2] & t_dataset[,ith_col]==1))/ length(which(t_dataset[,3]==which_direction[2]))
  direction3_rate = length(which(t_dataset[,3]==which_direction[3] & t_dataset[,ith_col]==1))/ length(which(t_dataset[,3]==which_direction[3]))
  direction4_rate = length(which(t_dataset[,3]==which_direction[4] & t_dataset[,ith_col]==1))/ length(which(t_dataset[,3]==which_direction[4]))
  direction5_rate = length(which(t_dataset[,3]==which_direction[5] & t_dataset[,ith_col]==1))/ length(which(t_dataset[,3]==which_direction[5]))
  direction6_rate = length(which(t_dataset[,3]==which_direction[6] & t_dataset[,ith_col]==1))/ length(which(t_dataset[,3]==which_direction[6]))
  if(num_of_direction==4){
    table <- data.frame(
      y = c(direction1_rate, direction2_rate, direction3_rate, direction4_rate),
      x = c(paste("  Enter dir",as.character(which_direction[1])),
            paste("  Enter dir",as.character(which_direction[2])),
            paste("  Enter dir",as.character(which_direction[3])),
            paste("  Enter dir",as.character(which_direction[4]))),
      trial_num = c(length(which(t_dataset[,3]==which_direction[1] & t_dataset[,ith_col]==1)),
                    length(which(t_dataset[,3]==which_direction[2] & t_dataset[,ith_col]==1)),
                    length(which(t_dataset[,3]==which_direction[3] & t_dataset[,ith_col]==1)),
                    length(which(t_dataset[,3]==which_direction[4] & t_dataset[,ith_col]==1))),
      total_trial_num = c(length(which(t_dataset[,3]==which_direction[1])),
                          length(which(t_dataset[,3]==which_direction[2])),
                          length(which(t_dataset[,3]==which_direction[3])),
                          length(which(t_dataset[,3]==which_direction[4]))))
  }
  return(table)
}



correct_rate_fov_each_subject <- function(dataset, ith_col){
  which_fov <- unique(dataset[,4])
  which_fov <- which_fov[order(which_fov)]
  fov1_rate = length(which(dataset[,4]==which_fov[1] & dataset[,ith_col]==1))/ length(which(dataset[,4]==which_fov[1]))
  fov2_rate = length(which(dataset[,4]==which_fov[2] & dataset[,ith_col]==1))/ length(which(dataset[,4]==which_fov[2]))
  table <- data.frame(
    y = c(fov1_rate, fov2_rate),
    x = c(which_fov[1],
          which_fov[2]),
    trial_num = c(length(which(dataset[,4]==which_fov[1] & dataset[,ith_col]==1)),
                  length(which(dataset[,4]==which_fov[2] & dataset[,ith_col]==1))),
    total_trial_num = c(length(which(dataset[,4]==which_fov[1])),
                        length(which(dataset[,4]==which_fov[2])))
  )
  return(table)
}
plot_correct_rate_fov <- function(t_dataset,num_table){   
  table <- data.frame(
    individual = c(t_dataset[,1],t_dataset[,2]),
    mean = c(rep(mean(t_dataset[,1]),length(t_dataset[,1])),
             rep(mean(t_dataset[,2]),length(t_dataset[,1]))),
    se = c(rep(std(t_dataset[,1]),length(t_dataset[,1])),
           rep(std(t_dataset[,2]),length(t_dataset[,1]))),
    trial_number = c(rep(mean(num_table[,1]),length(t_dataset[,1])),
                     rep(mean(num_table[,2]),length(t_dataset[,1]))),
    x_label = c(rep("fov_50",length(t_dataset[,1])),
                rep("fov_65",length(t_dataset[,1]))))
  table$x_label <- factor( table$x_label, levels=c("fov_50", "fov_65"))
  figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
    geom_bar(position=position_dodge(), stat="identity", width=.5) +
    geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.1)) +
    geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
    theme_classic() +    
    labs(x="",y="Performance")+ 
    ggtitle("field of view\n") + 
    theme(axis.text.x = element_text(size=20),axis.text.y =element_text(size=25),axis.title.y = element_text(size=25),plot.title = element_text(size=25)) +  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  return(figure)
}

correct_rate_facing_each_subject <- function(t_dataset, ith_col){
  facing1_rate = length(which(t_dataset[,8]==1 & t_dataset[,ith_col]==1))/ length(which(t_dataset[,8]==1))
  facing2_rate = length(which(t_dataset[,8]==2 & t_dataset[,ith_col]==1))/ length(which(t_dataset[,8]==2))
  facing3_rate = length(which(t_dataset[,8]==3 & t_dataset[,ith_col]==1))/ length(which(t_dataset[,8]==3))
  table <- data.frame(
    y = c(facing1_rate, facing2_rate, facing3_rate),
    x = c("Character 1","Character 2", "Character 3"),
    trial_num = c(length(which(t_dataset[,8]==1 & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,8]==2 & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,8]==3 & t_dataset[,ith_col]==1))),
    total_trial_num = c(length(which(t_dataset[,8]==1)),
                        length(which(t_dataset[,8]==2)),
                        length(which(t_dataset[,8]==3))))
  return(table)
}


correct_rate_facing_direction_each_subject <- function(t_dataset, ith_col){
  facing1_rate = length(which(t_dataset[,17]==1 & t_dataset[,ith_col]==1))/ length(which(t_dataset[,17]==1))
  facing2_rate = length(which(t_dataset[,17]==2 & t_dataset[,ith_col]==1))/ length(which(t_dataset[,17]==2))
  facing3_rate = length(which(t_dataset[,17]==3 & t_dataset[,ith_col]==1))/ length(which(t_dataset[,17]==3))
  table <- data.frame(
    y = c(facing1_rate, facing2_rate, facing3_rate),
    x = c("Facing direction 1","Facing direction 2", "Facing direction 3"),
    trial_num = c(length(which(t_dataset[,17]==1 & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,17]==2 & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,17]==3 & t_dataset[,ith_col]==1))),
    total_trial_num = c(length(which(t_dataset[,17]==1)),
                        length(which(t_dataset[,17]==2)),
                        length(which(t_dataset[,17]==3))))
  return(table)
}


correct_rate_targeting_each_subject <- function(t_dataset, ith_col){
  target1_rate = length(which(t_dataset[,11]==1 & t_dataset[,ith_col]==1))/ length(which(t_dataset[,11]==1))
  target2_rate = length(which(t_dataset[,11]==2 & t_dataset[,ith_col]==1))/ length(which(t_dataset[,11]==2))
  target3_rate = length(which(t_dataset[,11]==3 & t_dataset[,ith_col]==1))/ length(which(t_dataset[,11]==3))
  table <- data.frame(
    y = c(target1_rate, target2_rate, target3_rate),
    x = c("Character 1","Character 2", "Character 3"),
    trial_num = c(length(which(t_dataset[,11]==1 & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,11]==2 & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,11]==3 & t_dataset[,ith_col]==1))),
    total_trial_num = c(length(which(t_dataset[,11]==1)),
                        length(which(t_dataset[,11]==2)),
                        length(which(t_dataset[,11]==3))))
  return(table)
}

plot_correct_rate_facing <- function(t_dataset,num_table, title){   
  table <- data.frame(
    individual = c(t_dataset[,1],t_dataset[,2],t_dataset[,3]),
    mean = c(rep(mean(t_dataset[,1]),length(t_dataset[,1])),
             rep(mean(t_dataset[,2]),length(t_dataset[,1])),
             rep(mean(t_dataset[,3]),length(t_dataset[,1]))),
    se = c(rep(std(t_dataset[,1]),length(t_dataset[,1])),
           rep(std(t_dataset[,2]),length(t_dataset[,1])),
           rep(std(t_dataset[,3]),length(t_dataset[,1]))),
    trial_number = c(rep(mean(num_table[,1]),length(t_dataset[,1])),
                     rep(mean(num_table[,2]),length(t_dataset[,1])),
                     rep(mean(num_table[,3]),length(t_dataset[,1]))),
    x_label = c(rep("   person 1",length(t_dataset[,1])),
                rep("   person 2",length(t_dataset[,1])),
                rep("   person 3",length(t_dataset[,1]))))
  table$x_label <- factor( table$x_label, levels=c("   person 1", "   person 2", "   person 3"))
  figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
    geom_bar(position=position_dodge(), stat="identity", width=.26) +
    geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.1)) +
    geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
    theme_classic() +    
    labs(x="",y="Performance")+ 
    ggtitle(paste(title,"\n\n", sep = "")) +
    geom_hline(yintercept=0.33, linetype="dashed", color = "orange")+
    theme(axis.text.x = element_text(size=20),axis.text.y =element_text(size=25),axis.title.y = element_text(size=25),plot.title = element_text(size=25)) +  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) 
  return(figure)
}



plot_correct_rate_facing_direction <- function(rate_table,num_table, title){   
  table <- data.frame(
    individual = c(rate_table[,1],rate_table[,2],rate_table[,3]),
    mean = c(rep(mean(rate_table[,1]),length(rate_table[,1])),
             rep(mean(rate_table[,2]),length(rate_table[,1])),
             rep(mean(rate_table[,3]),length(rate_table[,1]))),
    se = c(rep(std(rate_table[,1]),length(rate_table[,1])),
           rep(std(rate_table[,2]),length(rate_table[,1])),
           rep(std(rate_table[,3]),length(rate_table[,1]))),
    trial_number = c(rep(mean(num_table[,1]),length(rate_table[,1])),
                     rep(mean(num_table[,2]),length(rate_table[,1])),
                     rep(mean(num_table[,3]),length(rate_table[,1]))),
    x_label = c(rep("Facing direction 1",length(rate_table[,1])),
                rep("Facing direction 2",length(rate_table[,1])),
                rep("Facing direction 3",length(rate_table[,1]))))
  table$x_label <- factor( table$x_label, levels=c("Facing direction 1", "Facing direction 2", "Facing direction 3"))
  figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
    geom_bar(position=position_dodge(), stat="identity", width=.26) +
    geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.1)) +
    geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
    theme_classic() +    
    labs(x="",y="Performance")+ 
    ggtitle(paste(title,"\n\n", sep = "")) +
    geom_hline(yintercept=0.33, linetype="dashed", color = "orange")+
    theme(axis.text.x = element_text(size=20),axis.text.y =element_text(size=25),axis.title.y = element_text(size=25),plot.title = element_text(size=25)) +  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) 
  return(figure)
}

plot_correct_rate_targeting <- function(rate_table,num_table, title){   
  table <- data.frame(
    individual = c(rate_table[,1],rate_table[,2],rate_table[,3]),
    mean = c(rep(mean(rate_table[,1]),length(rate_table[,1])),
             rep(mean(rate_table[,2]),length(rate_table[,1])),
             rep(mean(rate_table[,3]),length(rate_table[,1]))),
    se = c(rep(std(rate_table[,1]),length(rate_table[,1])),
           rep(std(rate_table[,2]),length(rate_table[,1])),
           rep(std(rate_table[,3]),length(rate_table[,1]))),
    trial_number = c(rep(mean(num_table[,1]),length(rate_table[,1])),
                     rep(mean(num_table[,2]),length(rate_table[,1])),
                     rep(mean(num_table[,3]),length(rate_table[,1]))),
    x_label = c(rep("Character 1",length(rate_table[,1])),
                rep("Character 2",length(rate_table[,1])),
                rep("Character 3",length(rate_table[,1]))))
  table$x_label <- factor( table$x_label, levels=c("Character 1", "Character 2", "Character 3"))
  figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
    geom_bar(position=position_dodge(), stat="identity", width=.26) +
    geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.1)) +
    geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
    theme_classic() +    
    labs(x="",y="Performance")+ 
    ggtitle(paste(title,"\n\n", sep = "")) +
    geom_hline(yintercept=0.33, linetype="dashed", color = "orange")+
    theme(axis.text.x = element_text(size=20),axis.text.y =element_text(size=25),axis.title.y = element_text(size=25),plot.title = element_text(size=25)) +  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) 
  return(figure)
}

correct_rate_ego_each_subject <- function(t_dataset, ith_col){
  ego1_rate = length(which(t_dataset[,10]=="left" & t_dataset[,ith_col]==1))/ length(which(t_dataset[,10]=="left"))
  ego2_rate = length(which(t_dataset[,10]=="right" & t_dataset[,ith_col]==1))/ length(which(t_dataset[,10]=="right"))
  ego3_rate = length(which(t_dataset[,10]=="back" & t_dataset[,ith_col]==1))/ length(which(t_dataset[,10]=="back"))
  table <- data.frame(
    y = c(ego1_rate, ego2_rate, ego3_rate),
    x = c("Left","Right","Back"),
    trial_num = c(length(which(t_dataset[,10]=="left" & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,10]=="right" & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,10]=="back" & t_dataset[,ith_col]==1))),
    total_trial_num = c(length(which(t_dataset[,10]=="left")),
                        length(which(t_dataset[,10]=="right")),
                        length(which(t_dataset[,10]=="back"))))
  return(table)
}
correct_rate_ego_each_subject_meg <- function(t_dataset, ith_col){
  ego1_rate = length(which(t_dataset[,11]=="left" & t_dataset[,ith_col]==1))/ length(which(t_dataset[,11]=="left"))
  ego2_rate = length(which(t_dataset[,11]=="right" & t_dataset[,ith_col]==1))/ length(which(t_dataset[,11]=="right"))
  ego3_rate = length(which(t_dataset[,11]=="back" & t_dataset[,ith_col]==1))/ length(which(t_dataset[,11]=="back"))
  table <- data.frame(
    y = c(ego1_rate, ego2_rate, ego3_rate),
    x = c("Left","Right","Back"),
    trial_num = c(length(which(t_dataset[,11]=="left" & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,11]=="right" & t_dataset[,ith_col]==1)),
                  length(which(t_dataset[,11]=="back" & t_dataset[,ith_col]==1))),
    total_trial_num = c(length(which(t_dataset[,11]=="left")),
                        length(which(t_dataset[,11]=="right")),
                        length(which(t_dataset[,11]=="back"))))
  return(table)
}

plot_correct_rate_ego <- function(rate_table,num_table,title){   
  table <- data.frame(
    individual = c(rate_table[,1],rate_table[,2],rate_table[,3]),
    mean = c(rep(mean(rate_table[,1]),length(rate_table[,1])),
             rep(mean(rate_table[,2]),length(rate_table[,1])),
             rep(mean(rate_table[,3]),length(rate_table[,1]))),
    se = c(rep(std(rate_table[,1]),length(rate_table[,1])),
           rep(std(rate_table[,2]),length(rate_table[,1])),
           rep(std(rate_table[,3]),length(rate_table[,1]))),
    trial_number = c(rep(mean(num_table[,1]),length(rate_table[,1])),
                     rep(mean(num_table[,2]),length(rate_table[,1])),
                     rep(mean(num_table[,3]),length(rate_table[,1]))),
    x_label = c(rep("          Left",length(rate_table[,1])),
                rep("         Right",length(rate_table[,1])),
                rep("          Back",length(rate_table[,1]))))
  table$x_label <- factor( table$x_label, levels=c("          Left", "         Right", "          Back"))
  figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
    geom_bar(position=position_dodge(), stat="identity", width=.26) +
    geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.1)) +
    geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
    theme_classic() +    
    labs(x="",y="Performance")+ 
    geom_hline(yintercept=0.33, linetype="dashed", color = "orange")+
    ggtitle(paste(title,"\n\n", sep = "")) +
    theme(axis.text.x = element_text(size=20),axis.text.y =element_text(size=25),axis.title.y = element_text(size=25),plot.title = element_text(size=25)) +  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) 
  return(figure)
}


add_absolute_location_facing_cha <- function(input_table){   
  counter = c()
  for (ith_row in 1:length(input_table[,1])) {
    head_trial_indicator = input_table[ith_row,5]
    ith_group = input_table[ith_row,1]
    ith_map = input_table[ith_row,2]
    ith_facing_cha = input_table[ith_row,8]
    if(is.na(head_trial_indicator)){
      if(    (ith_group==1 & ith_map == 1 & ith_facing_cha==1) | (ith_group==1 & ith_map == 2 & ith_facing_cha==1) | (ith_group==1 & ith_map == 3 & ith_facing_cha==2) |
             (ith_group==2 & ith_map == 1 & ith_facing_cha==1) | (ith_group==2 & ith_map == 2 & ith_facing_cha==1) | (ith_group==2 & ith_map == 3 & ith_facing_cha==3) |
             (ith_group==3 & ith_map == 1 & ith_facing_cha==1) | (ith_group==3 & ith_map == 2 & ith_facing_cha==2) | (ith_group==3 & ith_map == 3 & ith_facing_cha==2) |
             (ith_group==4 & ith_map == 1 & ith_facing_cha==1) | (ith_group==4 & ith_map == 2 & ith_facing_cha==3) | (ith_group==4 & ith_map == 3 & ith_facing_cha==3) |
             (ith_group==5 & ith_map == 1 & ith_facing_cha==2) | (ith_group==5 & ith_map == 2 & ith_facing_cha==2) | (ith_group==5 & ith_map == 3 & ith_facing_cha==3) |
             (ith_group==6 & ith_map == 1 & ith_facing_cha==3) | (ith_group==6 & ith_map == 2 & ith_facing_cha==2) | (ith_group==6 & ith_map == 3 & ith_facing_cha==3) 
      ){
        counter = c(counter, 1)
      }
      if(    (ith_group==1 & ith_map == 1 & ith_facing_cha==2) | (ith_group==1 & ith_map == 2 & ith_facing_cha==3) | (ith_group==1 & ith_map == 3 & ith_facing_cha==1) |
             (ith_group==2 & ith_map == 1 & ith_facing_cha==2) | (ith_group==2 & ith_map == 2 & ith_facing_cha==3) | (ith_group==2 & ith_map == 3 & ith_facing_cha==1) |
             (ith_group==3 & ith_map == 1 & ith_facing_cha==2) | (ith_group==3 & ith_map == 2 & ith_facing_cha==1) | (ith_group==3 & ith_map == 3 & ith_facing_cha==3) |
             (ith_group==4 & ith_map == 1 & ith_facing_cha==3) | (ith_group==4 & ith_map == 2 & ith_facing_cha==1) | (ith_group==4 & ith_map == 3 & ith_facing_cha==2) |
             (ith_group==5 & ith_map == 1 & ith_facing_cha==1) | (ith_group==5 & ith_map == 2 & ith_facing_cha==3) | (ith_group==5 & ith_map == 3 & ith_facing_cha==2) |
             (ith_group==6 & ith_map == 1 & ith_facing_cha==1) | (ith_group==6 & ith_map == 2 & ith_facing_cha==3) | (ith_group==6 & ith_map == 3 & ith_facing_cha==2) 
      ){
        counter = c(counter, 2)
      }
      if(    (ith_group==1 & ith_map == 1 & ith_facing_cha==3) | (ith_group==1 & ith_map == 2 & ith_facing_cha==2) | (ith_group==1 & ith_map == 3 & ith_facing_cha==3) |
             (ith_group==2 & ith_map == 1 & ith_facing_cha==3) | (ith_group==2 & ith_map == 2 & ith_facing_cha==2) | (ith_group==2 & ith_map == 3 & ith_facing_cha==2) |
             (ith_group==3 & ith_map == 1 & ith_facing_cha==3) | (ith_group==3 & ith_map == 2 & ith_facing_cha==3) | (ith_group==3 & ith_map == 3 & ith_facing_cha==1) |
             (ith_group==4 & ith_map == 1 & ith_facing_cha==2) | (ith_group==4 & ith_map == 2 & ith_facing_cha==2) | (ith_group==4 & ith_map == 3 & ith_facing_cha==1) |
             (ith_group==5 & ith_map == 1 & ith_facing_cha==3) | (ith_group==5 & ith_map == 2 & ith_facing_cha==1) | (ith_group==5 & ith_map == 3 & ith_facing_cha==1) |
             (ith_group==6 & ith_map == 1 & ith_facing_cha==2) | (ith_group==6 & ith_map == 2 & ith_facing_cha==1) | (ith_group==6 & ith_map == 3 & ith_facing_cha==1) 
      ){
        counter = c(counter, 3)
      }
      
    } else {
      counter = c(counter, NA)
    }
  }
  input_table = cbind(input_table, counter)
  return(input_table)
}


add_absolute_location_targeting_cha <- function(input_table){   
  counter = c()
  for (ith_row in 1:length(input_table[,1])) {
    head_trial_indicator = input_table[ith_row,5]
    ith_group = input_table[ith_row,1]
    ith_map = input_table[ith_row,2]
    ith_targeting_cha = input_table[ith_row,11]
    if(is.na(head_trial_indicator)){
      if(    (ith_group==1 & ith_map == 1 & ith_targeting_cha==1) | (ith_group==1 & ith_map == 2 & ith_targeting_cha==1) | (ith_group==1 & ith_map == 3 & ith_targeting_cha==2) |
             (ith_group==2 & ith_map == 1 & ith_targeting_cha==1) | (ith_group==2 & ith_map == 2 & ith_targeting_cha==1) | (ith_group==2 & ith_map == 3 & ith_targeting_cha==3) |
             (ith_group==3 & ith_map == 1 & ith_targeting_cha==1) | (ith_group==3 & ith_map == 2 & ith_targeting_cha==2) | (ith_group==3 & ith_map == 3 & ith_targeting_cha==2) |
             (ith_group==4 & ith_map == 1 & ith_targeting_cha==1) | (ith_group==4 & ith_map == 2 & ith_targeting_cha==3) | (ith_group==4 & ith_map == 3 & ith_targeting_cha==3) |
             (ith_group==5 & ith_map == 1 & ith_targeting_cha==2) | (ith_group==5 & ith_map == 2 & ith_targeting_cha==2) | (ith_group==5 & ith_map == 3 & ith_targeting_cha==3) |
             (ith_group==6 & ith_map == 1 & ith_targeting_cha==3) | (ith_group==6 & ith_map == 2 & ith_targeting_cha==2) | (ith_group==6 & ith_map == 3 & ith_targeting_cha==3) 
      ){
        counter = c(counter, 1)
      }
      if(    (ith_group==1 & ith_map == 1 & ith_targeting_cha==2) | (ith_group==1 & ith_map == 2 & ith_targeting_cha==3) | (ith_group==1 & ith_map == 3 & ith_targeting_cha==1) |
             (ith_group==2 & ith_map == 1 & ith_targeting_cha==2) | (ith_group==2 & ith_map == 2 & ith_targeting_cha==3) | (ith_group==2 & ith_map == 3 & ith_targeting_cha==1) |
             (ith_group==3 & ith_map == 1 & ith_targeting_cha==2) | (ith_group==3 & ith_map == 2 & ith_targeting_cha==1) | (ith_group==3 & ith_map == 3 & ith_targeting_cha==3) |
             (ith_group==4 & ith_map == 1 & ith_targeting_cha==3) | (ith_group==4 & ith_map == 2 & ith_targeting_cha==1) | (ith_group==4 & ith_map == 3 & ith_targeting_cha==2) |
             (ith_group==5 & ith_map == 1 & ith_targeting_cha==1) | (ith_group==5 & ith_map == 2 & ith_targeting_cha==3) | (ith_group==5 & ith_map == 3 & ith_targeting_cha==2) |
             (ith_group==6 & ith_map == 1 & ith_targeting_cha==1) | (ith_group==6 & ith_map == 2 & ith_targeting_cha==3) | (ith_group==6 & ith_map == 3 & ith_targeting_cha==2) 
      ){
        counter = c(counter, 2)
      }
      if(    (ith_group==1 & ith_map == 1 & ith_targeting_cha==3) | (ith_group==1 & ith_map == 2 & ith_targeting_cha==2) | (ith_group==1 & ith_map == 3 & ith_targeting_cha==3) |
             (ith_group==2 & ith_map == 1 & ith_targeting_cha==3) | (ith_group==2 & ith_map == 2 & ith_targeting_cha==2) | (ith_group==2 & ith_map == 3 & ith_targeting_cha==2) |
             (ith_group==3 & ith_map == 1 & ith_targeting_cha==3) | (ith_group==3 & ith_map == 2 & ith_targeting_cha==3) | (ith_group==3 & ith_map == 3 & ith_targeting_cha==1) |
             (ith_group==4 & ith_map == 1 & ith_targeting_cha==2) | (ith_group==4 & ith_map == 2 & ith_targeting_cha==2) | (ith_group==4 & ith_map == 3 & ith_targeting_cha==1) |
             (ith_group==5 & ith_map == 1 & ith_targeting_cha==3) | (ith_group==5 & ith_map == 2 & ith_targeting_cha==1) | (ith_group==5 & ith_map == 3 & ith_targeting_cha==1) |
             (ith_group==6 & ith_map == 1 & ith_targeting_cha==2) | (ith_group==6 & ith_map == 2 & ith_targeting_cha==1) | (ith_group==6 & ith_map == 3 & ith_targeting_cha==1) 
      ){
        counter = c(counter, 3)
      }
    } else {
      counter = c(counter, NA)
    }
  }
  input_table = cbind(input_table, counter)
  return(input_table)
}


angle <- function(input_map, input_enterDirection,  input_facingDirection){
  if(  ((input_map == "1" | input_map == "2") & input_enterDirection == "2"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "3"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "4"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "5"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "6"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "7"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "4"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "5"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "6"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "7"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "2"  & input_facingDirection == "facing_direction_4") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "3"  & input_facingDirection == "facing_direction_4") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "6"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "7"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "2"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "3"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "4"  & input_facingDirection == "facing_direction_4") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "5"  & input_facingDirection == "facing_direction_4")){
    return(157);
  }
  if(  ((input_map == "1" | input_map == "2") & input_enterDirection == "3"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "4"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "5"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "6"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "4"  & input_facingDirection == "facing_direction_4") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "5"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "6"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "3"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "3"  & input_facingDirection == "facing_direction_4") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "4"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "5"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "6"  & input_facingDirection == "facing_direction_4")){
    return(112);
  }
  if(  ((input_map == "1" | input_map == "2") & input_enterDirection == "1"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "1"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "1"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "1"  & input_facingDirection == "facing_direction_4") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "1"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "1"  & input_facingDirection == "facing_direction_3")){
    return(90);
  }
  if(  ((input_map == "1" | input_map == "2") & input_enterDirection == "2"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "4"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "5"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "7"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "4"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "5"  & input_facingDirection == "facing_direction_4") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "7"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "2"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "2"  & input_facingDirection == "facing_direction_4") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "4"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "5"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "7"  & input_facingDirection == "facing_direction_4")){
    return(68);
  }
  if(  ((input_map == "1" | input_map == "2") & input_enterDirection == "2"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "3"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "6"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "1" | input_map == "2") & input_enterDirection == "7"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "6"  & input_facingDirection == "facing_direction_4") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "7"  & input_facingDirection == "facing_direction_4") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "2"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "3"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "2"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "3"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "6"  & input_facingDirection == "facing_direction_3") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "7"  & input_facingDirection == "facing_direction_3")){
    return(23);
  }
  if(  ((input_map == "1" | input_map == "2") & input_enterDirection == "1"  & input_facingDirection == "facing_direction_2") |
       ((input_map == "3" | input_map == "4") & input_enterDirection == "1"  & input_facingDirection == "facing_direction_1") |
       ((input_map == "5" | input_map == "6") & input_enterDirection == "1"  & input_facingDirection == "facing_direction_4")){
    return(0);
  }
}

angle_new <- function(input_map, input_enterDirection,  input_facing){
  if(  (input_map == "1" & input_enterDirection == "1"  & input_facing == "3") |
       (input_map == "1" & input_enterDirection == "3"  & input_facing == "1") |
       (input_map == "1" & input_enterDirection == "4"  & input_facing == "2") |
       (input_map == "2" & input_enterDirection == "1"  & input_facing == "2") |
       (input_map == "2" & input_enterDirection == "3"  & input_facing == "1") |
       (input_map == "2" & input_enterDirection == "4"  & input_facing == "2") |
       (input_map == "3" & input_enterDirection == "1"  & input_facing == "3") |
       (input_map == "3" & input_enterDirection == "3"  & input_facing == "2") |
       (input_map == "3" & input_enterDirection == "4"  & input_facing == "1") |
       (input_map == "4" & input_enterDirection == "1"  & input_facing == "2") |
       (input_map == "4" & input_enterDirection == "3"  & input_facing == "3") |
       (input_map == "4" & input_enterDirection == "4"  & input_facing == "1") |
       (input_map == "5" & input_enterDirection == "1"  & input_facing == "1") |
       (input_map == "5" & input_enterDirection == "3"  & input_facing == "2") |
       (input_map == "5" & input_enterDirection == "4"  & input_facing == "3") |
       (input_map == "6" & input_enterDirection == "1"  & input_facing == "1") |
       (input_map == "6" & input_enterDirection == "3"  & input_facing == "3") |
       (input_map == "6" & input_enterDirection == "4"  & input_facing == "2")){
    return("left_45°")  
  } else if ((input_map == "1" & input_enterDirection == "1"  & input_facing == "2") |
             (input_map == "1" & input_enterDirection == "2"  & input_facing == "3") |
             (input_map == "1" & input_enterDirection == "4"  & input_facing == "1") |
             (input_map == "2" & input_enterDirection == "2"  & input_facing == "2") |
             (input_map == "2" & input_enterDirection == "4"  & input_facing == "1") |
             (input_map == "2" & input_enterDirection == "1"  & input_facing == "3") |
             (input_map == "3" & input_enterDirection == "1"  & input_facing == "1") |
             (input_map == "3" & input_enterDirection == "2"  & input_facing == "3") |
             (input_map == "3" & input_enterDirection == "4"  & input_facing == "2") |
             (input_map == "4" & input_enterDirection == "1"  & input_facing == "1") |
             (input_map == "4" & input_enterDirection == "2"  & input_facing == "2") |
             (input_map == "4" & input_enterDirection == "4"  & input_facing == "3") |
             (input_map == "5" & input_enterDirection == "1"  & input_facing == "3") |
             (input_map == "5" & input_enterDirection == "2"  & input_facing == "1") |
             (input_map == "5" & input_enterDirection == "4"  & input_facing == "2") |
             (input_map == "6" & input_enterDirection == "1"  & input_facing == "2") |
             (input_map == "6" & input_enterDirection == "2"  & input_facing == "1") |
             (input_map == "6" & input_enterDirection == "4"  & input_facing == "3")) {
    return("right_45°")  
  } else if ((input_map == "1" & input_enterDirection == "2"  & input_facing == "1") |
             (input_map == "1" & input_enterDirection == "3"  & input_facing == "2") |
             (input_map == "1" & input_enterDirection == "4"  & input_facing == "3") |
             (input_map == "2" & input_enterDirection == "2"  & input_facing == "1") |
             (input_map == "2" & input_enterDirection == "3"  & input_facing == "3") |
             (input_map == "2" & input_enterDirection == "4"  & input_facing == "3") |
             (input_map == "3" & input_enterDirection == "2"  & input_facing == "2") |
             (input_map == "3" & input_enterDirection == "3"  & input_facing == "1") |
             (input_map == "3" & input_enterDirection == "4"  & input_facing == "3") |
             (input_map == "4" & input_enterDirection == "2"  & input_facing == "3") |
             (input_map == "4" & input_enterDirection == "3"  & input_facing == "1") |
             (input_map == "4" & input_enterDirection == "4"  & input_facing == "2") |
             (input_map == "5" & input_enterDirection == "2"  & input_facing == "2") |
             (input_map == "5" & input_enterDirection == "3"  & input_facing == "3") |
             (input_map == "5" & input_enterDirection == "4"  & input_facing == "1") |
             (input_map == "6" & input_enterDirection == "2"  & input_facing == "3") |
             (input_map == "6" & input_enterDirection == "3"  & input_facing == "2") |
             (input_map == "6" & input_enterDirection == "4"  & input_facing == "1")) {
    return("left_135°")  
  } else if ((input_map == "1" & input_enterDirection == "1"  & input_facing == "1") |
             (input_map == "1" & input_enterDirection == "2"  & input_facing == "2") |
             (input_map == "1" & input_enterDirection == "3"  & input_facing == "3") |
             (input_map == "2" & input_enterDirection == "1"  & input_facing == "1") |
             (input_map == "2" & input_enterDirection == "2"  & input_facing == "3") |
             (input_map == "2" & input_enterDirection == "3"  & input_facing == "2") |
             (input_map == "3" & input_enterDirection == "1"  & input_facing == "2") |
             (input_map == "3" & input_enterDirection == "2"  & input_facing == "1") |
             (input_map == "3" & input_enterDirection == "3"  & input_facing == "3") |
             (input_map == "4" & input_enterDirection == "1"  & input_facing == "3") |
             (input_map == "4" & input_enterDirection == "2"  & input_facing == "1") |
             (input_map == "4" & input_enterDirection == "3"  & input_facing == "2") |
             (input_map == "5" & input_enterDirection == "1"  & input_facing == "2") |
             (input_map == "5" & input_enterDirection == "2"  & input_facing == "3") |
             (input_map == "5" & input_enterDirection == "3"  & input_facing == "1") |
             (input_map == "6" & input_enterDirection == "1"  & input_facing == "3") |
             (input_map == "6" & input_enterDirection == "2"  & input_facing == "2") |
             (input_map == "6" & input_enterDirection == "3"  & input_facing == "1")) {
    return("right_135°")  
  }
}

add_finger_index<- function(dataset){
  hand = c()
  for (i in 1:length(dataset[, 1])) {
    if(is.na(dataset[i,10])){
      if(dataset[i,5]=="Q1_right"){
        hand = c(hand,0)
      }
      if(dataset[i,5]=="Q1_left"){
        hand = c(hand,1)
      }
    } else {
      if(dataset[i,10]=="left"){
        if(substr(dataset[i,12],1,1)==1){
          hand = c(hand,1)
        } else {
          hand = c(hand,0)
        }
      }
      if(dataset[i,10]=="right"){
        if(substr(dataset[i,12],3,3)==1){
          hand = c(hand,1)
        } else {
          hand = c(hand,0)
        }
      }
      if(dataset[i,10]=="back"){
        if(substr(dataset[i,12],2,2)==1){
          hand = c(hand,1)
        } else {
          hand = c(hand,0)
        }
      }
    }
  }
  dataset <- cbind(dataset,hand)
  return(dataset)
}

add_rotation_marker_new <- function(t_dataset){
  rotate_angle = c()
  for (i in 1:length(t_dataset[, 1])) {
    if(!is.na(t_dataset[i,8])){
      rotate_angle = c(rotate_angle,angle_new(t_dataset[i,2],t_dataset[i,3],t_dataset[i,8]))
    } else{
      rotate_angle = c(rotate_angle,NA)
    }
  }
  t_dataset <- cbind(t_dataset,rotate_angle)
  return(t_dataset)
}

add_rotation_marker <- function(t_dataset){
  rotate_angle = c()
  for (i in 1:length(t_dataset[, 1])) {
    # rotate_angle = c(rotate_angle,angle_new(t_dataset[i,2],t_dataset[i,3],t_dataset[i,***(facing_direction)]))
  }
  t_dataset <- cbind(t_dataset,rotate_angle)
  return(t_dataset)
}

correct_rate_rotation_each_subject <- function(t_dataset) {  
  index_col <- grep("rotate_angle", colnames(t_dataset))
  R_0_dataset = t_dataset[which(t_dataset[,index_col]== 0),]
  R_23_dataset = t_dataset[which(t_dataset[,index_col]== 23),]
  R_68_dataset = t_dataset[which(t_dataset[,index_col]== 68),]
  R_90_dataset = t_dataset[which(t_dataset[,index_col]== 90),]
  R_112_dataset = t_dataset[which(t_dataset[,index_col]== 112),]
  R_157_dataset = t_dataset[which(t_dataset[,index_col]== 157),]
  table <- data.frame(
    total_trial_num = c(length(R_0_dataset[,index_col]),
                        length(R_23_dataset[,index_col]),
                        length(R_68_dataset[,index_col]),
                        length(R_90_dataset[,index_col]),
                        length(R_112_dataset[,index_col]),
                        length(R_157_dataset[,index_col])),
    y  =  c(length(which(R_0_dataset[,13]==1))/length(R_0_dataset[,13]),
            length(which(R_23_dataset[,13]==1))/length(R_23_dataset[,13]),
            length(which(R_68_dataset[,13]==1))/length(R_68_dataset[,13]),
            length(which(R_90_dataset[,13]==1))/length(R_90_dataset[,13]),
            length(which(R_112_dataset[,13]==1))/length(R_112_dataset[,13]),
            length(which(R_157_dataset[,13]==1))/length(R_157_dataset[,13])),
    trial_num = c(length(which(R_0_dataset[,13]==1)),
                  length(which(R_23_dataset[,13]==1)),
                  length(which(R_68_dataset[,13]==1)),
                  length(which(R_90_dataset[,13]==1)),
                  length(which(R_112_dataset[,13]==1)),
                  length(which(R_157_dataset[,13]==1))))
  return(table)
}


correct_rate_rotation_each_subject_new <- function(dataset, correctness_ith) {  
  index_col <- grep("rotate_angle", colnames(dataset))
  left_45_dataset = dataset[which(dataset[,index_col]== "left_45°"),]
  right_45_dataset = dataset[which(dataset[,index_col]== "right_45°"),]
  left_135_dataset = dataset[which(dataset[,index_col]== "left_135°"),]
  right_135_dataset = dataset[which(dataset[,index_col]== "right_135°"),]
  
  table <- data.frame(
    total_trial_num = c(length(left_45_dataset[,index_col]),
                        length(right_45_dataset[,index_col]),
                        length(left_135_dataset[,index_col]),
                        length(right_135_dataset[,index_col])),
    y  =  c(length(which(left_45_dataset[,correctness_ith]==1))/length(left_45_dataset[,correctness_ith]),
            length(which(right_45_dataset[,correctness_ith]==1))/length(right_45_dataset[,correctness_ith]),
            length(which(left_135_dataset[,correctness_ith]==1))/length(left_135_dataset[,correctness_ith]),
            length(which(right_135_dataset[,correctness_ith]==1))/length(right_135_dataset[,correctness_ith])),
    trial_num = c(length(which(left_45_dataset[,correctness_ith]==1)),
                  length(which(right_45_dataset[,correctness_ith]==1)),
                  length(which(left_135_dataset[,correctness_ith]==1)),
                  length(which(right_135_dataset[,correctness_ith]==1))))
  return(table)
}


plot_correct_rate_rotation <- function(rate_table,num_table){   
  table <- data.frame(
    individual = c(rate_table[,1],rate_table[,2],rate_table[,3],rate_table[,4],rate_table[,5],rate_table[,6]),
    mean = c(rep(mean(rate_table[,1]),length(rate_table[,1])),
             rep(mean(rate_table[,2]),length(rate_table[,1])),
             rep(mean(rate_table[,3]),length(rate_table[,1])),
             rep(mean(rate_table[,4]),length(rate_table[,1])),
             rep(mean(rate_table[,5]),length(rate_table[,1])),
             rep(mean(rate_table[,6]),length(rate_table[,1]))),
    se = c(rep(std(rate_table[,1]),length(rate_table[,1])),
           rep(std(rate_table[,2]),length(rate_table[,1])),
           rep(std(rate_table[,3]),length(rate_table[,1])),
           rep(std(rate_table[,4]),length(rate_table[,1])),
           rep(std(rate_table[,5]),length(rate_table[,1])),
           rep(std(rate_table[,6]),length(rate_table[,1]))),
    trial_number = c(rep(mean(num_table[,1]),length(rate_table[,1])),
                     rep(mean(num_table[,2]),length(rate_table[,1])),
                     rep(mean(num_table[,3]),length(rate_table[,1])),
                     rep(mean(num_table[,4]),length(rate_table[,1])),
                     rep(mean(num_table[,5]),length(rate_table[,1])),
                     rep(mean(num_table[,6]),length(rate_table[,1]))),
    x_label = c(rep("0 deg",length(rate_table[,1])), 
                rep("23 deg",length(rate_table[,1])),
                rep("68 deg",length(rate_table[,1])), 
                rep("90 deg",length(rate_table[,1])), 
                rep("112 deg",length(rate_table[,1])), 
                rep("157 deg",length(rate_table[,1]))))    
  
  table$x_label <- factor( table$x_label, levels=c("0 deg", "23 deg", "68 deg", "90 deg", "112 deg", "157 deg"))
  figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
    geom_bar(position=position_dodge(), stat="identity", width=.5) +
    geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.1)) +
    geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
    theme_classic() +    
    labs(x="",y="Performance")+ 
    ggtitle(paste("Angle","\n\n", sep = "")) +
    theme(axis.text.x = element_text(size=20),axis.text.y =element_text(size=25),axis.title.y = element_text(size=25),plot.title = element_text(size=25)) +  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) 
  return(figure)
}


plot_correct_rate_rotation_new <- function(rate_table,num_table){   
  table <- data.frame(
    individual = c(rate_table[,1],rate_table[,2],rate_table[,3],rate_table[,4]),
    mean = c(rep(mean(rate_table[,1]),length(rate_table[,1])),
             rep(mean(rate_table[,2]),length(rate_table[,1])),
             rep(mean(rate_table[,3]),length(rate_table[,1])),
             rep(mean(rate_table[,4]),length(rate_table[,1]))),
    se = c(rep(std(rate_table[,1]),length(rate_table[,1])),
           rep(std(rate_table[,2]),length(rate_table[,1])),
           rep(std(rate_table[,3]),length(rate_table[,1])),
           rep(std(rate_table[,4]),length(rate_table[,1]))),
    trial_number = c(rep(mean(num_table[,1]),length(rate_table[,1])),
                     rep(mean(num_table[,2]),length(rate_table[,1])),
                     rep(mean(num_table[,3]),length(rate_table[,1])),
                     rep(mean(num_table[,4]),length(rate_table[,1]))),
    x_label = c(rep("L 45 deg", length(rate_table[,1])),
                rep("R 45 deg",length(rate_table[,1])),
                rep("L 135 deg", length(rate_table[,1])),
                rep("R 135 deg",length(rate_table[,1]))))    
  
  table$x_label <- factor( table$x_label, levels=c("L 45 deg", "R 45 deg", "L 135 deg", "R 135 deg"))
  figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
    geom_bar(position=position_dodge(), stat="identity", width=.5) +
    geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.1)) +
    geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
    theme_classic() +    
    labs(x="",y="Performance")+ 
    ggtitle(paste("Angle","\n\n", sep = "")) +
    geom_hline(yintercept=0.33, linetype="dashed", color = "orange")+
    theme(axis.text.x = element_text(size=20),axis.text.y =element_text(size=25),axis.title.y = element_text(size=25),plot.title = element_text(size=25)) +  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) 
  return(figure)
}





























training_map_transformation <- function(input_table){
  if(input_table[1,1]==2){
    for(i in 1:length(input_table[,1])){
      if(input_table[i,4]==3){
        input_table[i,4]=4
      }
    }
  }else if(input_table[1,1]==3){
    for(i in 1:length(input_table[,1])){
      if(input_table[i,4]==2){
        input_table[i,4]=3
      }else if(input_table[i,4]==3){
        input_table[i,4]=5
      }
    }
  }else if(input_table[1,1]==4){
    for(i in 1:length(input_table[,1])){
      if(input_table[i,4]==1){
        input_table[i,4]=2
      }else if(input_table[i,4]==2){
        input_table[i,4]=4
      }else if(input_table[i,4]==3){
        input_table[i,4]=6
      }
    }
  }else if(input_table[1,1]==5){
    for(i in 1:length(input_table[,1])){
      if(input_table[i,4]==1){
        input_table[i,4]=3
      }else if(input_table[i,4]==2){
        input_table[i,4]=5
      }else if(input_table[i,4]==3){
        input_table[i,4]=6
      }
    }
  }else if(input_table[1,1]==6){
    for(i in 1:length(input_table[,1])){
      if(input_table[i,4]==1){
        input_table[i,4]=4
      }else if(input_table[i,4]==2){
        input_table[i,4]=5
      }else if(input_table[i,4]==3){
        input_table[i,4]=6
      }
    }
  }
  return(input_table)
}

training_correct_rate_map_each_subject <- function(data_frame){
  which_map <- unique(data_frame[,4])
  which_map <- which_map[order(which_map)]
  map1_rate = length(which(data_frame[,4]==which_map[1] & data_frame[,7]=="Correct"))/ length(which(data_frame[,4]==which_map[1]))
  map2_rate = length(which(data_frame[,4]==which_map[2] & data_frame[,7]=="Correct"))/ length(which(data_frame[,4]==which_map[2]))
  map3_rate = length(which(data_frame[,4]==which_map[3] & data_frame[,7]=="Correct"))/ length(which(data_frame[,4]==which_map[3]))
  table <- data.frame(
    y = c(map1_rate, map2_rate, map3_rate),
    x = c(which_map[1],
          which_map[2],
          which_map[3]),
    trial_num = c(length(which(data_frame[,4]==which_map[1] & data_frame[,7]=="Correct")),
                  length(which(data_frame[,4]==which_map[2] & data_frame[,7]=="Correct")),
                  length(which(data_frame[,4]==which_map[3] & data_frame[,7]=="Correct"))),
    total_trial_num = c(length(which(data_frame[,4]==which_map[1])),
                        length(which(data_frame[,4]==which_map[2])),
                        length(which(data_frame[,4]==which_map[3]))))
  return(table)
}

training_direction_transformation <- function(input_table){
  for(i in 1:length(input_table[,1])){
    if(input_table[i,4]==1 | input_table[i,4]==2){
      if(input_table[i,5]==1){
        input_table[i,5]=2
      } else if(input_table[i,5]==2){
        input_table[i,5]=3
      } else if(input_table[i,5]==3){
        input_table[i,5]=4
      } else if(input_table[i,5]==4){
        input_table[i,5]=5
      } else if(input_table[i,5]==5){
        input_table[i,5]=6
      } else if(input_table[i,5]==6){
        input_table[i,5]=7
      } else if(input_table[i,5]==7){
        input_table[i,5]=1
      }
    }
    if(input_table[i,4]==3 | input_table[i,4]==4){
      if(input_table[i,5]==1){
        input_table[i,5]=4
      } else if(input_table[i,5]==2){
        input_table[i,5]=5
      } else if(input_table[i,5]==3){
        input_table[i,5]=6
      } else if(input_table[i,5]==4){
        input_table[i,5]=7
      } else if(input_table[i,5]==5){
        input_table[i,5]=1
      } else if(input_table[i,5]==6){
        input_table[i,5]=2
      } else if(input_table[i,5]==7){
        input_table[i,5]=3
      }
    }
    if(input_table[i,4]==5 | input_table[i,4]==6){
      if(input_table[i,5]==1){
        input_table[i,5]=6
      } else if(input_table[i,5]==2){
        input_table[i,5]=7
      } else if(input_table[i,5]==3){
        input_table[i,5]=1
      } else if(input_table[i,5]==4){
        input_table[i,5]=2
      } else if(input_table[i,5]==5){
        input_table[i,5]=3
      } else if(input_table[i,5]==6){
        input_table[i,5]=4
      } else if(input_table[i,5]==7){
        input_table[i,5]=5
      }
    }
  }
  return(input_table)
}

training_correct_rate_direction_each_subject <- function(data_frame){
  which_direction <- unique(data_frame[,5])
  which_direction <- which_direction[order(which_direction)]
  direction1_rate = length(which(data_frame[,5]==which_direction[1] & data_frame[,7]=="Correct"))/ length(which(data_frame[,5]==which_direction[1]))
  direction2_rate = length(which(data_frame[,5]==which_direction[2] & data_frame[,7]=="Correct"))/ length(which(data_frame[,5]==which_direction[2]))
  direction3_rate = length(which(data_frame[,5]==which_direction[3] & data_frame[,7]=="Correct"))/ length(which(data_frame[,5]==which_direction[3]))
  direction4_rate = length(which(data_frame[,5]==which_direction[4] & data_frame[,7]=="Correct"))/ length(which(data_frame[,5]==which_direction[4]))
  table <- data.frame(
    y = c(direction1_rate, direction2_rate, direction3_rate, direction4_rate),
    x = c(paste("Direction",as.character(which_direction[1])),
          paste("Direction",as.character(which_direction[2])),
          paste("Direction",as.character(which_direction[3])),
          paste("Direction",as.character(which_direction[4]))),
    trial_num = c(length(which(data_frame[,5]==which_direction[1] & data_frame[,7]=="Correct")),
                  length(which(data_frame[,5]==which_direction[2] & data_frame[,7]=="Correct")),
                  length(which(data_frame[,5]==which_direction[3] & data_frame[,7]=="Correct")),
                  length(which(data_frame[,5]==which_direction[4] & data_frame[,7]=="Correct"))),
    total_trial_num = c(length(which(data_frame[,5]==which_direction[1])),
                        length(which(data_frame[,5]==which_direction[2])),
                        length(which(data_frame[,5]==which_direction[3])),
                        length(which(data_frame[,5]==which_direction[4]))))
  return(table)
}

meg_false_alarm_rate_each_subject <- function(dataset_hnd) {   
  
  hn_ith = unique(which(grepl("_1",dataset_hnd[,7]) | grepl("_1",dataset_hnd[,8]) | grepl("_1",dataset_hnd[,9])))
  no_hn_ith = unique(which(!grepl("_1",dataset_hnd[,7]) & !grepl("_1",dataset_hnd[,8]) & !grepl("_1",dataset_hnd[,9])))
  
  hit_rate_count <- length(which(dataset_hnd[,12]=="Correct" & (grepl("_1",dataset_hnd[,7]) | grepl("_1",dataset_hnd[,8]) | grepl("_1",dataset_hnd[,9]))))
  miss_rate_count <- length(which(dataset_hnd[,12]=="Incorrect" & (grepl("_1",dataset_hnd[,7]) | grepl("_1",dataset_hnd[,8]) | grepl("_1",dataset_hnd[,9]))))
  False_alarm_rate_count <- length(which(dataset_hnd[,12]=="Incorrect" & (!grepl("_1",dataset_hnd[,7]) & !grepl("_1",dataset_hnd[,8]) & !grepl("_1",dataset_hnd[,9]))))
  correct_rejection_rate_count <- length(which(dataset_hnd[,12]=="Correct" & (!grepl("_1",dataset_hnd[,7]) & !grepl("_1",dataset_hnd[,8]) & !grepl("_1",dataset_hnd[,9]))))
  
  hit_rate = hit_rate_count/length(hn_ith)
  miss_rate  = miss_rate_count/length(hn_ith)
  false_alarm_rate  = False_alarm_rate_count/length(no_hn_ith)
  correct_rejection_rate = correct_rejection_rate_count/length(no_hn_ith)
  
  table <- data.frame(
    y = c(hit_rate, miss_rate, false_alarm_rate, correct_rejection_rate),
    x = c("H","M","FA","CR"),
    trial_num = c(hit_rate_count,miss_rate_count,False_alarm_rate_count,correct_rejection_rate_count),
    total_trial_num = c(length(hn_ith),length(hn_ith),length(no_hn_ith),length(no_hn_ith)))
  table$x <- factor(table$x, levels=c("CR","FA","H","M"))
  return(table)
}




training_false_alarm_rate_each_subject <- function(dataset_hnd) {   
  
  hn_ith = unique(which(grepl("_1",dataset_hnd[,4]) | grepl("_1",dataset_hnd[,5]) | grepl("_1",dataset_hnd[,6])))
  no_hn_ith = unique(which(!grepl("_1",dataset_hnd[,4]) & !grepl("_1",dataset_hnd[,5]) & !grepl("_1",dataset_hnd[,6])))
  
  hit_rate_count <- length(which(dataset_hnd[,9]=="Correct" & (grepl("_1",dataset_hnd[,4]) | grepl("_1",dataset_hnd[,5]) | grepl("_1",dataset_hnd[,6]))))
  miss_rate_count <- length(which(dataset_hnd[,9]=="Incorrect" & (grepl("_1",dataset_hnd[,4]) | grepl("_1",dataset_hnd[,5]) | grepl("_1",dataset_hnd[,6]))))
  False_alarm_rate_count <- length(which(dataset_hnd[9]=="Incorrect" & (!grepl("_1",dataset_hnd[,4]) & !grepl("_1",dataset_hnd[,5]) & !grepl("_1",dataset_hnd[,6]))))
  correct_rejection_rate_count <- length(which(dataset_hnd[,9]=="Correct" & (!grepl("_1",dataset_hnd[,4]) & !grepl("_1",dataset_hnd[,5]) & !grepl("_1",dataset_hnd[,6]))))
  
  hit_rate = hit_rate_count/length(hn_ith)
  miss_rate  = miss_rate_count/length(hn_ith)
  false_alarm_rate  = False_alarm_rate_count/length(no_hn_ith)
  correct_rejection_rate = correct_rejection_rate_count/length(no_hn_ith)
  
  table <- data.frame(
    y = c(hit_rate, miss_rate, false_alarm_rate, correct_rejection_rate),
    x = c("H","M","FA","CR"),
    trial_num = c(hit_rate_count,miss_rate_count,False_alarm_rate_count,correct_rejection_rate_count),
    total_trial_num = c(length(hn_ith),length(hn_ith),length(no_hn_ith),length(no_hn_ith)))
  table$x <- factor(table$x, levels=c("CR","FA","H","M"))
  return(table)
}




formal_false_alarm_rate_each_subject <- function(dataset) {   
  hit_rate_count <- length(which(dataset[,4]==1 & dataset[,6]=="Correct"))
  miss_rate_count <- length(which(dataset[,4]==1 & dataset[,6]=="Incorrect"))
  False_alarm_rate_count <- length(which(dataset[,4]==0 & dataset[,6]=="Incorrect"))
  correct_rejection_rate_count <- length(which(dataset[,4]==0 & dataset[,6]=="Correct"))
  hit_rate  = hit_rate_count/length(which(dataset[,4]==1))
  miss_rate  = miss_rate_count/length(which(dataset[,4]==1))
  false_alarm_rate  = False_alarm_rate_count/length(which(dataset[,4]==0))
  correct_rejection_rate = correct_rejection_rate_count/length(which(dataset[,4]==0))
  table <- data.frame(
    y = c(hit_rate, miss_rate, false_alarm_rate, correct_rejection_rate),
    x = c("hit","miss","false_alarm","correct_rejection"),
    trial_num = c(hit_rate_count,miss_rate_count,False_alarm_rate_count,correct_rejection_rate_count),
    total_trial_num = c(length(which(dataset[,4]==1)),length(which(dataset[,4]==1)),length(which(dataset[,4]==0)),length(which(dataset[,4]==0))))
  table$x <- factor(table$x, levels=c("hit","miss","false_alarm","correct_rejection"))
  
  return(table)
}

plot_sm_trial_single_sub <- function(dataset) {   
  dataset_sc <- dataset[dataset[,2]=="sc",]
  dataset_dc <- dataset[dataset[,2]=="dc",]
  dataset_fpc <- dataset[dataset[,2]=="fpc",]
  dataset_tpc <- dataset[dataset[,2]=="tpc",]
  # dataset_hnd <- dataset[dataset[,2]=="hnd",]
  # length(which(dataset_hnd$counter==1))/length(dataset_hnd[,13])
  rate_table <- data.frame(
    x = c("dc","sc","cfp","ctp"),
    y = c(length(which(dataset_dc$counter==1))/length(dataset_dc[,13]),
          length(which(dataset_sc$counter==1))/length(dataset_sc[,13]),
          length(which(dataset_fpc$counter==1))/length(dataset_fpc[,13]),
          length(which(dataset_tpc$counter==1))/length(dataset_tpc[,13])),
    nt = c(length(dataset_dc[,13])-length(which(dataset_dc$counter==1)),
           length(dataset_sc[,13])-length(which(dataset_sc$counter==1)),
           length(dataset_fpc[,13])- length(which(dataset_fpc$counter==1)),
           length(dataset_tpc[,13])-length(which(dataset_tpc$counter==1)))
  )
  return(rate_table)
}

plot_sm_trial_single_sub_each_sess <- function(dataset, which_sess) {   
  dataset_cp = dataset
  dataset_cp = dataset[which(dataset_cp$SESS==which_sess),]
  dataset_sc <- dataset_cp[dataset_cp[,2]=="sc",]
  dataset_dc <- dataset_cp[dataset_cp[,2]=="dc",]
  dataset_fpc <- dataset_cp[dataset_cp[,2]=="fpc",]
  dataset_tpc <- dataset_cp[dataset_cp[,2]=="tpc",]
  # dataset_hnd <- dataset[dataset[,2]=="hnd",]
  # length(which(dataset_hnd$counter==1))/length(dataset_hnd[,13])
  rate_table <- data.frame(
    x = c("dc","sc","cfp","ctp"),
    y = c(length(which(dataset_dc$counter==1))/length(dataset_dc[,13]),
          length(which(dataset_sc$counter==1))/length(dataset_sc[,13]),
          length(which(dataset_fpc$counter==1))/length(dataset_fpc[,13]),
          length(which(dataset_tpc$counter==1))/length(dataset_tpc[,13])),
    nt = c(length(dataset_dc[,13])-length(which(dataset_dc$counter==1)),
           length(dataset_sc[,13])-length(which(dataset_sc$counter==1)),
           length(dataset_fpc[,13])- length(which(dataset_fpc$counter==1)),
           length(dataset_tpc[,13])-length(which(dataset_tpc$counter==1)))
  )
  return(rate_table)
}


plot_false_alarm_plot_single_sub <- function(rate_table) {   
  figure <- ggplot(data=rate_table, aes(x=rate_table$x, y=rate_table$y)) +  
    geom_bar(position=position_dodge(), stat="identity", width=.4) +
    theme_classic() +    
    labs(x="",y="Performance\n")+ 
    ggtitle(paste("","\n\n", sep = "")) +
    theme(axis.text.x = element_text(size=23), 
          axis.text.y =element_text(size=23), 
          axis.title.y = element_text(size=23),
          plot.title = element_text(size=0)) +  
    theme(plot.title = element_text(hjust = 0.5)) + 
    ylim(0,1) + 
    geom_text(aes(x=rate_table$x, y=0.05, , label=rate_table$trial_num, fontface=2, color="yellow"), size = 8) + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.4)) +
    theme(legend.position="none")
  figure
  return(figure)
}

plot_false_alarm_rate <- function(rate_table) {   
  table <- data.frame(
    individual = c(rate_table[,1],rate_table[,2],rate_table[,3],rate_table[,4]),
    mean = c(rep(mean(rate_table[,1],na.rm=TRUE),length(rate_table[,1])),
             rep(mean(rate_table[,2],na.rm=TRUE),length(rate_table[,1])),
             rep(mean(rate_table[,3],na.rm=TRUE),length(rate_table[,1])),
             rep(mean(rate_table[,4],na.rm=TRUE),length(rate_table[,1]))),
    se = c(rep(std(rate_table[,1]),length(rate_table[,1])),
           rep(std(rate_table[,2]),length(rate_table[,1])),
           rep(std(rate_table[,3]),length(rate_table[,1])),
           rep(std(rate_table[,4]),length(rate_table[,1]))),
    trial_number = c(rep(mean(num_table[,1]),length(rate_table[,1])),
                     rep(mean(num_table[,2]),length(rate_table[,1])),
                     rep(mean(num_table[,3]),length(rate_table[,1])),
                     rep(mean(num_table[,4]),length(rate_table[,1]))),
    x_label = c(rep("Hit",length(rate_table[,1])),
                rep("Miss",length(rate_table[,1])),
                rep("False alarm",length(rate_table[,1])),
                rep("Correct rej",length(rate_table[,1]))))
  figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
    geom_bar(position=position_dodge(), stat="identity", width=.5) +
    geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.1)) +
    geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
    theme_classic() +    
    labs(x="",y="Performance")+ 
    ggtitle(paste("","\n\n", sep = "")) +
    theme(axis.text.x = element_text(size=33),axis.text.y =element_text(size=37),axis.title.y = element_text(size=37),plot.title = element_text(size=37)) +  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + theme(axis.text.x = element_text(angle = 35, vjust = 0.5)) 
  
  return(figure)
}


plot_combined_false_alarm_rate <- function(rate_table_md, rate_table_fmri,font_size) {   
  table <- data.frame(
    individual = c(rate_table_md[,1],rate_table_md[,2],rate_table_md[,3],rate_table_md[,4],
                   rate_table_fmri[,1], rate_table_fmri[,2], rate_table_fmri[,3], rate_table_fmri[,4]),
    mean = c(rep(mean(rate_table_md[,1],na.rm=TRUE),12),
             rep(mean(rate_table_md[,2],na.rm=TRUE),12),
             rep(mean(rate_table_md[,3],na.rm=TRUE),12),
             rep(mean(rate_table_md[,4],na.rm=TRUE),12),
             rep(mean(rate_table_fmri[,1],na.rm=TRUE),10),
             rep(mean(rate_table_fmri[,2],na.rm=TRUE),10),
             rep(mean(rate_table_fmri[,3],na.rm=TRUE),10),
             rep(mean(rate_table_fmri[,4],na.rm=TRUE),10)),
    se = c(rep(std(rate_table_md[,1]),12),
           rep(std(rate_table_md[,2]),12),
           rep(std(rate_table_md[,3]),12),
           rep(std(rate_table_md[,4]),12),
           rep(std(rate_table_fmri[,1]),10),
           rep(std(rate_table_fmri[,2]),10),
           rep(std(rate_table_fmri[,3]),10),
           rep(std(rate_table_fmri[,4]),10)),
    x_label = c(rep("H",12),rep("M",12),rep("FA",12),rep("CR",12),
                rep("H",10),rep("M",10),rep("FA",10),rep("CR",10)),
    group = c(rep("Practice day", 12*4), rep("MEG day", 10*4)))
  
  
  table$group = factor(table$group, levels=c("Practice day", "MEG day"))
  
  figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean, fill=table$group)) +  
    
    geom_bar(position=position_dodge(), stat="identity", width=.7) +
    geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.7)) +
    geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1,position=position_dodge(0.7)) +
    theme_classic() +    
    scale_fill_manual("legend", values = c("Practice day" = "gray87","MEG day" = "gray40"))+
    labs(x="",y="Performance\n")+ 
    ggtitle(paste("","\n\n", sep = "")) +
    theme(axis.title.y = element_text(size=font_size))+
    theme(axis.title.x = element_text(size=font_size))+
    theme(legend.text=element_text(size=33))+
    theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + 
    theme(axis.text.x = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
    theme(axis.text.y = element_text(size=font_size, angle = 0, vjust = 0.5)) + 
    theme(legend.title=element_blank())+
    guides(fill = guide_legend(override.aes = list(shape = NA)))
  
  return(figure)
}







plot_correct_rate_map <- function(rate_table,map_table,num_table,title){   
  raw_data <- data.frame(
    rate = c(rate_table[,1],rate_table[,2],rate_table[,3]),
    map = factor(c(map_table[,1],map_table[,2],map_table[,3])),
    trial_num = c(num_table[,1],num_table[,2],num_table[,3]),
    sub = factor(rep(1:19,3))
  )
  
  library(nlme)
  m1 <- lme(rate~map,random=~1|sub,data=raw_data)
  anova(m1)
  
  
  library(lmerTest)
  lmm <- lmer(raw_data$rate~raw_data$map + (1|raw_data$sub))
  summary(lmm)
  anova(lmm)
  
  
  
  table <- data.frame(
    individual = raw_data$rate,
    mean = c(rep(mean(raw_data[which(raw_data$map==1),1]),length(which(raw_data$map==1))),
             rep(mean(raw_data[which(raw_data$map==2),1]),length(which(raw_data$map==2))),
             rep(mean(raw_data[which(raw_data$map==3),1]),length(which(raw_data$map==3))),
             rep(mean(raw_data[which(raw_data$map==4),1]),length(which(raw_data$map==4))),
             rep(mean(raw_data[which(raw_data$map==5),1]),length(which(raw_data$map==5))),
             rep(mean(raw_data[which(raw_data$map==6),1]),length(which(raw_data$map==6)))),
    se = c(rep(std(raw_data[which(raw_data$map==1),1]),length(which(raw_data$map==1))),
           rep(std(raw_data[which(raw_data$map==2),1]),length(which(raw_data$map==2))),
           rep(std(raw_data[which(raw_data$map==3),1]),length(which(raw_data$map==3))),
           rep(std(raw_data[which(raw_data$map==4),1]),length(which(raw_data$map==4))),
           rep(std(raw_data[which(raw_data$map==5),1]),length(which(raw_data$map==5))),
           rep(std(raw_data[which(raw_data$map==6),1]),length(which(raw_data$map==6)))),
    trial_number = c(rep(mean(raw_data[which(raw_data$map==1),3]),length(which(raw_data$map==1))),
                     rep(mean(raw_data[which(raw_data$map==2),3]),length(which(raw_data$map==2))),
                     rep(mean(raw_data[which(raw_data$map==3),3]),length(which(raw_data$map==3))),
                     rep(mean(raw_data[which(raw_data$map==4),3]),length(which(raw_data$map==4))),
                     rep(mean(raw_data[which(raw_data$map==5),3]),length(which(raw_data$map==5))),
                     rep(mean(raw_data[which(raw_data$map==6),3]),length(which(raw_data$map==6)))),
    x_label = c(rep("         Map 1",length(which(raw_data$map==1))),
                rep("         Map 2",length(which(raw_data$map==2))),
                rep("         Map 3",length(which(raw_data$map==3))),
                rep("         Map 4",length(which(raw_data$map==4))),
                rep("         Map 5",length(which(raw_data$map==5))),
                rep("         Map 6",length(which(raw_data$map==6)))))
  figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
    geom_bar(position=position_dodge(), stat="identity", width=.5) +
    geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.1)) +
    geom_hline(yintercept=0.33, linetype="dashed", color = "orange")+
    geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
    theme_classic() +    
    labs(x="",y="Performance")+ 
    ggtitle(paste(title,"\n\n", sep = "")) +
    theme(axis.text.x = element_text(size=20),axis.text.y =element_text(size=25),axis.title.y = element_text(size=25),plot.title = element_text(size=25)) +  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  return(figure)
}






plot_correct_rate_direction <- function(rate_table,num_table,num_of_direction, title) {   
  
  if(num_of_direction==4){
    table <- data.frame(
      individual = c(rate_table[,1],rate_table[,2],rate_table[,3],rate_table[,4]),
      mean = c(rep(mean(rate_table[,1]),length(rate_table[,1])),
               rep(mean(rate_table[,2]),length(rate_table[,1])),
               rep(mean(rate_table[,3]),length(rate_table[,1])),
               rep(mean(rate_table[,4]),length(rate_table[,1]))),
      se = c(rep(std(rate_table[,1]),length(rate_table[,1])),
             rep(std(rate_table[,2]),length(rate_table[,1])),
             rep(std(rate_table[,3]),length(rate_table[,1])),
             rep(std(rate_table[,4]),length(rate_table[,1]))),
      trial_number = c(rep(mean(num_table[,1]),length(rate_table[,1])),
                       rep(mean(num_table[,2]),length(rate_table[,1])),
                       rep(mean(num_table[,3]),length(rate_table[,1])),
                       rep(mean(num_table[,4]),length(rate_table[,1]))),
      x_label = c(rep("Direction 1",length(rate_table[,1])),
                  rep("Direction 2",length(rate_table[,1])),
                  rep("Direction 3",length(rate_table[,1])),
                  rep("Direction 4",length(rate_table[,1]))))
    figure <- ggplot(data=table, aes(x=table$x_label, y=table$mean)) +  
      geom_bar(position=position_dodge(), stat="identity", width=.35) +
      geom_errorbar(aes(ymin=table$mean-table$se, ymax=table$mean+table$se), colour="black", width=.2, position=position_dodge(0.1)) +
      geom_point(aes(x=table$x_label,y=table$individual),colour="red",size=6,shape=1) +
      theme_classic() +    
      geom_hline(yintercept=0.33, linetype="dashed", color = "orange")+
      labs(x="",y="Performance")+ 
      ggtitle(paste(title,"\n\n", sep = "")) +
      theme(axis.text.x = element_text(size=20),axis.text.y =element_text(size=25),axis.title.y = element_text(size=25),plot.title = element_text(size=25)) + theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) 
  }
  return(figure)
}
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}