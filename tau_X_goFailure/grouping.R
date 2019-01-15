# clear everything
rm(list = ls()) 

# load libraries
library (ggplot2) # to plot the data
library (doBy) # for data summaries
library (reshape) # to melt and cast data
library (reshape2) # to melt and cast data
library (plyr)
library (Hmisc)

# load the combined file with all processed data
load('./processed_data/combined_individual_data.Rdata')

# some predefined variables
NSubjects <- c(16, 32, 64, 128)
Nsamples <- 100 

# function to create random groups of subjects 
funcRandomGroup <- function(data){
  data.group <- data.frame()
  
  for (s in 1:Nsamples){
    for (n in NSubjects){
      # get a random group of subjects
      tmp.group <- subset(data, subjectNum %in% sample(1:1000, n)) 
      
      # calculate the mean differences for these groups
      tmp.SSRTall_diff <- mean(tmp.group$SSRTall_diff)
      tmp.SSRTresp_diff <- mean(tmp.group$SSRTresp_diff)
      tmp.SSRTadj_diff <- mean(tmp.group$SSRTadj_diff)
      tmp.SSRTmean_diff <- mean(tmp.group$SSRTmean_diff)
      
      # combine the means of the group
      tmp.combined <- data.frame(sample = s, 
                                 Nsubject = n, 
                                 SSRTall_diff = tmp.SSRTall_diff, 
                                 SSRTresp_diff = tmp.SSRTresp_diff, 
                                 SSRTadj_diff = tmp.SSRTadj_diff, 
                                 SSRTmean_diff = tmp.SSRTmean_diff)
      
      # combine the group means with the rest of the data
      data.group <- rbind(data.group, tmp.combined)
      
      # print progress
      text <- sprintf ("Ntrials: L %d, tau: L %d, pmiss: L %d; sample %d; Nsubjects %d", 
                       data$Ntrials[1], data$tau[1], data$pmissReq[1], s, n)
      print(text)  
    }
  }
  
  # return data.group
  return (data.group)
}

# we want groups for each combination of Ntrials, pmissReq and tau
group.data <- ddply(data, .(Ntrials, tau, pmissReq), funcRandomGroup)

# create new factors
group.data$sample <- factor(group.data$sample) 
group.data$Nsubject <- factor(group.data$Nsubject) 

# save the grouped data
save (group.data, file = './processed_data/combined_group_data.Rdata')
