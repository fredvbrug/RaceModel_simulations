# -----------------------------------------------------------------------------
# script to determine achieved power (based on simulations)
# -----------------------------------------------------------------------------

rm(list = ls()) # clear everything
library (ggplot2) # library to plot the data
library (reshape) # to melt and cast data
library (plyr) # library required for ddply 
library (cowplot) # library to combine plots
library (Hmisc) # library required for %nin% 
library (Cairo) # to save the plots 

# ---- auxiliary functions ----

# function to create random groups of subjects 
funcRandomGroup <- function(data, set, nsub){
  p.values <- data.frame()
  data$subjectNum <- as.numeric(data$subject)
  
  # split in different groups
  group1 <- subset(data, group == '1')
  group2 <- subset(data, group == '2')
  
  # exclude subjects if signal-respond > no-signal RT
  # this script can also be used for repeated measures, so all subjects with 
  # the same subject number have to be excluded from all groups
  violations <- subset(data, raceCheck < 0)
  included <- 1:nrow(group1) 
  included <- included[included %nin% violations$subject]
    
  # get random groups of subjects
  for (n in nsub){
    for (s in 1:Nsamples){
      sampled.subjects <- sample(included, n)
      group1.subset <- subset(group1, subjectNum %in% sampled.subjects) 
      group2.subset <- subset(group2, subjectNum %in% sampled.subjects) 
      
      # within-subjects (paired) or between-subjects comparison of estimated SSRTs using a simple t-test
      if (paired == TRUE){
        test <- t.test(group1.subset$SSRTall, group2.subset$SSRTall, paired=TRUE)
        }else{
        test <- t.test(group1.subset$SSRTall, group2.subset$SSRTall, paired=FALSE, var.equal = FALSE)
      }
      
      # extract the p-value 
      tmp <- data.frame(NSubjects = n, sample = s, p=round(test$p.value[[1]],4))
      
      # combine the data
      p.values <- rbind(p.values, tmp)
      
      # print progress
      text <- sprintf ("Sim %s, Nsubjects %d; Sample %d", set, n, s)
      print(text)  
    }
  }
      
  # return data.group
  return (p.values)
}

# determine for each sample size how often p < .05
funcSign <- function(data){
  ach.power <- nrow(subset(data, p < .05))/nrow(data)
  return(ach.power)
}

# ---- for each set of simulations, estimate achieved power ----

# some predefined variables
N1 <- c(1:6)*32 # group size can differ for medium effect
N2 <- c(1:6)*16 # group size can differ for large effect
NSubjects <- rbind(N1, N2)
paired <- FALSE # independent groups
Nsamples <- 1000 # how many samples (per group size and parameter combination) do we want? 

# get file names
files <- dir(path='./simulated_data/', pattern='simulation_')

# process the simulations
ach.power.combined <- data.frame()
for (i in files) {
  name <- paste("./simulated_data/", i, sep = "")
  load (name)

  # determine if there is a medium or very large SSRT difference
  tmp <- unlist(strsplit(i, "\\."))[1]
  tmp <- strsplit(tmp, "_")
  sim.number <- as.numeric(unname(sapply(tmp, `[[`, 3)))
  ssrt.diff <- ifelse(sim.number %% 2 == 1, 1, 2)
  
  # determine p-values of differences between samples
  p.values <- funcRandomGroup(simulation, i, NSubjects[ssrt.diff,])
  p.values$NSubjects <- factor(p.values$NSubjects)
  
  # determine how often the differences are significant 
  ach.power <- ddply(p.values, .(NSubjects), funcSign)  
  ach.power$simulation <- sim.number
  ach.power$ssrt.diff <- ifelse(sim.number %% 2 == 1, 15, 25)
  ach.power.combined <- rbind(ach.power.combined, ach.power)
}

# reorder simulations and rename
ach.power.combined <- arrange(ach.power.combined, simulation)
names(ach.power.combined)[2] <- 'achieved_power'

# save this for later reuse (so that we don't have to rerun everything when we want to update the figure)
save(ach.power.combined, file = './processed_data/power.tests.Rdata')
write.csv(ach.power.combined, './processed_data/power.tests.csv')

# ---- plot achieved power ----

# load file
rm(list = ls()) 
load('./processed_data/power.tests.Rdata')     


# get relevant information about the simulations (stored in a different file)
load('./simulated_data/information_simulations.Rdata')
sim.info$GORT.difference <- (sim.info$G2_Gmu + sim.info$G2_Gtau) - (sim.info$G1_Gmu + sim.info$G1_Gtau)
sim.info$SSRT.difference <- (sim.info$G2_Smu + sim.info$G2_Stau) - (sim.info$G1_Smu + sim.info$G1_Stau)
sim.info$miss.difference <- (sim.info$G2_Gmiss) - (sim.info$G1_Gmiss)

# now put this extra info in the achieved power data frame
ach.power.combined$ntrials <-  sim.info[ach.power.combined$simulation, 'ntrials']

# labels for figure (will put true SSRT difference in y-axis label)
temp <- c('a. ','b. ', 'c. ', 'd. ', 'e.', 'f. ', 'g. ', 'h. ', 'i. ')

ach.power.combined$code <- 1+floor((ach.power.combined$simulation-1)/2)
  
ach.power.combined$label <-  paste( temp[1+floor((ach.power.combined$simulation-1)/2)],
                                    'Total N = ', sim.info[ach.power.combined$simulation, 'ntrials'],
                                    '\n(stop signals = ', sim.info[ach.power.combined$simulation, 'ntrials']/4, ')',
                                   '\nGoRT ∆ = ', sim.info[ach.power.combined$simulation, 'GORT.difference'], ' ms',
                                   '\nP(miss) ∆ = .', sim.info[ach.power.combined$simulation, 'miss.difference'],
                                   sep="")

# relevel NSubject 
ach.power.combined$NSubjects <- ordered(ach.power.combined$NSubjects, levels = c("16", "32", "48", "64", "80", "96", 
                                                                                 "128", "160", "192"))


# show in a graph (but use different subpanels for small and large SSRT differences)
SSRT15 <- subset(ach.power.combined, simulation %% 2 == 1)
SSRT25 <- subset(ach.power.combined, simulation %% 2 == 0)

A <-  ggplot(data=SSRT15, aes(x=NSubjects, y=achieved_power)) +
    geom_bar(stat="identity")+
    geom_hline(yintercept = .80, colour = 'red')+
    facet_wrap(~label, nrow=3, ncol=3) +
    ylab('Achieved power when true SSRT ∆ = 15 ms (Cohen`s d ≈ .50-55)')+
    xlab('Number of subjects per group')+
    theme_bw()
    
B <- A %+% SSRT25 + ylab('Achieved power when true SSRT ∆ = 25 ms (Cohen`s d ≈ .85-90)')

# combine two plots  
AB <- plot_grid(A, B,labels=c("A", "B"), ncol = 1, nrow = 2, align = 'hv')
AB

ggsave(filename="./summary_data/power.plot.eps", plot=AB, device=cairo_ps,
       height = 297, width = 210, units = "mm", dpi = 1200)
  

