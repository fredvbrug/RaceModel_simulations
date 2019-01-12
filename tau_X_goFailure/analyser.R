# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# script to process the SSRT estimates
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# clear everything
rm(list = ls()) 

# load libraries
library (ggplot2) # to plot the data
library (doBy) # for data summaries
library (reshape) # to melt and cast data
library (reshape2) # to melt and cast data
library (plyr)
library (Hmisc)

# -----------------------------------------------------------------------------
# STEP: load all the relevant output files
# -----------------------------------------------------------------------------
# create a new data frame to combine all the data
data <- data.frame() 

# determine how many simulated data files we have 
files <- dir(path='./processed_data', pattern='*.Rdata') 

# load the files and combine them
for (i in files) {
  name <- paste("./processed_data/", i, sep = "") # path to the file
  load(name) # load the file

  # add the content to the data frame
  data <- rbind (data, signal.cast)
  rm(signal.cast) 
}

# reorder levels pmissReq and tau
data$pmissReq = factor(data$pmissReq,levels(data$pmissReq)[c(1,5,2,3,4)])
data$tau = factor(data$tau,levels(data$tau)[c(1,5,2,3,4)])

# create numeric version of the subject factor; will be used in loops below
data$subjectNum <- as.numeric(data$subject)

# also create a version with the p(miss) = .20 data excluded
# (as the bias becomes too extreme for some estimation methods; esp. the respOnly and Tannock method)
data.subset <- subset(data, pmissReq %nin% c('20'))

# -----------------------------------------------------------------------------
# STEP: put the data in long format to compare the methods etc. 
# -----------------------------------------------------------------------------
SSRT.diff.long <- melt(data.subset,
                  id.vars=c("subject", "Ntrials", "pmissReq", "tau"), # ID variables - all the variables to keep but not split apart on
                  measure.vars=c("SSRTall_diff", "SSRTresp_diff", "SSRTadj_diff", "SSRTmean_diff" ) # The source columns
)

# rename a few things
colnames(SSRT.diff.long) <- c('subject', 'ntrials', 'pmiss', 'tau', 'e_method', 'SSRT_diff') 
levels(SSRT.diff.long$e_method) <- c('standard', 'respOnly', 'prespAdj', 'mean')

# -----------------------------------------------------------------------------
# STEP: check the estimation methods show consistent biases:
# we use the averages for this (with N = 1000) and plot the outcome as heat maps
# -----------------------------------------------------------------------------
SSRT.diff.cast <- cast (SSRT.diff.long, ntrials * tau * pmiss * e_method ~ ., mean) 
names(SSRT.diff.cast)[5] <- "SSRT_diff"

# for the heat maps, we will manually set the range (to ensure that 0 is in the middle)
tmp <- abs(SSRT.diff.cast$SSRT_diff)
li <- c(-max(tmp), max(tmp))

h.plot <- ggplot(data = SSRT.diff.cast, aes(x = tau, y = pmiss)) +
  geom_raster(aes(fill = SSRT_diff), interpolate = TRUE) +
  scale_fill_gradientn(colours = c("blue", "white", "red"),limits =li)+
  facet_grid(ntrials~e_method, labeller = label_both)

h.plot 

# create pdf version as well
pdf("./summary_data/heat_plot.pdf", paper="a4")
h.plot
dev.off()

# -----------------------------------------------------------------------------
# STEP: compare true SSRT and the estimated SSRT for each individual 
# and show the distribution of the difference scores (with violin plots)
# -----------------------------------------------------------------------------
v <- ggplot(SSRT.diff.long, aes(x=tau, y=SSRT_diff, fill=pmiss)) +
  geom_violin(draw_quantiles= c(0.25, 0.5, 0.75), trim = T) +
  facet_grid(.~e_method, labeller = label_both) +
  ylim(-400, 500)

# create the figures for all estimation methods  
v.plots <- dlply(SSRT.diff.long, .(ntrials), function(x) v %+% x +
                    facet_grid(.~ntrials+e_method, labeller = label_both))
v.plots 

# create pdf version as well
pdf("./summary_data/violin_plots.pdf", paper="a4")
v.plots
dev.off()

# -----------------------------------------------------------------------------
# STEP: do some checks at a group level: create groups of different sizes,
# and check the 'true vs. estimated SSRT' difference (at group level)
# -----------------------------------------------------------------------------

# function to create random groups of subjects 
NSubjects <- c(16, 32, 64, 128)
Nsamples <- 100 

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

# we want this for each combination of Ntrials, pmissReq and tau
group.data <- ddply(data.subset, .(Ntrials, tau, pmissReq), funcRandomGroup)

# create new factors
group.data$sample <- factor(group.data$sample) 
group.data$Nsubject <- factor(group.data$Nsubject) 

# put this in long format again
group.diff.long <- melt(group.data,
                       id.vars=c("sample", "Nsubject", "Ntrials", "pmissReq", "tau"), # ID variables - all the variables to keep but not split apart on
                       measure.vars=c("SSRTall_diff", "SSRTresp_diff", "SSRTadj_diff", "SSRTmean_diff" ) # The source columns
)

# rename a few things
colnames(group.diff.long) <- c('sample', 'nsubjects', 'ntrials', 'pmiss', 'tau', 'e_method', 'SSRT_diff') 
levels(group.diff.long$e_method) <- c('standard', 'respOnly', 'prespAdj', 'mean')

# create another violon plot; we will create different plots for each level of p(miss)
# we will do this by creating a basic figure template, and then use this as a 'function' in dlply
  p <- ggplot(group.diff.long, aes(x=tau, y=SSRT_diff, fill=pmiss)) +
    geom_violin(draw_quantiles= c(0.25, 0.5, 0.75), trim = T) +
    facet_grid(ntrials~e_method, labeller = label_both)

# create the figures for all estimation methods  
  grp.plts <- dlply(group.diff.long, .(nsubjects), function(x) p %+% x +
          facet_grid(nsubjects+ntrials~e_method, labeller = label_both))

  # create pdfs as well
  pdf("./summary_data/group_plots.pdf", paper="a4")
  grp.plts
  dev.off()

# also put the quantile data in a csv file
  q.cast <- cast (group.diff.long, nsubjects * ntrials * tau * pmiss ~ e_method, 
                  function(x) round(quantile(x,c(0.025, 0.5, 0.975))))
  write.csv(q.cast, "./summary_data/groups_quantiles.csv")

# as the resp-only and Tannock method produce consistent biases, 
# we will exclude them from the next figure
  group.diff.long.subset <- subset(group.diff.long, e_method %nin% c('respOnly', 'prespAdj'))
  
  grp.plts.subset <- dlply(group.diff.long.subset, .(nsubjects), function(x) p %+% x +
          facet_grid(nsubjects+ntrials~e_method, labeller = label_both))
  
  # create pdfs as well
  pdf("./summary_data/group_plots_subset.pdf", paper="a4")
  grp.plts.subset
  dev.off()

  