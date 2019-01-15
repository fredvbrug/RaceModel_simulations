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


# load the combined file with all processed data
load('./processed_data/combined_individual_data.Rdata')

# also create a version with the p(miss) = .20 data excluded
# (as the bias becomes too extreme for some estimation methods; esp. the respOnly and Tannock method)
data.subset <- subset(data, pmissReq %nin% c('20'))

# -----------------------------------------------------------------------------
# STEP: put the individual data in long format to compare the methods etc. 
# -----------------------------------------------------------------------------
indiv.long <- melt(data.subset,
                  id.vars=c("subject", "Ntrials", "pmissReq", "tau"), # ID variables - all the variables to keep but not split apart on
                  measure.vars=c("SSRTall_diff", "SSRTresp_diff", "SSRTadj_diff", "SSRTmean_diff" ) # The source columns
)

# rename a few things
colnames(indiv.long) <- c('subject', 'ntrials', 'pmiss', 'tau', 'e_method', 'SSRT_diff') 
levels(indiv.long$e_method) <- c('standard', 'respOnly', 'prespAdj', 'mean')

# -----------------------------------------------------------------------------
# STEP: check the estimation methods show consistent biases:
# we use the averages for this (with N = 1000) and plot the outcome as heat maps
# -----------------------------------------------------------------------------
tic()
SSRT.diff.cast <- cast (indiv.long, ntrials * tau * pmiss * e_method ~ ., mean, value='SSRT_diff') 
names(SSRT.diff.cast)[5] <- "SSRT_diff"
toc()

tic()
test <- ddply(indiv.long, .(ntrials, tau, pmiss, e_method), summarise, SSRT_diff=mean(SSRT_diff))
toc()



text.tmp <- summaryBy(SSRT_diff ~ e_method*ntrials, data=as.data.frame(SSRT.diff.cast),  FUN=sd)
text.tmp$label <- paste("SD = ", round(text.tmp$SSRT_diff.sd), "ms")
  
# for the heat maps, we will manually set the range (to ensure that 0 is in the middle)
tmp <- abs(SSRT.diff.cast$SSRT_diff)
li <- c(-max(tmp), max(tmp))

h.plot <- ggplot(data = SSRT.diff.cast, aes(x = tau, y = pmiss)) +
  geom_raster(aes(fill = SSRT_diff), interpolate = TRUE) +
  scale_fill_gradientn(colours = c("blue", "white", "red"),limits =li)+
  facet_grid(ntrials~e_method, labeller = label_both)+
  geom_text(data=text.tmp, aes(x=2, y=1, label = label), 
            colour="black", inherit.aes=FALSE, parse=FALSE)

h.plot 

# create pdf version as well
pdf("./summary_data/heat_plot.pdf", paper="a4")
h.plot
dev.off()

# -------------------------------------------------------------------------------
# STEP: show for the distance between 0.025 and 0.975 quantiles (using heat maps)
# -------------------------------------------------------------------------------
# function calculate the 0.95 'quantile' interval (based on the 0.025 and 0.975 quantile)
funcQI <- function(data){
  tmp <- quantile(data$SSRT_diff, probs = c(0.025, 0.975))
  QI <- tmp[2]-tmp[1]
  names(QI) <- 'interval'
  return(QI)
}

# calculate the 0.95 'quantile' interval for each combo
QI.data <- ddply(SSRT.diff.long, .(ntrials, tau, pmiss, e_method), funcQI)

# use this for the heat maps 
ggplot(data = QI.data, aes(x = tau, y = pmiss)) +
  geom_raster(aes(fill = interval), interpolate = T) +
  scale_fill_gradientn(colours = c("white", "black"))+
  facet_grid(ntrials~e_method, labeller = label_both)+
  labs(fill = "0.95 QI")

# -----------------------------------------------------------------------------
# STEP: determine the reliability of each method by correlating true SSRT with 
# estimated SSRT
# -----------------------------------------------------------------------------
# function to calculate reliability coefficient 
funcRC <- function(data){
  tmp.standard <- cor(data$SSRTtrue,data$SSRTall)
  tmp.resp <- cor(data$SSRTtrue,data$SSRTresp)
  tmp.adj <- cor(data$SSRTtrue,data$SSRTadj)
  tmp.mean <-  cor(data$SSRTtrue,data$SSRTmean)
  
  RC <- data.frame(standard = tmp.standard, respOnly = tmp.resp, prespAdj = tmp.adj, mean = tmp.mean)
  return(RC)
}

# calculate the correlations for each combo
RC.data <- ddply(data.subset, .(Ntrials, pmissReq, tau), funcRC) 
RC.data.long <- melt(RC.data,
                       id.vars=c("Ntrials", "pmissReq", "tau"), # ID variables - all the variables to keep but not split apart on
                       measure.vars=c("standard", "respOnly", "prespAdj", "mean")) # The source columns
colnames(RC.data.long) <- c('ntrials', 'pmiss', 'tau', 'e_method', 'RC') 


ggplot(data = RC.data.long, aes(x = tau, y = pmiss)) +
  geom_raster(aes(fill = RC), interpolate = TRUE) +
  scale_fill_gradientn(colours = c("red", "yellow", "green"),limits =c(0,1))+
    facet_grid(ntrials~e_method, labeller = label_both)


# -----------------------------------------------------------------------------
# STEP: do some checks at a group level: create groups of different sizes,
# and check the 'true vs. estimated SSRT' difference (at group level)
# -----------------------------------------------------------------------------

# load the group data
load('./processed_data/combined_group_data.Rdata')

# put this in long format again
group.diff.long <- melt(group.data,
                       id.vars=c("sample", "Nsubject", "Ntrials", "pmissReq", "tau"), # ID variables - all the variables to keep but not split apart on
                       measure.vars=c("SSRTall_diff", "SSRTresp_diff", "SSRTadj_diff", "SSRTmean_diff" ) # The source columns
)

# rename a few things
colnames(group.diff.long) <- c('sample', 'nsubjects', 'ntrials', 'pmiss', 'tau', 'e_method', 'SSRT_diff') 
levels(group.diff.long$e_method) <- c('standard', 'respOnly', 'prespAdj', 'mean')

# summarise over samples
group.diff.cast <- cast (group.diff.long, nsubjects * ntrials * tau * pmiss * e_method ~ ., mean) 
names(group.diff.cast)[6] <- "SSRT_diff"

test <- subset(group.diff.cast, e_method %in% c('standard','mean') & pmiss %nin% c('20'))

tmp <- abs(test$SSRT_diff)
li <- c(-max(tmp), max(tmp))


# heat maps again
p <- ggplot(data = test, aes(x = tau, y = pmiss)) +
  geom_raster(aes(fill = SSRT_diff), interpolate = TRUE) +
  scale_fill_gradientn(colours = c("blue", "white", "red"),limits=li)+
  facet_grid(ntrials~e_method, labeller = label_both)

# create the figures for all estimation methods  
grp.plts <- dlply(test, .(nsubjects), function(x) p %+% x +
                    facet_grid(nsubjects+ntrials~e_method, labeller = label_both))
grp.plts


ggplot(data = test, aes(x = tau, y = pmiss)) +
  geom_raster(aes(fill = SSRT_diff), interpolate = TRUE) +
  scale_fill_gradientn(colours = c("blue", "white", "red"),limits=li)+
  facet_grid(ntrials~nsubjects+e_method, labeller = label_both)


# -----------------------------------------------------------------------------
# EXTRA: VIOLIN PLOTS TO SHOW THE DISTRIBUTIONS IN A DIFFERENT WAY  
# -----------------------------------------------------------------------------
  
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
  
  