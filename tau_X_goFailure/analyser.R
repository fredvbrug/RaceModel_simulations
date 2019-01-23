# -----------------------------------------------------------------------------
# script to process the SSRT estimates
# -----------------------------------------------------------------------------

rm(list = ls()) # clear everything
library (ggplot2) # library to plot the data
library (reshape) # library required for melting data
library (plyr) # library for ddply
library (cowplot) # library to combine plots

# load the file with all processed data
load('./processed_data/combined_individual_data.Rdata') # individual data

# exclude subjects when signal-respond RT > no-signal RT
data$violation <- ifelse(data$raceCheck > 0, 0, 1) 
data.old <- data # keep a copy of the full data set though
data <- subset(data, violation == 0)

# put the individual data in long format as well
indiv.long <- melt(data,
                   id.vars=c("subject", "Ntrials", "pmissReq", "tau"), # ID variables - all the variables to keep but not split apart on
                   measure.vars=c("SSRTall_diff", "SSRTresp_diff", "SSRTadj_diff", "SSRTmean_diff" )) # The source columns
colnames(indiv.long)[5:6] <- c('e_method', 'SSRT_diff') # change a few names
levels(indiv.long$e_method) <- c('standard', 'respOnly', 'prespAdj', 'mean')

# ---- auxiliary function ----
# function to calculate reliability coefficient (correlation between true SSRT and estimated SSRT)
funcRC <- function(data){
  tmp.standard <- data.frame(e_method="standard", RC = cor(data$SSRTtrue,data$SSRTall))
  tmp.resp <- data.frame(e_method="respOnly", RC = cor(data$SSRTtrue,data$SSRTresp))
  tmp.adj <- data.frame(e_method="prespAdj", RC = cor(data$SSRTtrue,data$SSRTadj))
  tmp.mean <-  data.frame(e_method="mean", RC = cor(data$SSRTtrue,data$SSRTmean))
  
  RC <- rbind(tmp.standard,tmp.resp,tmp.adj,tmp.mean)
  return(RC)
}

# not really a function, but we will use this for various plots
recurring.grid <- facet_grid(
  e_method~Ntrials,
  labeller = labeller(
    Ntrials = c(`100` = "Total N: 100 (25 signals)", 
               `200` = "Total N: 200 (50 signals)",
               `400` = "Total N: 400 (100 signals)",
               `800` = "Total N: 800 (200 signals)"),
    e_method = c(`standard` = "Integration (all)", 
               `mean` = "Mean"))
  )

# ---- get some more info about the subjects for whom signal-respond > no-signal RT ----
  # check for each parameter combination  
  violation.data <- ddply(data.old, .(Ntrials, tau, pmissReq), summarise, number=sum(violation)) 
  violation.data$e_method <- '' # only need this to keep layout the same as other plots
  
  # calculate correlation collapsed across p(miss and tau)
  collapsed.vd <- ddply(data.old, .(Ntrials), summarise, number=sum(violation))
  plotC.text <- collapsed.vd
  plotC.text$label <- paste("Total excl. = ", plotC.text$number) 
  
  # plot the data
  violation.plot <- ggplot(data = violation.data, aes(x = tau, y = pmissReq)) +
    geom_raster(aes(fill = number), interpolate = TRUE) +
    scale_fill_gradientn(colours = c("gold", "black"))+
    facet_grid(e_method~Ntrials, labeller = labeller(
      Ntrials = c(`100` = "Total N: 100 (25 signals)", 
                  `200` = "Total N: 200 (50 signals)",
                  `400` = "Total N: 400 (100 signals)",
                  `800` = "Total N: 800 (200 signals)")))+
    geom_text(data=plotC.text, aes(x=1, y=1, label = label), hjust = 0,
              colour="black", inherit.aes=FALSE, parse=FALSE)+
    theme(strip.text=element_text(margin = margin(0.1,0.1,0.1,0.1, "cm")),
          strip.background.y =element_rect(fill="white"))+
    xlab('Tau of the no-signal RT distribution')+
    ylab('Go failure (%)')+
    labs(fill = "Number of\nexcluded subjects")
  violation.plot
  
  # get an idea of how all this influences SSRT  
  # violations occur primarily when tau is low, p(miss) is high, and number of trials is low
  ddply(subset(data.old, Ntrials == 100 & tau %in% c('1', '50') & pmissReq %in% c(10, 15, 20)), 
        .(Ntrials, violation), summarise, 
        SSRTall_diff=mean(SSRTall_diff)) # bias for included vs. excluded subjects (standard integration method)
  
  ddply(subset(data.old, Ntrials == 100 & tau %in% c('1', '50') & pmissReq %in% c(10, 15, 20)), 
        .(Ntrials, violation), summarise, 
        SSRTmean_diff=mean(SSRTmean_diff)) # bias for included vs. excluded subjects (mean method) 


# ---- calculate the mean difference between true and estimated SSRT   ----
  # (to determine if there are consistent biases)
  indiv.mean <- ddply(indiv.long, .(Ntrials, tau, pmissReq, e_method), summarise, SSRT_diff=mean(SSRT_diff)) # calculate mean for all possible combinations
  indiv.mean <- subset(indiv.mean, e_method %in% c('standard', 'mean')) # only show the standard integration approach and mean approach
  plot.text <- ddply(indiv.mean, .(e_method,Ntrials), summarise, SSRT_diff.sd =sd(SSRT_diff)) # calculate standard deviations (to show as text in the plot)
  plot.text$label <- paste("SD: ", round(plot.text$SSRT_diff.sd), "ms") # adjust the text label
  plot.li <- c(-max(abs(indiv.mean$SSRT_diff)), max(abs(indiv.mean$SSRT_diff))) # ensure that 0 is in the middle of the range
  
  # plot the data
  indiv.mean.plot <- ggplot(data = indiv.mean, aes(x = tau, y = pmissReq)) +
    geom_raster(aes(fill = SSRT_diff), interpolate = TRUE) +
    scale_fill_gradientn(colours = c("blue", "white", "magenta"),limits =plot.li)+
    recurring.grid+
    geom_text(data=plot.text, aes(x=1, y=1, label = label), hjust = 0,
              colour="black", inherit.aes=FALSE, parse=FALSE)+
    xlab('Tau of the no-signal RT distribution')+
    ylab('Go failure (%)')+
    labs(fill = "Estimated -\n true SSRT")+
    theme(strip.text=element_text(margin = margin(0.1,0.1,0.1,0.1, "cm")))
  indiv.mean.plot
    
# ---- calculate the correlation between true and estimated SSRT ----
  # (to check reliability of the estimates)
  indiv.rc <- ddply(data, .(Ntrials, pmissReq, tau), funcRC)  # calculate correlation (via funcRC) for all possible combinations
  indiv.rc <- subset(indiv.rc, e_method %in% c('standard', 'mean')) # only show the standard integration approach and mean approach

  # calculate correlation collapsed across p(miss and tau)
  collapsed.rc <- ddply(data, .(Ntrials), funcRC) # use funcRC (see above) for the calculations
  collapsed.rc <- subset(collapsed.rc, e_method %in% c('standard', 'mean')) # only show the standard integration approach and mean approach
  plotB.text <- arrange(collapsed.rc, e_method) # same order as indiv.rc
  plotB.text$label <- paste(" Overall R = ", format(round(plotB.text$RC, 3), nsmall = 3)) # show the overall correlations as text

  # plot the data
  indiv.rc.plot <- ggplot(data = indiv.rc, aes(x = tau, y = pmissReq)) +
    geom_raster(aes(fill = RC), interpolate = TRUE) +
    scale_fill_gradientn(colours = c("red", "yellow", "green"),limits =c(0,1))+
    recurring.grid +
    geom_text(data=plotB.text, aes(x=1, y=1, label = label), hjust = 0,
              colour="black", inherit.aes=FALSE, parse=FALSE) +
    xlab('Tau of the no-signal RT distribution')+
    ylab('Go failure (%)')+
    labs(fill = "Estimated -\n true SSRT")+
    theme(strip.text=element_text(margin = margin(0.1,0.1,0.1,0.1, "cm")))
  indiv.rc.plot
  
# ---- combine some plots for the paper ----
  violation.plot <- violation.plot+theme(axis.title.x=element_text(color='white'))
  indiv.mean.plot <- indiv.mean.plot+theme(axis.title.x=element_text(color='white'))
  plot_grid(violation.plot, indiv.mean.plot, indiv.rc.plot, 
            labels=c("A", "B", "C"), ncol = 1, nrow = 3, align='v', axis='l', rel_heights=c(0.6,1,1))
  
# ---- plot the distribution of the true-vs-estimated SSRT differences ----
# create violin plots
  violin.plot <- ggplot(indiv.long, aes(x=e_method, y=SSRT_diff, fill=e_method)) +
    geom_hline(yintercept=0)+
    geom_violin(trim = T) +
    scale_fill_manual(labels = c("Integration\n(all)", "Integration\n(respond only)", "Integration\n(p adjusted)", "Mean"), 
                      values = c('red', 'green', 'blue', 'yellow'))+
    facet_grid(
      Ntrials~pmissReq*tau,
      labeller = labeller(
        Ntrials = c(`100` = "Total N: 100 (25 signals)", 
                    `200` = "Total N: 200 (50 signals)",
                    `400` = "Total N: 400 (100 signals)",
                    `800` = "Total N: 800 (200 signals)"),
        pmissReq = c(`0` = "p(m):0",
                     `5` = "p(m):5",
                     `10` = "p(m):10",
                     `15` = "p(m):15",
                     `20` = "p(m):20"),
        tau = c(`1` = "t:1",
                `50` = "t:50",
                `100` = "t:100",
                `150` = "t:150",
                `200` = "t:200")
      )
    ) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text=element_text(size=10,margin = margin(0.1,0.1,0.1,0.1, "cm")),
          legend.key.size = unit(2, 'lines'),
          legend.text=element_text(size=10)
          )+
    labs(fill = "Estimation\nmethod")+
    ylab('Estimated - true SSRT') 
  
  violin.plot   