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

# set theme for the first three plots
theme_size = 11
common_theme <- function(){
  theme(axis.text.x=element_text(angle=90, hjust=1, size = theme_size-1),
        axis.text.y=element_text(size = theme_size-1),
        strip.text=element_text(margin = margin(0.1,0.1,0.1,0.1, "cm")),
        legend.title=element_blank(),
        legend.position= "top", 
        # legend.box.background = element_rect(fill = "white", colour = "black"),
        legend.key.width=unit(2.2,"line"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(hjust=0, size = theme_size),
        text = element_text(size = theme_size)
        )
}


# not really a function, but we will use a similar grid for various plots
recurring.grid <- facet_grid(
  Ntrials~e_method,
  labeller = labeller(
    Ntrials = c(`100` = "Total N: 100\n(25 stop signals)", 
               `200` = "Total N: 200\n(50 stop signals)",
               `400` = "Total N: 400\n(100 stop signals)",
               `800` = "Total N: 800\n(200 stop signals)"),
    e_method = c(`standard` = "Integration\n(w. replacement)", 
               `mean` = "Mean"))
  )

# ---- get some more info about the subjects for whom signal-respond > no-signal RT ----
  # check for each parameter combination  
  violation.data <- ddply(data.old, .(Ntrials, tau, pmissReq), summarise, number=mean(violation)) 
  violation.data$number <- violation.data$number * 100
  violation.data$e_method <- 'standard' # only need this to keep layout the same as other plots
  
  # plot the data
  violation.plot <- ggplot(data = violation.data, aes(x = tau, y = pmissReq)) +
    geom_raster(aes(fill = number), interpolate = TRUE) +
    scale_fill_gradientn(colours = c("#E69F00", "#000000"))+
    recurring.grid +
    xlab('Tau of the\ngo RT distribution')+
    ylab('Go omission (%)')+
    labs(title = "Percentage of\nexcl. subjects") +
    common_theme()+
    theme(strip.background.x =element_rect(fill="white"),
          strip.text.x = element_text(color = "white"),
          legend.key.width=unit(1,"line"))
  violation.plot
  
  # get an idea of how all this influences SSRT  
  # violations occur primarily when tau is low, p(miss) is high, and number of trials is low
  ddply(subset(data.old, Ntrials == 100 & tau %in% c('1', '50') & pmissReq %in% c(10, 15, 20)), 
        .(Ntrials, violation), summarise, 
        SSRTall_diff=mean(SSRTall_diff)) # bias for included vs. excluded subjects (integration method w. replacement)
  
  ddply(subset(data.old, Ntrials == 100 & tau %in% c('1', '50') & pmissReq %in% c(10, 15, 20)), 
        .(Ntrials, violation), summarise, 
        SSRTall_diff=mean(SSRTresp_diff)) # bias for included vs. excluded subjects (integration method without replacement)
  
  ddply(subset(data.old, Ntrials == 100 & tau %in% c('1', '50') & pmissReq %in% c(10, 15, 20)), 
        .(Ntrials, violation), summarise, 
        SSRTall_diff=mean(SSRTadj_diff)) # bias for included vs. excluded subjects (integration method with adjusted p(respond|signal))
  
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
    scale_fill_gradientn(colours = c("#0072B2", "white", "#CC79A7"),limits =plot.li)+
    recurring.grid +
    geom_text(data=plot.text, aes(x=1, y=1, label = label), hjust = 0,
              colour="black", inherit.aes=FALSE, parse=FALSE, size = 3)+
    xlab('Tau of the\ngo RT distribution')+
    ylab('Go omission (%)')+
    labs(title = "Difference (in ms)\nestimated - true SSRT")+
    common_theme()
    
  indiv.mean.plot
    
# ---- calculate the correlation between true and estimated SSRT ----
  # (to check reliability of the estimates)
  indiv.rc <- ddply(data, .(Ntrials, pmissReq, tau), funcRC)  # calculate correlation (via funcRC) for all possible combinations
  indiv.rc.all <- indiv.rc # before excluding the two other methods, keep a copy of the original data frame
  indiv.rc <- subset(indiv.rc, e_method %in% c('standard', 'mean')) # only show the standard integration approach and mean approach

  # calculate correlation collapsed across p(miss and tau)
  collapsed.rc <- ddply(data, .(Ntrials), funcRC) # use funcRC (see above) for the calculations
  collapsed.rc <- subset(collapsed.rc, e_method %in% c('standard', 'mean')) # only show the standard integration approach and mean approach
  plotB.text <- arrange(collapsed.rc, e_method) # same order as indiv.rc
  plotB.text$label <- paste("Overall\nR: ", format(round(plotB.text$RC, 3), nsmall = 3)) # show the overall correlations as text

  # plot the data
  indiv.rc.plot <- ggplot(data = indiv.rc, aes(x = tau, y = pmissReq)) +
    geom_raster(aes(fill = RC), interpolate = TRUE) +
    scale_fill_gradientn(colours = c("#D55E00", "white", "#009E73"),limits =c(0,1))+
    recurring.grid +
    geom_text(data=plotB.text, aes(x=1, y=1.5, label = label), hjust = 0,
              colour="black", inherit.aes=FALSE, parse=FALSE, size = 3) +
    xlab('Tau of the\ngo RT distribution')+
    ylab('Go omission (%)')+
    labs(title = "Correlation\nestimated - true SSRT")+
    common_theme()
  indiv.rc.plot
  
  
# ---- combine some plots for the paper ----
  # violation.plot <- violation.plot+theme(strip.text.y = element_blank())
  # indiv.mean.plot <- indiv.mean.plot+theme(strip.text.y = element_blank(),
  #                                          axis.title.y=element_blank(),
  #                                          axis.text.y=element_blank())
  # indiv.rc.plot <- indiv.rc.plot + theme(axis.title.y=element_blank(),
  #                                        axis.text.y=element_blank())
    
  combined.plot <- plot_grid(violation.plot, indiv.mean.plot, indiv.rc.plot, 
            labels=c("A", "B", "C"), ncol = 3, nrow = 1, align='v', axis='l', rel_widths=c(0.65,1,1)) 
  combined.plot
  
    ggsave(filename="./summary_data/main_summary.eps", plot=combined.plot,
         height = 175, width = 200, units = "mm", dpi = 1200)
  
  
# ---- plot the distribution of the true-vs-estimated SSRT differences for the appendix ----
# create violin plots
    ntrials1 <- subset(indiv.long, Ntrials == "100")
    ntrials2 <- subset(indiv.long, Ntrials == "200")
    ntrials3 <- subset(indiv.long, Ntrials == "400")
    ntrials4 <- subset(indiv.long, Ntrials == "800")
    
    violin.ntrials1 <- ggplot(ntrials1, aes(x=pmissReq, y = SSRT_diff, fill=e_method)) +
      geom_violin(trim = TRUE) +
      geom_hline(yintercept=0, colour="grey30")+
      scale_fill_manual(labels = c("Integration\nomissions\nreplaced", 
                                   "Integration\nomissions\nexcluded", 
                                   "Integration\np(respond|signal)\nadjusted", 
                                   "Mean"), 
                        values = c('#E69F00', '#56B4E9', '#009E73', '#CC79A7'))+
      facet_grid(.~tau,
        labeller = labeller(
          tau = c(`1` = "tau go = 1",
                  `50` = "tau go = 50",
                  `100` = "tau go = 100",
                  `150` = "tau go = 150",
                  `200` = "tau go = 200")
        )
      ) +
      theme(strip.text=element_text(margin = margin(0.1,0.1,0.1,0.1, "cm")),
            axis.text.x=element_text(angle=90, hjust=1, size = theme_size-1),
            axis.text.y=element_text(size = theme_size-1),
            legend.key.width = unit(1, 'lines'),
            legend.key.height=unit(3,"line"),
            legend.spacing.x = unit(0.5, 'cm'),
            legend.position="bottom",
            legend.title = element_blank(),
            plot.background = element_rect(fill = "white"),
            plot.title = element_text(hjust=0, size = theme_size),
            text = element_text(size = theme_size)
            )+
      labs(fill = "Estimation\nmethod") +
      xlab('Go omission (%)') +
      ylab('Difference estimated - true SSRT (in ms)') +
      coord_flip() +
      labs(title = "A. Total N: 100 (25 stop signals)")

    violin.ntrials2 <- violin.ntrials1 %+% ntrials2 + labs(title = "B. Total N: 200 (50 stop signals)")
    violin.ntrials3 <- violin.ntrials1 %+% ntrials3 + labs(title = "C. Total N: 400 (100 stop signals)")
    violin.ntrials4 <- violin.ntrials1 %+% ntrials4 + labs(title = "D. Total N: 800 (200 stop signals)")
    
    violin.ntrials1   
    violin.ntrials2
    violin.ntrials3
    violin.ntrials4
    
  ggsave(filename="./summary_data/violin.ntrials1.eps", plot=violin.ntrials1,
         height = 16, width = 16, units = "cm", dpi = 1200)
  ggsave(filename="./summary_data/violin.ntrials2.eps", plot=violin.ntrials2,
         height = 16, width = 16, units = "cm", dpi = 1200)
  ggsave(filename="./summary_data/violin.ntrials3.eps", plot=violin.ntrials3,
         height = 16, width = 16, units = "cm", dpi = 1200)
  ggsave(filename="./summary_data/violin.ntrials4.eps", plot=violin.ntrials4,
         height = 16, width = 16, units = "cm", dpi = 1200)
  
# ---- calculate all correlations for the appendix ----
funcRC(data)