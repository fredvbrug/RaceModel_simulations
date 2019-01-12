# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# script to generate the simulated data
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

rm(list = ls()) # clear everything
library (gamlss) # we will use gamlss library for the simulation of the RT distributions

#------------------------------------------
# define some variables
#------------------------------------------
NSIM <- 1000 # number of simulated 'subjects'
SSD.start <- 300 # start value of SSD
PSIGNALS <- 0.25 # proportion of stop signals
STRIALS <- c(25, 50, 100, 200) # number of stop-signal trials per simulation
NTRIALS <- STRIALS/PSIGNALS

pSIM <- data.frame(NTRIALS = NTRIALS, NSIM = NSIM, SSD.start = SSD.start, PSIGNALS = PSIGNALS)

#------------------------------------------
# determine shape of RT distribution
#------------------------------------------
RT.mu <- 500 # mean of normal part of RT distribution
RT.sdpop <- 50 # to allow within-group differences, mean(simulation) = random value from normal distribution with SD
RT.sigma <- 50 # standard deviation of normal part of RT distribution
RT.tau <- c(1, 50, 100, 150, 200) # mean of exponential part of the RT distribution
RT.miss <- c(0, 5, 10, 15, 20) # percentage missed go responses
RT.cutoff <- 1500 # cut-off value for go RTs

RT.tau_c <- rep(RT.tau, each = length(RT.miss)) # we want co combine every tau with every miss
RT.miss_c <- rep(RT.miss, times = length(RT.tau)) # we want co combine every miss with every tau

pRT <- data.frame (mu = RT.mu, sdpop = RT.sdpop, sigma = RT.sigma, tau = RT.tau_c, miss = RT.miss_c, cutoff = RT.cutoff) # put everything in a single data frame

#------------------------------------------
# determine shape of SSRT distrbution
#------------------------------------------
SSRT.mu = 200 #mean of normal part of RT distribution
SSRT.sdpop = 20 # to allow within-group differences, mean(simulation) = random value from normal distribution with SD
SSRT.sigma = 15 #standard deviation of normal part of SSRT distribution
SSRT.tau = 30 # mean of exponential part of the SSRT distribition

pSSRT <- data.frame (mu = SSRT.mu, sdpop = SSRT.sdpop, sigma = SSRT.sigma, tau = SSRT.tau) # put everything in a single data frame 

#------------------------------------------
# function to run the simulation
#------------------------------------------
simulator <- function (pSIM, pRT, pSSRT) {
  
  simulation <- data.frame() # create empty data frame

	for (s in 1:pSIM$NSIM){
		#reset the SSD value at the beginning of each simulation
		SSD = pSIM$SSD.start +  pRT$tau
	  
    #-----------------------------------------
		# get an ex-gaussian RT distribution 
    #-----------------------------------------
	  # step 1: determine mean of the subject; this mean is allowed to differ from population mean
    # but the lowest possible subject mean = 300 ms
    sub.mean <- 0 
		while (sub.mean < 300)
			sub.mean <- rnorm(1, pRT$mu, pRT$sdpop) # subject mean comes from a normal distribution  
		
    # step 2: get for this subject an ex-gaussian go-RT distribution 
		RT.true <- rexGAUS(pSIM$NTRIALS, mu=sub.mean, sigma=pRT$sigma, nu= pRT$tau) # ex-gaussian RT distribution
		
		# step 3a: check for missed trials... 
		RT.true <- ifelse(RT.true > pRT$cutoff, 9999, RT.true) # label these with 9999
		
		# step 3b: ...and insert extra ones if needed
		RT.true <- sort(RT.true)  # start with putting the slow trials ('original' go failures) at the end 
		missed <- sum(RT.true>pRT$cutoff) # determine how many 'go failures' occur
		extra_missed <- (pSIM$NTRIAL * pRT$miss/100) - missed # determine how many extra 'go failures' we need
		
		if (extra_missed > 0){
		  trial_missed <- sample(1:(pSIM$NTRIAL-missed), extra_missed) # get a random selection of trials that we will turn into a 'go failure'
		  RT.true[trial_missed] <- 9999 # replace the values for these trials
		}
		
    # step 4: randomize the RTs again
		RT.true <- sample(RT.true)
		
    #-----------------------------------------
		# get an ex-gaussian SSRT distribution 
    #-----------------------------------------
		# step 1: determine mean of the subject; this mean is allowed to differ from population mean
		# but the lowest possible subject mean = 100 ms
    sub.mean <- 0
		while (sub.mean < 100)
			 sub.mean <- rnorm(1, pSSRT$mu, pSSRT$sdpop)
		
    # step 2: get a ex-gaussian distribution
    SSRT.true <- rexGAUS(pSIM$NTRIALS, mu=sub.mean, sigma=pSSRT$sigma, nu= pSSRT$tau) # ex-gaussian SSRT distribution
		
    #-----------------------------------------
  	# make a stop-signal list & adjust SSRTs
    #-----------------------------------------
    signals <- sample(1:(pSIM$NTRIAL)) # make a stop-signal list
		signals <- ifelse(signals > (pSIM$NTRIAL * pSIM$PSIGNALS), 0, 1)
		SSRT.true = SSRT.true * signals

    #-----------------------------------------
		# do a very simple race simulation for NTRIALS
		# if RT > SSRT + SSD, then signal-inhibit; else signal-respond
		# SSD will be adjusted accordingly
    #-----------------------------------------
		used.ssd <- rep(0, pSIM$NTRIALS) #define/reset variable 'used.ssd'
		outcome <-  rep(0, pSIM$NTRIALS) #define/reset variable 'outcome'
		race <- rep(0, pSIM$NTRIALS)  #define/reset variable 'race'' (difference relative finishing times)

    # start a loop, required for the tracking
		for (i in 1:pSIM$NTRIALS){
			used.ssd[i] <- signals[i] * SSD
	
			if (signals[i] == 0)
				outcome[i] <- "no-signal" 	
			else {
				race[i] <- RT.true[i] - (SSRT.true[i] + SSD)
		
				if (race[i] > 0){	
					outcome[i] <- "s-inhibit"
					SSD <- SSD + 50
				} else{
					outcome[i] <- "s-respond"
					SSD <- SSD - 50
					}
				}
			}

		# adjust RT, put everything in a single data frame & determine simulation/block number
    tmp <- data.frame(signals, RT.true, SSRT.true, used.ssd, outcome)
		tmp$subject <- s
		tmp$trial <- c(1:pSIM$NTRIALS)
 
		# combine the simulation with the other simulations
		simulation <- rbind(simulation, tmp)
		
		# print the progress
		text <- sprintf ("simulation nmbr. %d (NTRIALS: %d; RT: tau = %d, miss = %d)", s, pSIM$NTRIALS, pRT$tau, pRT$miss)
		print(text)
		}
		
	# save the data when all simulations are finished (Rdata for R-users; csv copy for others)
	filename <- sprintf ("./simulated_data/Ntrials_%d_RTtau_%d_RTmiss_%d.Rdata", pSIM$NTRIALS, pRT$tau, pRT$miss)
	save (simulation, file = filename)
	
	# csvname <- sprintf ("./simulated_data/Ntrials_%d_RTtau_%d_RTmiss_%d.csv", pSIM$NTRIALS, pRT$tau, pRT$miss)
	# write.csv (simulation, file = csvname)
}

#------------------------------------------------------------------------------------------
# do the simulations for all possible combinations of parameters
#------------------------------------------------------------------------------------------
for (n in 1:nrow(pSIM))
  for (i in 1:nrow(pRT)){
    for (j in 1:nrow(pSSRT)){
          # check if the data file already exists
          filename <- sprintf ("./simulated_data/Ntrials_%d_RTtau_%d_RTmiss_%d.Rdata", pSIM$NTRIALS[n], pRT$tau[i], pRT$miss[i])
          exist <- try(load(filename), silent = TRUE)

          # if the data file doesn't exist yet, do the simulation
          if (class(exist) == "try-error")
            simulator(pSIM[n,], pRT[i,], pSSRT[j,])
    }
}
