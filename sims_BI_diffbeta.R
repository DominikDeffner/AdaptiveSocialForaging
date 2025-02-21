
#Social foraging (full) Bayesian Inference simulations

#Set working directory
setwd("~/GitHub/AdaptiveSocialForaging")

beta_mean <- function(a,b){
  return(a/(a+b))
}

beta_sd <- function(a,b){
  return(sqrt( (a*b)/ ( (a+b)^2 * (a +b+1) )) )
}


#Numerical simulation of difference between two beta distributions
diff_beta <- function(a1,b1,a2,b2, N_sample){
  p1 <- rbeta(N_sample, a1, b1)
  p2 <- rbeta(N_sample, a2, b2)
  return(p2-p1)
}


#Each second, agents decide between different patches to forage at

Sim_fct <- function(Tmax  = 90,        #Seconds per round
                    N      = 5,        #Number of foragers
                    N_options = 2,     #Number of patches
                    max_catch = 0.5,   #Maximum catch probability
                    pond_ratio = 0.95  #Catch prob ratios
                    )       
  {    
  
  p_catch    <- max_catch*pond_ratio^(0:(N_options-1))    #Assign initial catch probabilities for each lake
  Beta_counts<- array(1, dim = c(N,N_options,2 ) )        #Counts for beta distribution
  
  patch      <- sample(1:N_options, N, replace = TRUE)    #Current patch
  
  #Create output object to record choices and payoffs participants received as well as beta counts for belief distribution
  
  Result <- list(id=rep(1:N, each  = Tmax),
                 time = rep(1:Tmax, N),
                 Patch=NA,
                 Payoff=NA,
                 Best = NA,
                 beta_counts = array(NA, dim = c(Tmax, N, N_options,2)) )
  
  # Start simulation loop
  for (t in 1:Tmax) {
    
    #Store expectations
    Result$beta_counts[t,,,] <- Beta_counts
    
    #Record best patch
    Result$Best[which(Result$time==t)]  <- which.max(p_catch)
    
    #Loop over all individuals
    for (id in 1:N){
      
        #Record patch
        Result$Patch[which(Result$id==id & Result$time==t)]  <- patch[id]

        #Catch fish with certain probability  
        payoff <- rbinom(1, 1, p_catch[patch[id]])
          
        #Record Payoff
        Result$Payoff[which(Result$id==id & Result$time==t)]  <- payoff
        
        #Update counts based on payoff
        if (payoff ==1){
          Beta_counts[id,patch[id],1] <- Beta_counts[id,patch[id],1] + 1
        }else{
          Beta_counts[id,patch[id],2] <- Beta_counts[id,patch[id],2] + 1
        }
          
        #Agents switch probabilistically based on the posterior distribution of the difference between both patches
        p_diff <- diff_beta(Beta_counts[id,1,1],
                            Beta_counts[id,1,2],
                            Beta_counts[id,2,1],
                            Beta_counts[id,2,2], 
                            1e3)
        choice <- ifelse(sample(p_diff, 1)>0, 2, 1)

      #Update patch
      patch[id] <- choice
      
    }#individual id
    
  }#t
  
  return(Result)
  
}#sim_funct

#Run simulation

Result <- Sim_fct()


#color stuff
require(RColorBrewer)#load package
library(scales)
col.pal <- brewer.pal(4, "Dark2")

par(mfrow = c(5,5), 
    mar = c(1,1,1,1), 
    oma = c(2,1,1,0))

for (id in 1:5) {
  idx <- Result$id==id
  
  plot(Result$Patch[idx], ylim = c(1, 2), col = col.pal[1+ Result$Patch[idx]], pch = ifelse(Result$Payoff[idx]==1,16,1), ylab = "Patch")
  text(sum(Result$Payoff[idx]), x = 90, y= 1.25)
  if ( id == 1) mtext("Choices", side = 3,line = 0)

  #Expected values
  plot(beta_mean( Result$beta_counts[,id,1,1],  Result$beta_counts[,id,1,2] ), type = "l", pch = 1, ylim = c(0,1), ylab = "EV", col = col.pal[2])
  lines(beta_mean( Result$beta_counts[,id,2,1],  Result$beta_counts[,id,2,2] ), col = col.pal[3])
  mtext("Time [s]", side = 1,line = 1, outer = TRUE)
  if ( id == 1) mtext("EV", side = 3,line = 0)
  
  #Uncertainties
  plot(beta_sd( Result$beta_counts[,id,1,1],  Result$beta_counts[,id,1,2] ), type = "l", pch = 1, ylim = c(0,0.3), ylab = "SD", col = col.pal[2])
  lines(beta_sd( Result$beta_counts[,id,2,1],  Result$beta_counts[,id,2,2] ), col = col.pal[3])
  if ( id == 1) mtext("SD", side = 3,line = 0)
  
  #Beliefs over time
  for (t in 1:90) {
    curve(dbeta(x, Result$beta_counts[t,id,1,1], Result$beta_counts[t,id,1,2]), from = 0, to = 1, ylab = "Belief", xlab = "P(catch)", yaxt = "n", add = ifelse(t==1, FALSE, TRUE), col = alpha("black", alpha = 0.1+ (t/100)), ylim = c(0,10))
  }
  if ( id == 1) mtext("Belief lake I", side = 3,line = 0)
  
  for (t in 1:90) {
    curve(dbeta(x, Result$beta_counts[t,id,2,1], Result$beta_counts[t,id,2,2]), from = 0, to = 1, ylab = "Belief", xlab = "P(catch)",  yaxt = "n", add = ifelse(t==1, FALSE, TRUE), col = alpha("black", alpha = 0.1+ (t/100)), ylim = c(0,10))
  }  
  if ( id == 1) mtext("Belief lake II", side = 3,line = 0)
  
  
}

