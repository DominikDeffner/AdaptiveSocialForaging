
#Social foraging (full) Bayesian Inference simulations

#Set working directory
setwd("~/GitHub/AdaptiveSocialForaging")
source("functions.R")

beta_mean <- function(a,b){
  return(a/(a+b))
}

beta_sd <- function(a,b){
  return(sqrt( (a*b)/ ( (a+b)^2 * (a +b+1) )) )
}

beta_UCB <- function(a,b,beta){
  E <- beta_mean(a,b)
  SD <- beta_sd(a,b)
  return(E + beta*SD)
}

prob_2_better <- function(a1,b1,a2,b2){
  return(sum(sapply(0:(a2-1), function(i) (beta(a1+i, b1+b2))/( (b2+i)*beta(1+i,b2)*beta(a1,b1) ) ) ))
}



#Each second, agents decide between different patches to forage at

#Social foragers

Sim_fct <- function(Tmax  = 90,        #Seconds per round
                    N      = 5,        #Number of foragers
                    N_options = 2,     #Number of patches
                    max_catch = 0.6,   #Maximum catch probability
                    pond_ratio = 0.8,  #Catch prob ratios
                    beta       = 0.2,  #Uncertainty-directed exploration
                    lambda  = 6,        #Inverse temperature
                    tau  = 3,         #Choice autocorrelation
                    sigma = 0.3, 
                    theta = 3,       #Conformity exponent
                    social_int = "DB"  #Decision-biasing ("DB"), value-shaping ("VS") or no ("NO") social influence
                    )       
  {    
  
  p_catch    <- max_catch*pond_ratio^(0:(N_options-1))    #Assign initial catch probabilities for each lake
  Beta_counts<- array(1, dim = c(N,N_options,2 ) )        #Counts for beta distribution
  patch      <- sample(1:N_options, N, replace = TRUE)    #Current or previous patch (in case agent is travelling)
  
  #Create output object to record choices and payoffs participants received as well as observed choices and latent Q_values
  
  Result <- list(id=rep(1:N, each  = Tmax), time = rep(1:Tmax, N), Patch=NA, Payoff=NA, Best = NA, beta_counts = array(NA, dim = c(Tmax, N, N_options,2)), Choice_probs = array(NA, dim = c(Tmax, N, N_options)))

  # Start simulation loop
  for (t in 1:Tmax) {
    
    #Store expectations
    Result$beta_counts[t,,,] <- Beta_counts
    
    #Record best patch
    Result$Best[which(Result$time==t)]  <- which.max(p_catch)
    
    #Social information
    numbers_at_patch <- sapply(1:N_options, function(i) length( which(patch==i) ) )
    
    #Loop over all individuals
    for (id in 1:N){
      
        #Record patch
        Result$Patch[which(Result$id==id & Result$time==t)]  <- patch[id]
        Result$Traveling[Result$id==id & Result$time==t] <- 0
      
        #Social info
        others_at_patch <- numbers_at_patch
        others_at_patch[patch[id]] <- others_at_patch[patch[id]] - 1
        
        #Catch fish with certain probability  
        payoff <- rbinom(1, 1, p_catch[patch[id]])
          
        #Record Payoff
        Result$Payoff[which(Result$id==id & Result$time==t)]  <- payoff
        
        #Update counts based on payoff and calculate UCBs
        if (payoff ==1){
          Beta_counts[id,patch[id],1] <- Beta_counts[id,patch[id],1] + 1
        }else{
          Beta_counts[id,patch[id],2] <- Beta_counts[id,patch[id],2] + 1
        }
          
        UCBs <- sapply(1:N_options, function(x) beta_UCB(Beta_counts[id,x,1],Beta_counts[id,x,2], beta ) )
        
        #In case of value-shaping social influence, we update Q-values based on location of others
        if (social_int == "VS" & sum(others_at_patch > 0)){
          Q_values[id,] <- Q_values[id,] + sigma * ((others_at_patch^theta)/sum(others_at_patch^theta) -  Q_values[id,])
        }
        
       #Individual choice probabilities based on values and choice autocorrelation
       #Construct variable for choice trace  
       C <- rep(0, N_options)
       C[patch[id]] <- 1
      
       p = softmax(lambda * UCBs + tau * C)
       
       #In case of decision-biasing social influence, we update policy based on location of others
       if (social_int == "DB" & sum(others_at_patch > 0)){
         p <- p + sigma * ((others_at_patch^theta)/sum(others_at_patch^theta) -  p)
       }
       
       #Record choice probabilities
       Result$Choice_probs[t, id, ] <- p
      
      #Make choices
      choice <- sample(1:N_options, size = 1, prob = p)
      
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
col.pal <- brewer.pal(4, "Dark2")

  par(mfrow = c(5,4), 
      mar = c(1.5,4,1,1), 
      oma = c(2,0,0,0))
  
  for (id in 1:5) {
  idx <- Result$id==id
  
  plot(Result$Patch[idx], ylim = c(1, 2), col = col.pal[1+ Result$Patch[idx]], pch = ifelse(Result$Payoff[idx]==1,16,1), ylab = "Patch")
  text(sum(Result$Payoff[idx]), x = 90, y= 1.25)

#Expected values
  plot(beta_mean( Result$beta_counts[,id,1,1],  Result$beta_counts[,id,1,2] ), type = "l", pch = 1, ylim = c(0,1), ylab = "EV", col = col.pal[2])
  lines(beta_mean( Result$beta_counts[,id,2,1],  Result$beta_counts[,id,2,2] ), col = col.pal[3])
  mtext("Time [s]", side = 1,line = 0.75, outer = TRUE)
  
  #Uncertainties
  plot(beta_sd( Result$beta_counts[,id,1,1],  Result$beta_counts[,id,1,2] ), type = "l", pch = 1, ylim = c(0,0.3), ylab = "SD", col = col.pal[2])
  lines(beta_sd( Result$beta_counts[,id,2,1],  Result$beta_counts[,id,2,2] ), col = col.pal[3])
  mtext("Time [s]", side = 1,line = 0.75, outer = TRUE)

  plot(Result$Choice_probs[,id,1], type = "l", pch = 1, ylim = c(0,1), ylab = "P(choice)", col = col.pal[2])
  lines(Result$Choice_probs[,id,2], col = col.pal[3])
  mtext("Time [s]", side = 1,line = 0.75, outer = TRUE)

  }
  
  
  


