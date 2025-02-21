
#Adaptive social foraging simulations

#Set working directory
setwd("~/GitHub/AdaptiveSocialForaging")
source("functions.R")

#Each second, agents decide between different patches to forage at

#Social foragers

Sim_fct <- function(Tmax  = 270,        #Seconds per round
                    N      = 6,        #Number of foragers
                    N_options = 3,     #Number of patches
                    t_travel  = 2,     #Travel time between patches
                    max_catch = 0.8,   #Maximum catch probability
                    pond_ratio = 0.7, #Catch prob ratios
                    alphaG = 0.2,      #Learning rate for gains
                    alphaL = 0.4,      #Learning rate for losses
                    betaQ  = 20,        #Inverse temperature
                    betaC  = 2,       #Choice autocorrelation
                    Env_change = TRUE,  #Does the environment change as in experiment, 
                    sigma = 0.1, 
                    theta = 2,       #Conformity exponent
                    social_int = "DB"  #Decision-biasing ("DB"), value-shaping ("VS") or no ("NO") social influence
                    )       
  {    
  
  p_catch    <- max_catch*pond_ratio^(0:(N_options-1))    #Assign initial catch probabilites for each lake
  Q_values   <- matrix(max_catch*pond_ratio, nrow = N, ncol = N_options )    #Latent value estimates
  R_counter  <- rep(0, N)                                 #Traveling or not
  patch      <- sample(1:N_options, N, replace = TRUE)    #Current or previous patch (in case agent is travelling)
  
  #Create output object to record choices and payoffs participants received as well as observed choices and latent Q_values
  
  Result <- list(id=rep(1:N, each  = Tmax), time = rep(1:Tmax, N), Traveling = NA, Patch=NA, Payoff=NA, Best = NA, Q_values = array(NA, dim = c(Tmax, N, N_options)), Choice_probs = array(NA, dim = c(Tmax, N, N_options)))
  
  #Times of change
  times_of_change <- sapply(1:2, function(t) round(runif(1, t*90-10, t*90+10 )  ) )
  
  # Start simulation loop
  for (t in 1:Tmax) {
    
    #Store Q values
    Result$Q_values[t,,] <- Q_values
    
    #Record best patch
    Result$Best[which(Result$time==t)]  <- which.max(p_catch)
    
    #Social information
    numbers_at_patch <- sapply(1:N_options, function(i) length( which(patch[R_counter == 0]==i) ) )
    
    #Loop over all individuals
    for (id in 1:N){
      
      #Check if individual is currently relocating
      if (R_counter[id] > 0){
        Result$Traveling[Result$id==id & Result$time==t] <- 1
        Result$Patch[which(Result$id==id & Result$time==t)]  <- 0
        Result$Payoff[which(Result$id==id & Result$time==t)]  <- 0
        
        R_counter[id] <- R_counter[id]-1
      
        } else {
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
        
        #Update Q values based on payoff
        pay <- rep(0, N_options)
        pay[patch[id]] <- payoff
          
        #Updating only chosen option
        Q_values[id,patch[id]] <- Q_values[id,patch[id]] + ifelse(payoff==1, alphaG, alphaL) * (payoff-Q_values[id,patch[id]])
        
        #In case of value-shaping social influence, we update Q-values based on location of others
        if (social_int == "VS" & sum(others_at_patch > 0)){
          Q_values[id,] <- Q_values[id,] + sigma * ((others_at_patch^theta)/sum(others_at_patch^theta) -  Q_values[id,])
        }
        
       #Individual choice probabilities based on values and choice autocorrelation
       #Construct variable for choice trace  
       C <- rep(0, N_options)
       C[patch[id]] <- 1
      
       p = softmax(betaQ * Q_values[id, ] + betaC * C)
       
       #In case of decision-biasing social influence, we update policy based on location of others
       if (social_int == "DB" & sum(others_at_patch > 0)){
         p <- p + sigma * ((others_at_patch^theta)/sum(others_at_patch^theta) -  p)
       }
       
       #Record choice probabilities
       Result$Choice_probs[t, id, ] <- p
      
      #Make choices
      choice <- sample(1:N_options, size = 1, prob = p)
      
      if (choice != patch[id]){
          R_counter[id] <- t_travel
      }
      
      #Update patch
      patch[id] <- choice
      
      }#else
      
    }#individual id
    
    #Change catch probabilities (randomly rotate clockwise or anti-clockwise)
    if (Env_change == TRUE & t %in% times_of_change){
      rot_dir <- sample(c(1,-1),1)
      p_catch <- c(tail(p_catch, -rot_dir), head(p_catch, rot_dir))
    } 
      
    
  }#t
  
  return(list (Result = Result, 
         times_of_change = times_of_change))
  
}#sim_funct

#Run simulation

Result <- Sim_fct()




#color stuff
require(RColorBrewer)#load package
col.pal <- brewer.pal(4, "Dark2")

#We have different plotting scripts depending on whether we're only running a single simulation
no_agent <- Result[[1]]$id |> unique() |> length()

if( no_agent == 1 ){

  par(mfrow = c(2,1), 
      mar = c(1.5,4,1,1), 
      oma = c(2,0,0,0))
  
plot(Result[[1]]$Patch, ylim = c(0, 3.25), col = col.pal[1+ Result[[1]]$Patch], pch = ifelse(Result[[1]]$Payoff==1,16,1), ylab = "Patch")
abline(v = Result$times_of_change, lty = 2, col = "grey")

segments(0,3.25,Result$times_of_change[1], 3.25, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[1:Result$times_of_change[1]]) ] )
segments(Result$times_of_change[1],3.25,Result$times_of_change[2], 3.25, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[(Result$times_of_change[1]+1):Result$times_of_change[2]]) ] )
segments(Result$times_of_change[2],3.25,270, 3.25, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[(Result$times_of_change[2]+1):270]) ] )


plot(Result[[1]]$Q_values[,1,1], type = "l", ylim = c(0,1), ylab = "Q values", col = col.pal[2])
lines(Result[[1]]$Q_values[,1,2], col = col.pal[3])
lines(Result[[1]]$Q_values[,1,3], col = col.pal[4])
mtext("Time [s]", side = 1,line = 0.75, outer = TRUE)
abline(v = Result$times_of_change, lty = 2, col = "grey")

} else if (no_agent == 6){
  
  par(mfrow = c(6,3), 
      mar = c(1.5,4,1,1), 
      oma = c(2,0,0,0))
  
  for (id in 1:no_agent) {
  idx <- Result[[1]]$id==id
  
  plot(Result[[1]]$Patch[idx], ylim = c(0, 3.25), col = col.pal[1+ Result[[1]]$Patch[idx]], pch = ifelse(Result[[1]]$Payoff[idx]==1,16,1), ylab = "Patch")
  abline(v = Result$times_of_change, lty = 2, col = "grey")
  text(sum(Result[[1]]$Payoff[idx]), x = 270, y= 0.25)
  if (id == 1){
  segments(0,3.25,Result$times_of_change[1], 3.25, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[1:Result$times_of_change[1]]) ] )
  segments(Result$times_of_change[1],3.25,Result$times_of_change[2], 3.25, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[(Result$times_of_change[1]+1):Result$times_of_change[2]]) ] )
  segments(Result$times_of_change[2],3.25,270, 3.25, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[(Result$times_of_change[2]+1):270]) ] )
  }#if id == 1
  
  plot(Result[[1]]$Q_values[,id,1], type = "l", ylim = c(0,1), ylab = "Q values", col = col.pal[2])
  lines(Result[[1]]$Q_values[,id,2], col = col.pal[3])
  lines(Result[[1]]$Q_values[,id,3], col = col.pal[4])
  mtext("Time [s]", side = 1,line = 0.75, outer = TRUE)
  abline(v = Result$times_of_change, lty = 2, col = "grey")
  if (id == 1){
    segments(0,1,Result$times_of_change[1], 1, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[1:Result$times_of_change[1]]) ] )
    segments(Result$times_of_change[1],1,Result$times_of_change[2], 1, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[(Result$times_of_change[1]+1):Result$times_of_change[2]]) ] )
    segments(Result$times_of_change[2],1,270, 1, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[(Result$times_of_change[2]+1):270]) ] )
  }#if id == 1
  
  plot(Result[[1]]$Choice_probs[,id,1], type = "p", pch = 1, ylim = c(0,1), ylab = "P(choice)", col = col.pal[2])
  points(Result[[1]]$Choice_probs[,id,2], pch = 1, col = col.pal[3])
  points(Result[[1]]$Choice_probs[,id,3], pch = 1, col = col.pal[4])
  mtext("Time [s]", side = 1,line = 0.75, outer = TRUE)
  abline(v = Result$times_of_change, lty = 2, col = "grey")
  if (id == 1){
    segments(0,1,Result$times_of_change[1], 1, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[1:Result$times_of_change[1]]) ] )
    segments(Result$times_of_change[1],1,Result$times_of_change[2], 1, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[(Result$times_of_change[1]+1):Result$times_of_change[2]]) ] )
    segments(Result$times_of_change[2],1,270, 1, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[(Result$times_of_change[2]+1):270]) ] )
  }#if id == 1
  
  }
  
  
  
} else{
  
  par(mfrow = c(3,1), 
      mar = c(1.5,4,1,1), 
      oma = c(2,0,0,0))
  
  for (id in 1:no_agent) {
    idx <- Result[[1]]$id==id
    
    plot(Result[[1]]$Patch[idx], ylim = c(0, 3.25), col = col.pal[1+ Result[[1]]$Patch[idx]], pch = ifelse(Result[[1]]$Payoff[idx]==1,16,1), ylab = "Patch")
    abline(v = Result$times_of_change, lty = 2, col = "grey")
    
    if (id == 1){
      segments(0,3.25,Result$times_of_change[1], 3.25, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[1:Result$times_of_change[1]]) ] )
      segments(Result$times_of_change[1],3.25,Result$times_of_change[2], 3.25, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[(Result$times_of_change[1]+1):Result$times_of_change[2]]) ] )
      segments(Result$times_of_change[2],3.25,270, 3.25, lwd = 5, col = col.pal[1+ unique(Result[[1]]$Best[(Result$times_of_change[2]+1):270]) ] )
    }#if id == 1
    
  }
  
  
}


