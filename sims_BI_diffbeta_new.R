
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

#Analytical solution for the probability that patch 2 is better
#Following the derivation found here:https://www.evanmiller.org/bayesian-ab-testing.html

prob_2_better <- function(a1,b1,a2,b2){
  return(sum(sapply(0:(a2-1), function(i) (beta(a1+i, b1+b2))/( (b2+i)*beta(1+i,b2)*beta(a1,b1) ) ) ))
}

#Each second, agents decide between different patches to forage at

Sim_fct <- function(Tmax  = 75,        #Seconds per round
                    N      = 1,        #Number of foragers
                    N_options = 2,     #Number of patches
                    max_catch = 0.7,   #Maximum catch probability
                    pond_ratio = 0.65  #Catch prob ratios
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
        # p_diff <- diff_beta(Beta_counts[id,1,1],
        #                     Beta_counts[id,1,2],
        #                     Beta_counts[id,2,1],
        #                     Beta_counts[id,2,2], 
        #                     1e3)
        # choice <- ifelse(sample(p_diff, 1)>0, 2, 1)
        
        #Alternatively, agents infer which patch is better 
        p2 <- prob_2_better(Beta_counts[id,1,1],
                            Beta_counts[id,1,2],
                            Beta_counts[id,2,1],
                            Beta_counts[id,2,2])
        
        choice <- sample(c(1,2), 1, prob = c(1-p2,p2))
        
      #Update patch
      patch[id] <- choice
      
    }#individual id
    
  }#t
  
  return(Result)
  
}#sim_funct

#Run simulation

Result <- Sim_fct()



graphics.off()
pdf("BayesianAgents.pdf", width = 8, height = 2.5)
par(mfrow = c(1,4), 
    mar = c(1,1,1,1), 
    oma = c(2.75,1,1.25,0))

plot(Result$Patch, ylim = c(0.5, 2.5), type = "n", yaxt = "n", pch = ifelse(Result$Payoff==1,16,1), xlab = "", ylab = "")
segments(1:75, Result$Patch-0.25, 1:75, Result$Patch+0.25, col = ifelse(Result$Payoff==1,"black","grey"))

axis(side = 2, at = 1:2, labels = c("Lake I", "Lake II"), cex.axis = 1.5)
mtext("Choices", side = 3,line = 0)
mtext("Time (s)", side = 1,line = 2.5)


#Beliefs over time
for (t in 1:75) {
  curve(dbeta(x, Result$beta_counts[t,1,1,1], Result$beta_counts[t,1,1,2]), from = 0, to = 1, ylab = "", xlab = "", yaxt = "n", add = ifelse(t==1, FALSE, TRUE), col = alpha("black", alpha = 0.1+ (t/200)), ylim = c(0,10), lwd = 0.5)
  if(t == 1) abline(v = 0.7, lty = 2, col = "grey")
  }
mtext("Belief Lake I", side = 3,line = 0)
mtext("P(Catch)", side = 1,line = 2.5)

for (t in 1:75) {
  curve(dbeta(x, Result$beta_counts[t,1,2,1], Result$beta_counts[t,1,2,2]), from = 0, to = 1, ylab = "Belief", xlab = "P(catch)",  yaxt = "n", add = ifelse(t==1, FALSE, TRUE), col = alpha("black", alpha = 0.1+ (t/200)), ylim = c(0,10),lwd = 0.5)
  if(t == 1) abline(v = 0.7*0.65, lty = 2, col = "grey")
}

mtext("Belief Lake II", side = 3,line = 0)
mtext("P(Catch)", side = 1,line = 2.5)

p_better <- sapply(1:75, function(t) prob_2_better(Result$beta_counts[t,1,1,1], Result$beta_counts[t,1,1,2],Result$beta_counts[t,1,2,1], Result$beta_counts[t,1,2,2])  )
p_better <- 1-p_better

plot(p_better, type = "l", xlab = "", ylab = "", ylim = c(0,1))
mtext("P(Lake I better)", side = 3,line = 0)
mtext("Time (s)", side = 1,line = 2.5)
dev.off()
