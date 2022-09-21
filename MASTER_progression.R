#====================================== eGFR-C ECONOMIC MODEL 2021======================================#

# This *draft* model and associated files were developed by Bethany Shinkins and Alison Smith from the University of Leeds as part of a NIHR funded project - eGFR-C. It was developed in anticipation of the study requiring a full economic evaluation to evaluate the cost-effectiveness of using different eGFR equations to monitor disease progression. As explained in the report, this evaluation was not required in the end, but we are sharing our code to date in case of use for future evaluations.#
# This code will need to be further developed and adapted in line with the researcher's research question. The references for the model parameters have been cited within the code or report, but these are likely to require updating based on the latest evidence.#
#A full systematic review (up to 2021) of economic evaluations focused on the cost-effectiveness of monitoring for kidney progression can be found in the main report (link to be added once published).#
#All code will need to be thoroughly checked for errors. The authors, the University of Leeds, and NIHR hold no responsibility or liability for any errors in this code#
#If you use this code, please acknowledge by citing the HTA report (link to be added)#


#====== Set working directory and load packages=======#

# set working directory (location where model files are stored: alter as required)
setwd("") #insert the location where model files are stored

# upload required library packages. You may need to install the packages if you don't have them already e.g. install.packages("gtools", type="source") 
library(MASS)     
library(gtools)  #install.packages("gtools", type="source") 
library(plyr)   


#=========================================================================================================================#
#============================================== ECONOMIC DECISON MODEL ===================================================#
#=========================================================================================================================#

#======= Load Global Model Data =======# 

# remove any stored data
rm(list=ls(all=TRUE))

# set global variables
seed       <- 10             # Set seed for random number generator. Keeping the same seed throughout analyses ensures any differences between results are not due to stochastic variation. 
disc       <- 0.035          # Discount rate for future costs and benefits (as per NICE 2013 methods guide)
Nsim       <- 10000          # Number of Monte Carlo simulations. Reduce to 1 or 10 when testing code, then 10000 in final version.    
threshold  <- 20000          # NICE willingness to pay per QALY threshold (NICE 2013 methods guide)

# load parameters and model files
source("Parameters_progression delay.R")  

# Draw model parameters Nsim times, to use in the baseline model     
parameter_draw()

#Store parameters
Params <<- Params
#write.csv(Params, "Params.csv")


#====================== Standard care ('XXX') arm =========================#

#Set intervention parameter to zero
int <- "referred"

# Load model  
source("Model progression delay.R")  

# Run model
sim.base <- array(NA,c(Nsim,2))                  # Define array table in which costs and QALYs from model will be stored 
for (i in 1:Nsim) {                              # Run model Nsim times 
  sim.base[i,]     <- model(i, int)     
}
costs.base  <- sim.base[,2]                      # Store the cost results           
QALYs.base  <- sim.base[,1]                      # Store the QALY results    

#Check costs and QALYs
head(costs.base)
head(QALYs.base)

#Store tps
tps_base <<- tps
#write.csv(tps_base, "tps_base.csv")

#Check trace validity (rowSums should all add to 100)        
trace.base   <- apply(trace*100, c(1,2), mean)   #*100 to report full % rather than proportions (avoid e+ scientific notation)    
rowSums(trace.base) 
#write.csv(trace.base, "trace_base.csv")

#Save results for lifetime analysis
CEresults_base <- sim.base[,1:2]


#Plot trace over time (use as a validity check) - NEEDS EDITING
# 
# plot(rowSums(trace.base[,c(1:8,46,47)]), ylim=c(0, 100), xlim=c(0,60), col="black")   #
# lines(rowSums(trace.base[,9:13]), col="red")                                          #
# lines(rowSums(trace.base[,14:35]), col="green")                                       #         
# lines(rowSums(trace.base[,c(36:39, 48:51)]), col="blue")                              #
# lines(rowSums(trace.base[,40:43]), col="orange")                                      #
# lines(rowSums(trace.base[,44:45]), col="purple")                                      #
# lines(rowSums(trace.base[,c(3,4,11,15, 18, 22, 26, 29, 33)]), col="deepskyblue")      #
# lines(rowSums(trace.base[,c(6,7,13,16, 19, 24, 27, 30, 35)]), col="magenta")          #
#plot(trace.base[T,])
#plot(trace.base[,36], col="green")


#====================== Intervention ('XXX') arm =========================#

#Set test value
int  <- 1

#Run model
sim.int <- array(NA,c(Nsim,2))
for (i in 1:Nsim) {
  sim.int[i,]     <- model(i, int)
}
costs.int <- sim.int[,2]
QALYs.int <- sim.int[,1]

head(costs.int)
head(QALYs.int)

tps_int <<-tps
#write.csv(tps_int, "tps_int.csv")

trace.int   <- apply(trace*100, c(1,2), mean)     
rowSums(trace.int)  
#write.csv(trace.int, "trace_int.csv")

CEresults_int <-  sim.int[,1:2]  

#Plot trace over time (use as a validity check)
# 
# plot(rowSums(trace.int[,c(1:8,46,47)]), ylim=c(0, 100), xlim=c(0,60), col="black")   #
#lines(rowSums(trace.int[,c(1:8,46,47)]), col="red")
# lines(rowSums(trace.int[,9:13]), col="red")                                          #
# lines(rowSums(trace.int[,14:35]), col="red")                                         #          
# lines(rowSums(trace.int[,c(36:39, 48:51)]), col="red")                               #
# lines(rowSums(trace.int[,40:43]), col="red")                                         #
# lines(rowSums(trace.int[,44:45]), col="red")                                         #
# lines(rowSums(trace.int[,c(3,4,11,15, 18, 22, 26, 29, 33)]), col="red")              #
# lines(rowSums(trace.int[,c(6,7,13,16, 19, 24, 27, 30, 35)]), col="red")              #


#Compare tps- check differences in transitions across arms makes sense 
#tps_base[51,,1]
#tps_int[37,,1]
#sum(tps_base[51,,1])

#Compare traces
#trace_final  <- data.frame(trace.base[T,], trace.int[T,])
#plot(trace_final[,2]-trace_final[,1])



#=========================================================================================================================#
#===================================================  ANALYSIS ===========================================================#
#=========================================================================================================================#

#============ INTERVENTION VS. STANDARD CARE ===========#
threshold  <- 20000
NB.base    <- CEresults_base[,1] - CEresults_base[,2]/threshold
NB.int     <- CEresults_int[,1] - CEresults_int[,2]/threshold
INB_int    <- NB.int-NB.base
mean(INB_int)

#=== Compile and save results ===# 
int_results       <-  data.frame(CEresults_base, CEresults_int, NB.base, NB.int, INB_int)
int_results[,8]   <-  int_results[,3]- int_results[,1]
int_results[,9]   <-  int_results[,4]- int_results[,2]
int_results[,10]  <-  ifelse(int_results[,8]>0,1,0)
int_results[,11]  <-  ifelse(int_results[,9]<0,1,0)
int_results[,12]  <-  ifelse(int_results[,7]>0,1,0)
colnames(int_results) <- c("QALYs_base","Costs_base","QALYs_int","Costs_int","NB_base","NB_int","INB", "Incr.QALY", "Incr.Cost", "p_eff", "p_save", "p_costeff")
write.csv(int_results, "int.csv")
mean(int_results$Incr.QALY)
mean(int_results$Incr.Cost)
mean(int_results$Incr.Cost)/mean(int_results$Incr.QALY)

mean(int_results$QALYs_base)
mean(int_results$Costs_base)
mean(int_results$QALYs_int)
mean(int_results$Costs_int)

#==== 95% CIs ===#
Conf_ints <- array(NA, c(length(int_results[1,]),2))
for (i in 1:(length(int_results[1,]))){
  Conf_ints[i,1] <- quantile(int_results[,i],c(0.025))
  Conf_ints[i,2] <- quantile(int_results[,i], c(0.975))
}
write.csv(Conf_ints, "Conf_ints_basecase1.csv")

#=== Traces ====#
# write.csv(trace, "trace.csv")
# write.csv(trace_int, "trace_int.csv")


#=== Scatter Plot ===#
plot(int_results[c(8,9)], xlim=c(-0.2,0.2), ylim=c(-5000, 5000),  xlab = "Incremental QALY",
     ylab = "Incremental Cost (£)", yaxs="i", col="turquoise3")
abline(a=0, b=20000, col = "blue4", lwd = 2)
abline(a=0, b=0, col = "black", lwd = 0.5)
abline(v=0)
points(x=mean(int_results$Incr.QALY), y=mean(int_results$Incr.Cost),pch=24, bg="white", col="darkorchid", cex=1.4, lwd=2)
labels <- c("£20,000/QALY threshold", "Mean value")
cols <- c("blue4","darkorchid")
legend(-0.185,3000, labels, col = cols, cex = 0.8, lwd=2, lty=c(1,NA), pch=c(NA, 24))

# mean(int_results$Incr.QALY)
# mean(int_results$p_eff)
# mean(int_results$Incr.Cost)
# mean(int_results$p_save)
# mean(int_results$p_costeff)



#==================== CEAF  ======================#                   
OUT <- list()
thresh <- seq(100,200000,100)

for (i in 1:length(thresh)) {
  
  NB.int                   <- int_results$QALYs_int - int_results$Costs_int/thresh[i]  #Note NB in QALYS (i.e. net health benefit, NOT monetary)
  NB.base                  <- int_results$QALYs_base - int_results$Costs_base/thresh[i]  
  maxNB                    <- ifelse(NB.base>= NB.int, NB.base, NB.int)   
  CE.base                  <- ifelse(NB.base == maxNB, 1, 0)
  CE.int                   <- ifelse(NB.int == maxNB, 1, 0)
  prob.CE.base             <- mean(CE.base) 
  prob.CE.int              <- mean(CE.int) 
  mean.NB.base             <- mean(NB.base)
  mean.NB.int              <- mean(NB.int)
  mean.NB.max              <- max(mean.NB.base,mean.NB.int)    
  ceaf                     <- ifelse(mean.NB.max==mean.NB.base, prob.CE.base,prob.CE.int)
  
  OUT[[i]] <- list(NB.base=NB.base,NB.int=NB.int,maxNB=maxNB,
                   CE.base=CE.base,CE.int=CE.int, prob.CE.base=prob.CE.base, prob.CE.int=prob.CE.int,
                   mean.NB.base = mean.NB.base, mean.NB.int= mean.NB.int,ceaf=ceaf)
}

CEAF <- data.frame(threshold=thresh)
for (i in 1:length(thresh)){
  CEAF$Probability.CE[i]  <- OUT[[i]]$ceaf
  CEAF$CE.base[i]         <- OUT[[i]]$prob.CE.base
  CEAF$CE.int[i]       <- OUT[[i]]$prob.CE.int
  }


plot(CEAF$threshold,CEAF$CE.base,type="l", ylim = c(0,1), xlim = c(0,50000), xaxp  = c(0, 50000, 5), col = "blue4", lwd=2, xlab = "Willingness-to-pay threshold (£ per QALY)", ylab = "Probability cost-effective")
lines(CEAF$threshold,CEAF$CE.int,type="l", ylim = c(0,1), xlim = c(0,50000), col = "turquoise3",lwd=2)
lines(CEAF$threshold[seq(1,500,6)],CEAF$Probability.CE[seq(1,500,6)],type="p", ylim = c(0,1))
labels <- c("Standard Care","Standard Care + int", "Cost-effectiveness frontier")
cols <- c("blue4","turquoise3", "black")
legend(26000,0.60,labels, col = cols, lty = c(1,1,NA), pch=c(NA, NA, 1), cex = 1, lwd=c(2,2,1))
write.csv(CEAF, "CEAF_int_basecase1.csv", row.names=FALSE)
