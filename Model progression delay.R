
#======================= eGFR-C 2021 ==================================================#
# This *draft* model and associated files were developed by Bethany Shinkins and Alison Smith from the University of Leeds as part of a NIHR funded project - eGFR-C. It was developed in anticipation of the study requiring a full economic evaluation to evaluate the cost-effectiveness of using different eGFR equations to monitor disease progression. As explained in the report, this evaluation was not required in the end, but we are sharing our code to date in case of use for future evaluations.#
# This code will need to be further developed and adapted in line with the researcher's research question. The references for the model parameters have been cited within the code or report, but these are likely to require updating based on the latest evidence.#
#A full systematic review (up to 2021) of economic evaluations focused on the cost-effectiveness of monitoring for kidney progression can be found in the main report (link to be added once published).#
#All code will need to be thoroughly checked for errors. The authors, the University of Leeds, and NIHR hold no responsibility or liability for any errors in this code#
#If you use this code, please acknowledge by citing the HTA report (link to be added)#


# Model health states:        1=  CKD Stage 3a
#                             2=  CKD Stage 3b
#                             3=  CKD Stage 4
#                             4=  CKD Stage 5
#                             5=  ESRD dialysis
#                             6=  ESRD Kidney Transplant
#                             7=  Death (all cause)
#                             8=  Spare (Leave these in case need additional health states)
#                             9=  Spare
#                             10= Spare


# Global parameters of model 
S  <- 7                                # Number of health states in model = 7
H  <- 49 #110-startage                 # Time Horizon = 110 years - startage 
T  <- H + 1                            # Number of model cycles. Add one for half-cycle correction.                   


#create empty objects
tps   <- array(NA, dim=c(S,S,T))        # Model Transition Probability Matrix (defines probability of moving between health states over time)
trace <- array(NA, dim=c(T,S,Nsim))     # Trace: to record proportion of cohort in each health state over time, for each simulation
qtime <- rep(NA,len=T)                  # Record QALYs over time 
cost  <- rep(NA,len=T)                  # Record costs over time 


model <- function(i, int) {                                

  
#======================= TRANSITION PROBABILITY MATRIX=======================================#
  
  for (t in 1:T) {   
    
    #Health State 1=  CKD Stage 3a     [Remain Stage 3a, progress to Stage 3b, die]
    tps[1,1,t] <-0  #see at end
    tps[1,2,t] <-if(int==0) {pt_3ato3b[i]*HR_stage3to4} else {pt_3ato3b[i]}  #have applied general stage 3 reduction - needs checking
    tps[1,3,t] <-0          
    tps[1,4,t] <-0
    tps[1,5,t] <-0
    tps[1,6,t] <-if(deathother[startage+t,2] > pt_3todeath[i]) {deathother[startage+t,2]} else {pt_3todeath[i]}
    tps[1,7,t] <-0                                                          
    #tps[1,8,t] <-0
    #tps[1,9,t] <-0
    #tps[1,10,t] <-0
    tps[1,1,t] <-1-sum(tps[1,,t])
   
    #Health state 2= CKD Stage 3b   [Remain Stage 3b, progress to Stage 4, die]
    tps[2,1,t] <-0
    tps[2,2,t] <-0 #see at end
    tps[2,3,t] <-pt_3bto4[i]
    tps[2,4,t] <-0                    
    tps[2,5,t] <-0
    tps[2,6,t] <-if(deathother[startage+t,2] > pt_3todeath[i]) {deathother[startage+t,2]} else {pt_3todeath[i]}  
    tps[2,7,t] <-0                                                          
    #tps[2,8,t] <-0
    #tps[2,9,t] <-0
    #tps[2,10,t] <-0
    tps[2,2,t]  <-1-sum(tps[2,,t])

        
    #Health state 3= CKD Stage 4     [Remain Stage 4, progress to Stage 5, progress to ESRD, die]
    tps[3,1,t] <-0
    tps[3,2,t] <-0
    tps[3,3,t] <-0 #see at end
    tps[3,4,t] <-pt_4to5[i] #ckd5
    tps[3,5,t] <-0 
    tps[3,6,t] <-0
    tps[3,7,t] <-if(deathother[startage+t,2] > pt_4todeath[i]) {deathother[startage+t,2]} else {pt_4todeath[i]}
    #tps[3,8,t] <-0
    #tps[3,9,t] <-0
    #tps[3,10,t] <-0
    tps[3,3,t]<-(1-sum(tps[3,,t]))
  
    #Health state 4=  CKD Stage 5 [Remain Stage 5, progress to ESRD, die] # this needs fixing, adds up to more than 1
    tps[4,1,t] <-0
    tps[4,2,t] <-0
    tps[4,3,t] <-0
    tps[4,4,t] <-0 #see end
    tps[4,5,t] <-pt_5todialysis[i]
    tps[4,6,t] <-if(t<21) {pt_5totransplant[i]} else {0} #need to cap age can go to transplant
    tps[4,7,t] <-(if (deathother[startage+t,2] > pt_5todeath[i]) (deathother[startage+t,2]) else pt_5todeath[i])
    #tps[4,8,t] <-0
    #tps[4,9,t] <-0
    #tps[4,10,t] <-0
    tps[4,4,t] <-1-sum(tps[4,,t])     

      
    # Health state 5= ESRD Dialysis [Remain on dialysis,kidney transplant, die]     
    tps[5,1,t] <-0
    tps[5,2,t] <-0
    tps[5,3,t] <-0
    tps[5,4,t] <-0
    tps[5,5,t] <-0
    tps[5,6,t] <-if(t<21) {pt_dialysistotransplant[i]} else {0}
    tps[5,7,t] <-if (deathother[startage+t,2] > pt_dialysistodeath[i]) {deathother[startage+t,2]} else {pt_dialysistodeath[i]}
    #tps[5,8,t] <-0
    #tps[5,9,t] <-0
    #tps[5,10,t] <-0
    tps[5,5,t] <-1-sum(tps[5,,t])         


#Transitions from state 6=  ESRD Kidney Transplant [Remain with transplant,switch to dialysis,die] 
    tps[6,1,t] <-0
    tps[6,2,t] <-0
    tps[6,3,t] <-0
    tps[6,4,t] <-0
    tps[6,5,t] <-pt_transplanttodialysis[i]
    tps[6,6,t] <-0 #see end
    tps[6,7,t] <-if(deathother[startage+t,2] > pt_transplanttodeath[i]) {deathother[startage+t,2]} else {pt_transplanttodeath[i]}
    #tps[6,8,t] <-0
    #tps[6,9,t] <-0
    #tps[6,10,t] <-0
    tps[6,6,t] <-1-sum(tps[6,,t])

#Transitions from state 7=  Death [End State]
    tps[7,1,t] <-0
    tps[7,2,t] <-0
    tps[7,3,t] <-0
    tps[7,4,t] <-0
    tps[7,5,t] <-0
    tps[7,6,t] <-0
    tps[7,7,t] <-1
    #tps[7,8,t] <-0
    #tps[7,9,t] <-0
    #tps[7,10,t] <-0
 
    
#Transitions from state 8=  Spare
    #tps[8,1,t] <-0
    #tps[8,2,t] <-0
    #tps[8,3,t] <-0
    #tps[8,4,t] <-0
    #tps[8,5,t] <-0
    #tps[8,6,t] <-0
    #tps[8,7,t] <-0
    #tps[8,8,t] <-0
    #tps[8,9,t] <-0
    #tps[8,10,t] <-0
 
#Transitions from state 9=  Spare
    #tps[9,1,t] <-0
    #tps[9,2,t] <-0
    #tps[9,3,t] <-0
    #tps[9,4,t] <-0
    #tps[9,5,t] <-0
    #tps[9,6,t] <-0
    #tps[9,7,t] <-0
    #tps[9,8,t] <-0
    #tps[9,9,t] <-0
    #tps[9,10,t] <-0
 
#Transitions from state 10= Spare
    #tps[10,1,t] <-0
    #tps[10,2,t] <-0
    #tps[10,3,t] <-0
    #tps[10,4,t] <-0
    #tps[10,5,t] <-0
    #tps[10,6,t] <-0
    #tps[10,7,t] <-0
    #tps[10,8,t] <-0
    #tps[10,9,t] <-0
    #tps[10,10,t] <-0



#Error test
#Produces list of any states for which tps does not sum to 1, has negative or NA values
 for (s in 1:S){
     if ((round(sum(tps[s,,t]))==1)==FALSE)   {print(c(1,s,t))} #& print("tps for state does not sum to 1")
     if (any(tps[s,,t]< 0) == TRUE)           {print(c(2,s,t))} #& print("Value in tps for state <0") 
     if (any(is.na(tps[s,,t]))==TRUE)         {print(c(3,s,t))} #& print("NA in tps for state")
   }
}
    

tps <<- tps*100


#==================== MARKOV trace ========================================================================================================================

# Define what proportion of the patient cohort start in each health state 

# trace for t = 1  
# Patients health state at point they meet NICE referral critera based on mGFR

trace [1,1,i]   <-  p_3a_start[i]              
trace [1,2,i]   <-  p_3b_start[i]
trace [1,3,i]   <-  p_4_start[i]
trace [1,4,i]   <-  0
trace [1,5,i]   <-  0 
trace [1,6,i]   <-  0
trace [1,7,i]   <-  0  
trace [1,8,i]   <-  0
trace [1,9,i]   <-  0
trace [1,10,i]  <-  0  


# Define movement of cohort after first cycle according to transition probability matrix 

# trace for t >= 2

for (t in 2:T) {
  trace[t,,i] <- trace[t-1,,i] %*% tps[,,t]                     
}


# Store trace
trace<<-trace      

# test.NA  <- NA
# for (i in 1:S) {
#   test.NA[i]  <- any(is.na(trace[,i,]))
# }
# trace_NAs <- test.NA
# if (any(trace_NAs==TRUE)) print(trace_NAs) & stop("NA in trace") 


#======================= Outputs ============================================================================================================================

#========= QALYs ==========#

# Calculate QALYs: trace*utility

for (t in 1:T) {
  
  qtime[t] <-   (   trace[t,1,i]  * u_stage3a[i]
                  + trace[t,2,i]  * u_stage3b[i]
                  + trace[t,3,i]  * u_stage4[i]
                  + trace[t,4,i]  * u_stage5[i]
                  + trace[t,5,i]  * u_dialysis[i]
                  + trace[t,6,i]  * u_kidneytransplant[i]
                  + trace[t,7,i]  * 0 #death
                  + trace[t,8,i]  * 0
                  + trace[t,9,i]  * 0
                  + trace[t,10,i] * 0
                 )/((1+disc)^(t-1))    #Discount rate
}

#Total QALYs, applying half cycle correction
QALYs <- qtime[1]/2+sum(qtime[2:(T-1)])+qtime[T]/2               

# Store QALYs
qtime<<-qtime

#Test for errors
# test.NA  <- NA
# for (i in 1:T) {
#   test.NA[i]  <- any(is.na(qtime[i]))
# }
# qtime_NAs <- test.NA
# if (any(qtime_NAs==TRUE)) print(qtime_NAs) & stop("NA in model qtime vector") 


#========= COSTs ============#     

#Costs in first year & NOT to be included in half cycle correction
#Only includes those states which patients can occupy in the first year (#1,2,3)

cost_initial <- ( trace [1,1,i]  * (if (int=="referred") (cost_referral[i])  else    
                                                (0))    
                + trace [1,2,i]  * (if (int=="referred") (cost_referral[i])  else    
                                                (0))  
                + trace [1,3,i]  * (if (int=="referred") (cost_referral[i])  else    
                                                (0))   
                )
  

#Costs to be included in the half-cycle correction

for (t in 1:T) {
  
  cost[t] <-   (   trace [t,1,i]  * 10 #stage 3a - needs updating
                 + trace [t,2,i]  * 0 #stage 3b - needs updating
                 + trace [t,3,i]  * 0 #stage 4 - needs updating
                 + trace [t,4,i]  * 0 #stage 5 - needs updating
                 + trace [t,5,i]  * 0 #ESRD dialysis
                 + trace [t,6,i]  * 0 #ESRD transplant
                 + trace [t,7,i]  * 0 #death
                 + trace [t,8,i]  * 0
                 + trace [t,9,i]  * 0
                 + trace [t,10,i] * 0) /((1+disc)^(t-1))    #Discount rate
 }

#Total costs, applying half-cycle correction
cost <- cost_initial + cost[1]/2 + sum(cost[2:(T-1)])+cost[T]/2    


# Store costs
cost<<-cost 

# test.NA  <- NA
# for (i in 1:T) {
#   test.NA[i]  <- any(is.na(qtime[i]))
# }
# cost_NAs <- test.NA
# if (any(cost_NAs==TRUE)) print(cost_NAs) & stop("NA in model cost vector") 


#Model return:
return(c(QALYs,cost))

}
