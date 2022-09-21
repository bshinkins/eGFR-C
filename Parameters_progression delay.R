#Dialysis Model Parameters


parameter_draw  <- function() { 
  
  set.seed(seed)
  
  
  #====================== CODE FOR FUNCTIONS USED IN PARAMETERS ========================# 
  
  beta.draw <- function(n, m, s) {
    v <- s ^ 2
    rbeta(n, (m ^ 2 - m ^ 3 - m * v) / v,
          (m - 2 * m ^ 2 + m ^ 3 - v + m * v) /v)
  }
  
  gamma.draw <- function(n, m, s) {
    v <- s ^ 2
    rgamma(n, m ^ 2 / v, m / v)
  }
  
  lnorm.draw  <- function(n, m, s){    #CHECK THIS
    v <- s ^ 2
    rlnorm(n, log(m) - 0.5 * (log (1 + v / (m ^ 2))), 
           sqrt (log (1 + v / (m ^ 2))))
  }
  
  
  #============================== PATIENT CHARACTERISTICS ===================================#
  
  #Model cohort starting age           
  startage   <<-50  #needs amending once longer-term all cause mortality data added
  
  #Proportion of cohort male
  p_male     <- 0.5 #update with true proportion from study

  #=============================== STRUCTURAL PARAMETERS =====================================#
    
  H  <-50               #110-startage     # Time Horizon = 110 years - startage 
  
  T  <- H + 1           # Add one to allow for half cycle correction

  
  #============================== STARTING PROBABILITIES ==================================#
    
  p_3a_start<-0.4 #to be updated with data from main study  
  p_3b_start<-0.4 #to be updated with data from main study    
  p_4_start<-0.2  #to be updated with data from main study    
  p_5_start_start<-0
  p_dialysis_start<-0
  p_transplant_start<-0
    
  #============================== TRANSITION PROBABILITIES ==================================#
    
  pt_3ato3b<-0.096 #Elbasha
  pt_3bto4<-0.137
  pt_4to5<-0.081
  pt_4toESRD<-0
  pt_5todialysis<-0.626
  pt_5totransplant<-0.009
  pt_dialysistotransplant<-0.019
  pt_transplanttodialysis<-0.046
  
  pt_3todeath<-0.041 #Sugrue syst rev
  pt_4todeath<-0.080 #Sugrue syst rev
  pt_5todeath<-0.108 #Elbasha
  pt_dialysistodeath<-0.177 #Sugrue syst rev
  pt_transplanttodeath<-0.053 #Sugrue syst rev
  
  #reduction in progression if referred
  HR_stage3to4<-0.8 #Orlando et al
  HR_stage4to5<-0.75 #Orlando et al
    
  #Hazard rate of risk of all-cause mortality at each stage
  HR_3atodeath<-1.2
  HR_3btodeath<-1.8
  HR_4todeath<-3.2
  HR_5todeath<-5.9
  HR_dialysistodeath<-0.167
  HR_transplanttodeath<-0.028
  
  #Post-dialysis transplant rates by age
  #transplant_less50<-c(0.0459,0.05577,0.05616,0.05238,0.05081,0.04209,0.03312,0.02274,0.02449,0.01782)
  #transplant_50to64<-c(0.01731,0.02256,0.02664,0.02788,0.02832,0.02367,0.01903,0.01248,0.0135,0.005)
  #transplant_65to74<-c(0.0056,0.0074,0.0081,0.0076,0.006,0.0083,0.0036,0.0017,0.001,0)
  #transplant_75plus<-c(0.0005,0.0004,0.0004,0.0004,0,0,0,0,0,0)

  ##---TRANSITIONS TO DEATH---##
  
  #==== ALL-CAUSE MORTALITY ====#
  
  #All-cause mortality from: https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/adhocs/005944pastandprojectedmortalityratesqxfromthe2014basedlifetablesforengland1981to20800to125yrs
  #cohort statistics
  
  mort     <- read.csv("All cause mortality - to update.csv", stringsAsFactors = FALSE)
  mort     <- mort[-c(1,2),-c(2,4)]
  names(mort)[1:3] <- c("Age", "Male" , "Female")
  mort[,1] <- as.numeric(mort[,1])
  mort[,2] <- as.numeric(mort[,2])
  mort[,3] <- as.numeric(mort[,3])
  
  mort$pop    <- mort$Male*p_male + mort$Female*(1-p_male)
  deathother  <<- data.frame(mort$Age, mort$pop)    #In model ckd states, need to subtract proportion dying from ckd in each cycle from deathother (will need if statement so that subtraction does not lead to negative numbers)
  deathother  <<-deathother[-1,]
  
  #Hazard ratio of increased risk of all-cause mortality for each ckd stage
  HR_3a_allcause<-1.2
  HR_3b_allcause<-1.8
  HR_4_allcause<-3.2
  HR_5_allcause<-5.9

  }

#
  
  