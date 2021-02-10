##########################################################################################################
####################                        WB Model + Group Dispersal caused by hunting           ##################

#Input data and WBModel adapted to France - 9km sq
#Change the probability of movement by hunting manually
#setwd("X:/Luis/Model")
source("Initialization.R")
source("AnimalProcesses.R")
library(dplyr)
library(flux)
#library(doParallel)

rpert <- function(n, a, l, b) {
  mu <- (a+4*l+b)/6
  if (mu==l) v <- w <- 3 else {
    v <- (mu-a)*(2*l-a-b)/(l-mu)/(b-a)
    w <- v*(b-mu)/(mu-a)
  }
  a+(b-a)*rbeta(n, v, w)
}

WBModel <- function( MaxIterations    = 10,
                     MaxDays          = 365*5,
                     HuntingSeason    = (9*7+365)-(40*7),             # Hunting season goes from October to February
                     Detailed         = TRUE,
                     InitialOccCells  = 500,
                     MaxMortProbAd    = 1 - ((0.4)^(1/365)),          # Minimum Survival Prob for adults Lange et al. (2015), we recalculate per day 
                     MaxMortProbSAd   = 1 - ((0.4)^(1/365)),          # Minimum Survival Prob for Sub-adults Lange et al. (2015) we recalculate per day 
                     MaxMortProbPig   = 1 - ((0.1)^(1/365)),          # Minimum Survival Prob for piglets Lange et al. (2015) we recalculate per day
                     SurvivalProbAdF  = 0.9986^365,                   # Yearly Survival Prob for adults Kueling et al. (2013)    
                     SurvivalProbAdM  = 0.9986^365,                   # Yearly Survival Prob for adults Kueling et al. (2013) 
                     SurvivalProbSAdF = 0.9988^365,                   # Yearly Survival Prob for sub-adults (we use adults values as they are closer to Lange et al. 2015) Kueling et al. (2013)
                     SurvivalProbSAdM = 0.9986^365,                   # Yearly Survival Prob for sub-adults(we use adults values as they are closer to Lange et al. 2015) Kueling et al. (2013)    
                     SurvivalProbPigF = 0.50,                         # Yearly Survival Prob for piglets Lange et al. (2015) 
                     SurvivalProbPigM = 0.50,                         # Yearly Survival Prob for piglets Lange et al. (2015), reduced to consider higher survival of females
                     ProbHarvest      = 0.40,                         # A wild boar has 40% chance to be harvested
                     ProbHarvestAM    = 0.70,                         # An adult male has 70% of being harvested Toigo (2008)
                     ProbMovHunt      = 0.30,                         # Probability to move due to hunting
                     AgeProbab        = c(0.39, 0.24, 0.15, 0.09, 0.06, 0.03, 0.02, 0.01, 0.01), # Age probabilities from Lange et al. (2015)
                     ProbSplitSDS     = 1/7,  #cumsum(rep(1/7,7))     # Probability of short distance split if all conditions are satisfied for a group to split per day during the week of splitting (female splitting)
                     ProbSplitMSA     = 1/(6*7), #cumsum(rep(1/(6*7),6*7)) # Probability of male splitting per group per day of the week (7days) of the total splitting period (6 weeks)
                     SettlingProb     = c(0.75, 1),                   # probability to settle for a male in a habitat cell if it was only good habitat (0.5) or good habitat with females (1.0)
                     ProbMoveHS       = 1/HuntingSeason,              # daily probability of group splitting during the whole Hunting season
                     FemPropGroup     = mean(c(14.7,18.1,12.8,14.6)), # Proportion of females within a group Merta et al. (2015)  
                     MalPropGroup     = mean(c(9.2,9,5.7,7.9)),       # Proportion of males within a group Merta et al. (2015)
                     SubAPropGroup    = mean(c(33.4,34.3,22.8,24.1)), # Proportion of sub-adults within a group Merta et al. (2015)
                     PigletPropGroup  = mean(c(42.7,38.6,58.7,53.4)), # Proportion of piglets within a group Merta et al. (2015)
                     
                     ### Distance threshold for dispersion in meters
                     DistThreshold = 9000,
                     
                     ### ASF related parameters
                     WGProbInf        = 1 - (1 - 0.05)^(1/7),         # Direct contact Lange et al. (2015EFSA) and Lange and Thulke (2017)
                     CarcProbInf      = 1 - (1 - 0.2) ^(1/7),         # Carcass Lange et al. (2015EFSA) and Lange and Thulke (2017)
                     ProbAccesCarc    = 0.9,                          # Lange and Thulke (2017)
                     ProbDeathInf     = 0.95,                         # Blome (2012)  
                     CarcDaySur       = 28,                           # Lange and Thulke (2017)
                     DieInConEdge     = 0.8,                          # Lange et al. (2015EFSA)
                     TimeSeedInf      = 730,                          # Day 730 (2nd year)
                     NumbGSeedInf     = 1,                            # seed in one individual randomly
                     MinBCap          = 1,
                     ModeBCap         = 3,
                     MaxBCap          = 6,
                     AREA             = 1,
                     
                     #Reproduction Parameters
                     # RepProbList      = list(Week = 1:26,
                     #                         Prob = c(0, 0, 0, 0.01, 0.05, 0.1, 0.2, 0.23, 0.28, 0.325, 0.35, 0.325, 0.3, 0.28, 0.25, 0.22, 0.18, 0.17, 0.16, 0.15,
                     #                                0.14, 0.13, 0.12, 0.08, 0.04, 0.0)),
                     NumOfSprProbList = list(Number = 0:8,
                                             Prob = c(0.01, 0.07, 0.16, 0.25, 0.25, 0.16, 0.07, 0.02, 0.001)),
                     DirectionList    =  list(c(6, 7, 8), c(8, 10, 13), c(11, 12, 13), c(6, 9, 11)), ## North, ## East, ## South, and ## West
                     
                     HabitatProb      = c(0.75, 0.20, 0.05),         # The probability to select a new pixel while male walking (we model that they prefer to walk in a habitat cell - Hab Cat c(1,2,3)
                     FileWildBoarMat = "Inputs/DEPA64/Input_dep64.csv",
                     runID            = paste0("WB_Model", "Dep64")
){
  
  

# Raster and Distance Definition ------------------------------------------
  
  WBMat    <- DefineWbMat(FileWildBoarMat)
  coords   <- as.matrix(WBMat[ ,c('Lon', 'Lat')])
  Distance <- as.matrix(dist(coords))
  
  #Gamma distribution in reproduction parameters
  probRepro <- 0
  for(i in 1:51) probRepro <- c(probRepro, trapz(i:(i+1), dgamma(i:(i+1), shape = 6.4,scale = 2.1)))
  probRepro[1] = 1-sum(probRepro)
  RepProbList      = list(Week = 1:52,
                          Prob = probRepro)
  
  

# Initialize Model Outcomes.  ---------------------------------------------

  DayOutInfMat  <- matrix(0, ncol = MaxIterations, nrow = MaxDays)
  DayOutPopMat  <- matrix(0, ncol = MaxIterations, nrow = MaxDays)
  DayOutAniMat  <- matrix(0, ncol = MaxIterations, nrow = MaxDays)
  NewInfGroups  <- matrix(numeric(0), ncol = 4)
  NewInfAnimals <- matrix(numeric(0), ncol = 3)
  NewInfCarcass <- matrix(numeric(0), ncol = 3)
  EpDuration    <- rep(0, MaxIterations)
  FreqRelapse   <- rep(0, MaxIterations)
  cumDeath      <- rep(0, MaxIterations)
  Fpopcount     <- matrix(0, nrow = MaxDays, ncol = nrow(WBMat))
  Finfcount     <- matrix(0, nrow = MaxDays, ncol = nrow(WBMat))
  Fcarcount     <- matrix(0, nrow = MaxDays, ncol = nrow(WBMat))
  Fimmunecount  <- matrix(0, nrow = MaxDays, ncol = nrow(WBMat))
  

#  Initialize Mortality probabilities -------------------------------------
 
  #ProbMortFact deleted and hunting probability added
  ProbMortAdF  <- 1 - (SurvivalProbAdF)^(1/365)
  ProbMortAdM  <- 1 - (SurvivalProbAdM)^(1/365)
  ProbMortSAdF <- 1 - (SurvivalProbSAdF)^(1/365)
  ProbMortSAdM <- 1 - (SurvivalProbSAdM)^(1/365)
  ProbMortPigF <- 1 - (SurvivalProbPigF)^(1/365)
  ProbMortPigM <- 1 - (SurvivalProbPigM)^(1/365)
  ProbHunted   <- 1 - (1 - ProbHarvest)^(1/HuntingSeason)
  ProbHuntedAM <- 1 - (1 - ProbHarvestAM)^(1/HuntingSeason)  
  

# Start Iterations  -------------------------------------------------------

  for(Iter in 1:MaxIterations){
    
    set.seed(Iter)

    ## Add the habitat capacity
    

# Initialize the Population (Matrix of individuals) -----------------------
    PopMatWB <- InitPop(InitialOccCells, WBMat)

    
    GroupsToSplit   <- matrix(numeric(0), ncol = 4)     
    ## Store the groups that have split this year
    SplittedGroups  <- numeric(0)
    MovedGroups     <- numeric(0)
    cumDeathPar     <- 0
    cumDeathInf     <- 0
    Year            <- 1
    gTime           <- 0
    Criteria        <- TRUE
    TMPOutYesterday <- FALSE
    OnlyOnce        <- TRUE
    

# Start of daily loop -----------------------------------------------------
    while(Criteria & (gTime < MaxDays)){

      gTime <- gTime + 1
      if(gTime > 365) Year <- ceiling(gTime/365)

# Ageing ------------------------------------------------------------------
      PopMatWB <- Ageing(gTime, PopMatWB)

#  initiate variables and matrix that have to be initiated on year --------

      
      if(gTime %in% c(1, (366 * 1:(MaxDays/365)))){

# Reproduction ------------------------------------------------------------
  # Select females that will reproduce during the year
        PopMatWB        <- PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
        AnimalsWInG     <- unlist(sapply(unique(PopMatWB[ ,2]), function(x) 1:sum(PopMatWB[ ,2] == x)))
        BreedCapCells   <- c(WBMat[PopMatWB[ ,7], 15], rep(0, sum(PopMatWB[ ,7] == 0)))
        # Animals are allowed to breed when they are females and the group size is less than the maximum allowed size (breeding capacity)
        #(size that can be maintained by the pixel habitat)
        # and the animals above the breeding capacity based on age.
        AniAllowBreed   <- AnimalsWInG <= BreedCapCells & PopMatWB[ ,3] == 1 & PopMatWB[ ,9] < 2
        PopMatWB[AniAllowBreed, 6] <- (sample(RepProbList$Week, sum(AniAllowBreed), rep = T, prob = RepProbList$Prob)*7) + ((Year-1)*365)
       
        
        ### Reset Male Splitting
        
        
        #Males keep moving
        
        #PopMatWB[ ,12]          <- 0 #!PopMatWB[,2]%in%GroupsToSplit
        #IndexRem               <- which(PopMatWB[ ,2]%in%GroupsToSplit[ ,1])
        #PopMatWB[IndexRem, ]    <- 0
        #PopMatWB               <- PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3],PopMatWB[ ,5], decreasing = T), ]
        #cumDeathPar            <- cumDeathPar + length(IndexRem)
        #GroupsToSplit          <- matrix(numeric(0), ncol = 4) 
      }
      

# Mortality ---------------------------------------------------------------


  ToDieNormal <- Mortality(PopMatWB, ProbMortPigF, ProbMortPigM,
                       ProbMortSAdF, ProbMortSAdM, ProbMortAdF, ProbMortAdM)
      #Mortality caused by Hunting in population and male adults
      #Take out Piglets from being hunted
     
      #to check if there are no piglets and male adults
      #table(data.frame(PopMatWB[which(!(PopMatWB[ ,3] ==  0 & PopMatWB[ ,4] == 3)  &  PopMatWB[ ,4] != 1),c("Sex", "Age_Cat")]))
      
      
# Hunting Season ----------------------------------------------------------
      
      if(gTime %in% ((40*7 + ((Year - 1)*365)):((9*7 + ((Year)*365))))){
        
        #ToDieHunted <- which(rbinom(nrow(PopMatWB), 1, ProbHunted) == 1) #Without Hunting Pro for AM
        ToDieHunted   <- which(!(PopMatWB[ ,3] ==  0 & PopMatWB[ ,4] == 3)  &  PopMatWB[ ,4] != 1)
        ToDieHuntedAM <- which(PopMatWB[ ,3] == 0 & PopMatWB[ ,4] == 3)
        ## Hunting without piglets and with a different hunting probability for adult males
        ToDieHunted <- ToDieHunted[rbinom(length(ToDieHunted), 1, ProbHunted) == 1]
        ToDieHuntedAM <- ToDieHuntedAM[rbinom(length(ToDieHuntedAM), 1, ProbHuntedAM) == 1]
        Hunted <- c(ToDieHuntedAM, ToDieHunted)
        
        ## Movement caused by hunting during Hunting Season
        ## Groups that will split due to hunting
        GroupsMoveHunt <- sort(unique(PopMatWB[Hunted ,2]))
        #GroupsMoveHunt <- GroupsMoveHunt[GroupsMoveHunt[ ,2] >= 1, , drop = FALSE]
        GroupsMoveHunt <- GroupsMoveHunt[rbinom(length(GroupsMoveHunt), 1, ProbMovHunt) == 1]
        
        GroupsMoveShD  <- sapply(PopMatWB[match(GroupsMoveHunt, PopMatWB[ ,2]), 7], function(x){
          tempo1     <- which(Distance[x,] <= DistThreshold)
          any(WBMat[tempo1, 5] == 1)})
        
        GroupsMoveNumSD <- cbind(GroupsMoveHunt,  GroupsMoveShD)
        GroupsMoveNumSD <- GroupsMoveNumSD[GroupsMoveNumSD[ ,2] == 1, , drop = FALSE]
        
        ## Keep track of the groups that have moved due to hunting
        PopMatWB[PopMatWB[ ,2] %in% GroupsMoveNumSD[ ,1], 16] <- 2
        
        ## Make Short distance happened
        if(dim(GroupsMoveNumSD)[1] > 0){
          OriginalPixel     <- PopMatWB[match(GroupsMoveNumSD[ ,1], PopMatWB[ ,2]), 7]
          
          TargetPixel <- numeric(0)
          for(xx in OriginalPixel){
            tempo1     <- which(Distance[xx,] <= DistThreshold)
            tempo2      <- tempo1[WBMat[tempo1, 5] == 1 & !(tempo1 %in% PopMatWB[ ,7])]
            if(length(tempo2) > 1)  tempo3 <- sample(tempo2, 1)
            if(length(tempo2) == 1) tempo3 <- tempo2
            if(length(tempo2) == 0) tempo3 <- 0
            TargetPixel <- c(TargetPixel, tempo3)
          }
          
          GroupsMoveNumSD <- GroupsMoveNumSD[TargetPixel > 0, , drop = FALSE]
          OriginalPixel     <- OriginalPixel[TargetPixel > 0]
          TargetPixel     <- TargetPixel[TargetPixel > 0]
          GroupsMoveNumSD <- cbind(GroupsMoveNumSD, TargetPixel) 
          
          ## Allow Movement
          if(dim(GroupsMoveNumSD)[1] > 0){
            ## Create a list to keep track of pixels where the group have been. 
            ## Exported on a daily basis, the list will be re-initiate every day.
            ## How many times does the group tries to move?
            PixelsMovedHS <- as.list(matrix(0, ncol = length(TargetPixel)))
            CurrentPosition <- OriginalPixel
            
            for(i in 1:length(TargetPixel)){
              Trail <- 0
              previousEdge <- CurrentPosition[i]
              while(CurrentPosition[i] != TargetPixel[i] & Trail < 10){
                Trail      <- Trail + 1
                Edges      <- unlist(WBMat[CurrentPosition[i], 6:13])
                Edges      <- Edges[Edges > 0 & Edges != previousEdge]
                if(length(Edges) > 0){
                  previousEdge <- CurrentPosition[i]
                  NewPosition  <- Edges[which.min(Distance[TargetPixel[i],Edges])]
                  if(length(NewPosition) > 1) NewPosition <- sample(NewPosition, 1)
                  CurrentPosition[i] <- NewPosition
                  
                  # Here we keep track of the pixels where the group has moved to.
                  PixelsMovedHS[[i]] <- c(PixelsMovedHS[[i]], NewPosition)
                  #print(c(i,OriginPixel[i],prevEdge,Trail,CurrentPos[i],TargetPixel[i]))
                }
              }
            }
            
            names(PixelsMovedHS) <- GroupsMoveNumSD[ ,1]
            MovedGroupsHS <- PopMatWB[ ,2] %in% GroupsMoveNumSD[ ,1]
            #IndexNotSplit    <- which(CurrentPosition != TargetPixel)
            #GroupsMoveNumSD <- GroupsMoveNumSD[-IndexNotSplit, , drop = FALSE]
            
            ### Update Home Pixel
            
            #UPDATE IN MATHEU'S CODE
            for(b in 1:dim(GroupsMoveNumSD)[1]){
              GroupsCanMove <-  PopMatWB[ ,2] %in% GroupsMoveNumSD[b, 1]
              NewHomePixel     <- rep(GroupsMoveNumSD[b, 3], sum(GroupsCanMove))
              PopMatWB[GroupsCanMove, 7]  <- NewHomePixel
              PopMatWB[GroupsCanMove, 8]  <- NewHomePixel
              
            }
          }
        }
        
        
        ## Include hunted individuals in the death toll. Else is outside Hunting Season
        
        ToDieAll <- unique(c(ToDieNormal, ToDieHunted, ToDieHuntedAM)) 
      } else {
        ToDieAll <- ToDieNormal
      }

      ## Include death toll in the cumulative deaths (not infected)
      
      if(length(ToDieAll) > 0){
        PopMatWB <- PopMatWB[-ToDieAll, ]
        cumDeathPar <- cumDeathPar + length(ToDieAll)
      }
      
      ## Reset the movement due to hunting each year
      
      if(gTime %in% ((40*7 + ((Year - 1)*365)):((9*7 + ((Year)*365))) + 1)){
        PopMatWB[ ,16]    <- 0
        #SplittedGroups   <- numeric(0)
        PixelsMovedHS      <- as.list(matrix(0, ncol=1))
      }
      
      

# Reproduction ------------------------------------------------------------

      
      ## On daily basis, from Jan to end June, check if there are animals to deliver and make them deliver
        if(gTime %in% PopMatWB[ ,6]){
        DelIndex        <- which(PopMatWB[ ,6] %in% gTime & PopMatWB[ ,9] != 3)
        NumOfSpring     <- sample(NumOfSprProbList$Number, length(DelIndex), rep = T, prob = NumOfSprProbList$Prob)
        DelIndex        <- DelIndex[NumOfSpring > 0]
        NumOfSpring     <- NumOfSpring[NumOfSpring > 0]
        if(length(NumOfSpring) > 0){
          newBorns <- cbind((max(PopMatWB[ ,1]) + 1):(max(PopMatWB[ ,1]) + sum(NumOfSpring)),
                          rep(PopMatWB[DelIndex, 2], NumOfSpring),
                          rbinom(sum(NumOfSpring), 1, 0.5),
                          1,
                          1,
                          0,
                          rep(PopMatWB[DelIndex, 7], NumOfSpring),
                          rep(PopMatWB[DelIndex, 7], NumOfSpring),
                          0,
                          0,
                          rep(PopMatWB[DelIndex, 1], NumOfSpring),
                          0,
                          0,
                          0,
                          0,
                          0)
          colnames(newBorns) = colnames(PopMatWB)
          PopMatWB <- rbind(PopMatWB, newBorns)         
          PopMatWB <- PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
          
        }
      }


# Female Group Splitting --------------------------------------------------

      
      ## Occurs only once in week 28 Kramer-Schadt et al. (2009)
      if(gTime %in% ((28*7 + ((Year - 1)*365)) + 0:6)){
        
        ## Calculate the day in the week
        #DayInWeekSSplit <- (1-(ceiling((gTime-(365*(Year-1)))/7) - ((gTime-(365*(Year-1)))/7)))*7
        PopMatWB        <- PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
        AnimalsWInG     <- unlist(sapply(unique(PopMatWB[ ,2]), function(x) 1:sum(PopMatWB[ ,2] == x)))
        BreedCapCells   <- c(WBMat[PopMatWB[ ,7] ,15], rep(0, sum(PopMatWB[ ,7] == 0)))
        ## Identify the animals within the groups that will split
        ## Females sub-adults > breading capacity, did not split before and are not sick
        GroupSplitNum   <- cbind(sort(unique(PopMatWB[ ,2])), (tapply((AnimalsWInG > BreedCapCells &
                                                                        PopMatWB[ ,4] == 2 &
                                                                        PopMatWB[ ,3] == 1 &
                                                                        PopMatWB[,10] == 0 &
                                                                        PopMatWB[,9] < 2),  PopMatWB[ ,2], sum)))
        GroupSplitNum   <- GroupSplitNum[GroupSplitNum[ ,2] >= 2, , drop = FALSE]
        
        ## Identify groups where there are only males
        ## To allow females to join groups where there are only males.
        PixelMalesOnly  <- cbind(sort(unique(PopMatWB[ ,2])), (tapply((PopMatWB[ ,3] == 0), PopMatWB[ ,2], all)))
        PixelMalesOnly  <- PixelMalesOnly[PixelMalesOnly[ ,2] == 1, , drop = FALSE]
        
        GroupsSplitShD  <- sapply(PopMatWB[match(GroupSplitNum[ ,1], PopMatWB[ ,2]), 7], function(x) {
          tmp1 <-which(Distance[x,] <= DistThreshold)
          any(WBMat[tmp1, 5] == 1 & !(tmp1 %in% PopMatWB[!PopMatWB[ ,2] %in% PixelMalesOnly[ ,1], 7]))})
        GroupSplitNumSD <- cbind(GroupSplitNum, GroupsSplitShD)
        GroupSplitNumSD <- GroupSplitNumSD[GroupSplitNumSD[ ,3] == 1, , drop=FALSE]
        
        ## Groups to split, may split during any day of the week. No need for a probability here the group splits only once a year.
        #if(dim(GroupSplitNumSD)[1]>0) GroupSplitNumSD <- GroupSplitNumSD[rbinom(dim(GroupSplitNumSD)[1],1,prob=ProbSplitSDS)==1,,drop=FALSE]
        PopMatWB[PopMatWB[ ,2] %in% GroupSplitNumSD[ ,1], 10] <- 1  #Make all group split status 1??
        
        ## Make short term split to happen
        if(dim(GroupSplitNumSD)[1] > 0){
          OriginPixel     <- PopMatWB[match(GroupSplitNumSD[ ,1], PopMatWB[ ,2]), 7]
          
          TargetPixel <- numeric(0)
          for(xx in OriginPixel){
            tmp1 <-which(Distance[xx,] <= DistThreshold)
            tmp2      <- tmp1[WBMat[tmp1, 5] == 1 & !(tmp1 %in% PopMatWB[ ,7])]
            if(length(tmp2) > 1)  tmp3 <- sample(tmp2, 1)
            if(length(tmp2) == 1) tmp3 <- tmp2
            if(length(tmp2) == 0) tmp3 <- 0
            TargetPixel <- c(TargetPixel, tmp3)
          }
          
          GroupSplitNumSD <- GroupSplitNumSD[TargetPixel>0, , drop = FALSE]
          OriginPixel     <- OriginPixel[TargetPixel > 0]
          TargetPixel     <- TargetPixel[TargetPixel > 0]
          GroupSplitNumSD <- cbind(GroupSplitNumSD, TargetPixel) 
          
          ## This part of the code to allow movement
          if(dim(GroupSplitNumSD)[1] > 0){
            
            ## A list to keep track for the pixels where the pigs have been. 
            ## Make sure that this information is exported on daily basis, because the 
            ## list will be re-initiate every day splitting may happen.
            PixelsMoved <- as.list(matrix(0, ncol = length(TargetPixel)))
            
            
            CurrentPos <- OriginPixel
            for(i in 1:length(TargetPixel)){
              Trail <- 0
              prevEdge <- CurrentPos[i]
              while(CurrentPos[i] != TargetPixel[i] & Trail < 10){
                Trail      <- Trail + 1
                Edges      <- unlist(WBMat[CurrentPos[i], 6:13])
                Edges      <- Edges[Edges > 0 & Edges != prevEdge]
                if(length(Edges) > 0){
                  prevEdge    <- CurrentPos[i]
                  NewPosition <- Edges[which.min(Distance[TargetPixel[i],Edges])]
                  if(length(NewPosition) > 1) NewPosition <- sample(NewPosition, 1)
                  CurrentPos[i] <- NewPosition
                  # Here we keep track of the pixels where the pigs moved to.
                  PixelsMoved[[i]]<- c(PixelsMoved[[i]], NewPosition)
                  #print(c(i,OriginPixel[i],prevEdge,Trail,CurrentPos[i],TargetPixel[i]))
                }
              }
            }
            
            ## Here we assume that groups that did not find the way, did not actually split as the edges were not connected.
            IndexNotSplit   <- which(CurrentPos != TargetPixel)
            #GroupSplitNumSD <- GroupSplitNumSD[-IndexNotSplit, ,drop = FALSE]
            
            ## Make females from different groups moving to a new pixel to form a new group.
            if(dim(GroupSplitNumSD)[1] > 0){
              newGroupIDs <- (max(PopMatWB[,2]) + 1):(max(PopMatWB[ ,2]) + dim(GroupSplitNumSD)[1])
              if(sum(duplicated(GroupSplitNumSD[ ,4])) > 0){ #Target pixel
                for(l in unique(GroupSplitNumSD[ ,4])){
                  TEMP  <- which(GroupSplitNumSD[ ,4] == l)
                  if(length(TEMP) > 1){
                    newGroupIDs[TEMP] <- newGroupIDs[TEMP[1]]
                  }
                }
              }
              
              names(PixelsMoved) <- GroupSplitNumSD[ ,1]
              for(b in 1:dim(GroupSplitNumSD)[1]){
                
                ## If there is a male group already in the pixel and is not splitting, then make them one group with the new females.
                ## Notice the code above allows them to come into a new pixel only if it is empty or there are only males. 
                ## Dont worry tat IndexMalComb1 is matching with all PopMatWB[,7]
                IndexMalComb1  <- GroupSplitNumSD[b, 4] %in% PopMatWB[ ,7]
                if(IndexMalComb1){
                  IndexMalComb2  <- unique(PopMatWB[PopMatWB[ ,7] == GroupSplitNumSD[b, 4], 2]) %in% GroupsToSplit[ ,1]
                  if(!IndexMalComb2){
                    AnimalsCanSplitF <- AnimalsWInG>BreedCapCells & 
                                          PopMatWB[ ,3] == 1 & 
                                          PopMatWB[ ,4] == 2 & 
                                          PopMatWB[ ,9] < 2  & 
                                          PopMatWB[ ,2] %in% GroupSplitNumSD[b, 1]
                    MalesInPixel     <- PopMatWB[ ,3] == 0 & 
                                        PopMatWB[ ,9] != 3 & 
                                        PopMatWB[ ,7] == GroupSplitNumSD[b, 4] & 
                                       !PopMatWB[ ,2] %in% GroupsToSplit[ ,1]
                    newGroupIDsAn    <- rep(newGroupIDs[b], sum(AnimalsCanSplitF) + sum(MalesInPixel)) 
                    NewHomePixel     <- rep(GroupSplitNumSD[b, 4], sum(AnimalsCanSplitF) + sum(MalesInPixel))
                    PopMatWB[AnimalsCanSplitF|MalesInPixel, 2]  <- newGroupIDsAn
                    PopMatWB[AnimalsCanSplitF|MalesInPixel, 7]  <- NewHomePixel
                    PopMatWB[AnimalsCanSplitF|MalesInPixel, 8]  <- NewHomePixel
                  }
                  if(IndexMalComb2){
                    AnimalsCanSplitF <- AnimalsWInG > BreedCapCells & 
                                        PopMatWB[ ,3] == 1 & 
                                        PopMatWB[ ,4] == 2 & 
                                        PopMatWB[ ,9] < 2  & 
                                        PopMatWB[ ,2] %in% GroupSplitNumSD[b, 1]
                    newGroupIDsAn    <- rep(newGroupIDs[b], sum(AnimalsCanSplitF)) 
                    NewHomePixel     <- rep(GroupSplitNumSD[b, 4], sum(AnimalsCanSplitF))
                    PopMatWB[AnimalsCanSplitF, 2]  <- newGroupIDsAn
                    PopMatWB[AnimalsCanSplitF, 7]  <- NewHomePixel
                    PopMatWB[AnimalsCanSplitF, 8]  <- NewHomePixel
                  }   
                }
                if(!IndexMalComb1){
                  AnimalsCanSplitF <- AnimalsWInG > BreedCapCells & 
                                      PopMatWB[ ,3] == 1 & 
                                      PopMatWB[ ,4] == 2 & 
                                      PopMatWB[ ,9] < 2  & 
                                      PopMatWB[ ,2] %in% GroupSplitNumSD[b, 1]
                  newGroupIDsAn    <- rep(newGroupIDs[b], sum(AnimalsCanSplitF)) 
                  NewHomePixel     <- rep(GroupSplitNumSD[b, 4], sum(AnimalsCanSplitF))
                  PopMatWB[AnimalsCanSplitF, 2]  <- newGroupIDsAn
                  PopMatWB[AnimalsCanSplitF, 7]  <- NewHomePixel
                  PopMatWB[AnimalsCanSplitF, 8]  <- NewHomePixel
                }   
              }
            }
          } 
        }
      }#Closes short distance splitting
      
      ## After the end of female splitting, we reset female splitting.
      if(gTime %in% ((28*7 + ((Year - 1)*365)) + 7)){
        PopMatWB[ ,10]    <- 0
        SplittedGroups   <- numeric(0)
        PixelsMoved      <- as.list(matrix(0, ncol = 1))
      }
      
      
      ### MALE GROUP SPLITTING
     
      
      ## First we determine the period of splitting of males, as defined in Lange et al., (2012). 
      ## Males may find a pixel to live in otherwise they will keep Wondering around 
      ## until they either die or find a place to live in
      if(gTime%in%((25*7 + ((Year - 1)*365)):(30*7 + ((Year - 1)*365)))){
        
        ## Sort the matrix
        PopMatWB    <-  PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
        ## Define the males as adults and sub-adults
        ## Identify the groups where sub-adult males may split
        ## Male sub-adults from groups where males have not yet split this year.
        #GroupSplittedThisYear <- cbind(sort(unique(PopMatWB[,2])),tapply(PopMatWB[,12],PopMatWB[,2],sum))
        #GroupSplittedThisYear <- GroupSplittedThisYear[GroupSplittedThisYear[,2]>0,,drop=FALSE]
        
        ## Calculate number of adult males
        NumbAdultMInG <-  cbind(sort(unique(PopMatWB[ ,2])), tapply(PopMatWB[ ,3] == 0 & 
                                                                    PopMatWB[ ,4] == 3, 
                                                                    PopMatWB[ ,2], sum))
        ## No need for sub-adult males to split if there is need for them in the group
        tmpIndexSG    <- NumbAdultMInG[match(PopMatWB[ ,2], NumbAdultMInG[ ,1]), 2]
        ## Groups that may split must not be group already splitting and have more than 2 adult males
        GroupSubAM    <- unique(PopMatWB[!(PopMatWB[ ,2] %in% SplittedGroups) & 
                                           PopMatWB[ ,3]  == 0 & 
                                           PopMatWB[ ,4]  == 2 & 
                                           PopMatWB[ ,12] == 0 & 
                                           tmpIndexSG > 2, 2])
        if(length(GroupSubAM) > 0) 
          GroupSubAM <- GroupSubAM[rbinom(length(GroupSubAM), 1, ProbSplitMSA) == 1]
        if(length(GroupSubAM) > 0){
         
        ## Identify the sub-adult males that may split.
        ## We must ensure that the groups have adult males before the sub-adults are split.
        ## No need for the sub-adults to find new areas if they can nourish and breed in their home
          #TEMPMat <- cbind(sort(unique(PopMatWB[,2])),tapply((PopMatWB[,3]==0 & PopMatWB[,4]==3 & PopMatWB[,9] < 2),  PopMatWB[,2], sum))
          tmpSAM <- which(PopMatWB[ ,2] %in% GroupSubAM & 
                            PopMatWB[ ,3] == 0 & 
                            PopMatWB[ ,4] == 2 & 
                            PopMatWB[ ,9] < 2)#& TEMPMat[match(PopMatWB[,2],TEMPMat[,1]),2]>=2)
          ## Check which ones that would split
          if(length(tmpSAM) > 0) {
            ## Make subadults to split based on probability from 25 to 75% (Truve et al., 2004).
            tmpSAM           <- tmpSAM[rbinom(length(tmpSAM), 1, prob = runif(length(tmpSAM), 0.25, 0.75)) == 1]
            AinmFromGToSplit <- cbind(sort(unique(PopMatWB[tmpSAM, 2])), tapply(tmpSAM, PopMatWB[tmpSAM, 2], length))
            AinmFromGToSplit <- AinmFromGToSplit[AinmFromGToSplit[ ,2]>1, ,drop = FALSE]
            tmpSAM           <- tmpSAM[PopMatWB[tmpSAM ,2] %in% AinmFromGToSplit[ ,1]]
            
            if(length(tmpSAM) > 0) {# do not worry, if the dimension is > 0 then there will be at least 2 animals from one group to split ;-) check above.
              PopMatWB[tmpSAM,12] <- 1
              ## Store the group number of the groups that have splited
              SplittedGroups <- c(SplittedGroups, unique(PopMatWB[tmpSAM, 2]))
              ## Number of splitting males per group
              AinmFromGToSplit   <- cbind(AinmFromGToSplit, round(runif(dim(AinmFromGToSplit)[1], 1, 4)), unique(PopMatWB[tmpSAM, 8]))
              ## The third column, informs about the direction
              ## Which should be depicted from the DirectionList. each group would walk in a specific direction
              ## If they cannot find a connecting edge in that direction, a random cell is then selected.
              ## Make them to prefer to move through habitat cells rather than fields.
              indexsNGroups  <- (max(PopMatWB[ ,2]) + 1):(max(PopMatWB[ ,2]) + dim(AinmFromGToSplit)[1])
              indexNGPigs    <- rep(indexsNGroups, AinmFromGToSplit[ ,2])
              indexNewGNum   <- sort(which(PopMatWB[ ,2] %in% AinmFromGToSplit[ ,1] & PopMatWB[ ,12] == 1), decreasing = T)
              PopMatWB[indexNewGNum, 2] <- indexNGPigs
              AinmFromGToSplit[ ,1] <- indexsNGroups
              GroupsToSplit <- rbind(GroupsToSplit, AinmFromGToSplit)      
            }
          }
        }
      }
      
      if(dim(GroupsToSplit)[1] > 0){
        ## First we check that all groups still exist in the population
        Checktmp         <- which(!GroupsToSplit[ ,1] %in% PopMatWB[ ,2])
        if(length(Checktmp) > 0) GroupsToSplit[-Checktmp, , drop = FALSE]
        ## Determine the number of pixels they may move per day and the direction
        NumbPixMovesTod  <- round(runif(dim(GroupsToSplit)[1], 0, 1)) ## Lange et al. (2015)
        MovedPixMale <- as.list(matrix(0, ncol = length(NumbPixMovesTod > 0)))
        ToRemovejj <- numeric(0)
        if(any(NumbPixMovesTod > 0)){
          for(i in 1:max(NumbPixMovesTod)){
            if(any(PopMatWB[ ,12] == 1) & any(NumbPixMovesTod >= i)){
              IndexGroup   <- which(NumbPixMovesTod >= i)
              if(length(IndexGroup) > 0){
                for(jj in IndexGroup){
                  #Find the directions for all neighbors
                  NeighboursIndex <- unlist(WBMat[GroupsToSplit[jj, 4], DirectionList[[GroupsToSplit[jj, 3]]]])
                  #Keep The ones in the selected area :exclude 0's
                  NeighboursIndex <- NeighboursIndex[NeighboursIndex > 0]
                  #Check the habitat cat
                  NbHabCat <- WBMat[NeighboursIndex, 5]
                  #Keep accessible and suitable ones only
                  SuitableNb <- which(NbHabCat > 0 & NbHabCat <= length(HabitatProb))
                  TSuitableNb <- NeighboursIndex[SuitableNb]
                  
                  # if(length(TSuitableNb)>0){
                    #if only 1 suitable neighbour => go there
                    if(length(TSuitableNb) == 1) Destination <- TSuitableNb
                    #if more thqn one select one neighbour randomly
                    if(length(TSuitableNb) > 1) Destination <- sample(TSuitableNb, 1 , prob = HabitatProb[NbHabCat[SuitableNb]])
                  # }
                  if(length(TSuitableNb) == 0){
                    NeighboursIndex <- unlist(WBMat[GroupsToSplit[jj, 4], 6:13])
                    NeighboursIndex <- NeighboursIndex[NeighboursIndex > 0]
                    NbHabCat <- WBMat[NeighboursIndex, 5]
                    #SuitableNb <- which(NbHabCat>0 & NbHabCat<=3)
                    SuitableNb <- which(NbHabCat > 0 & NbHabCat <= length(HabitatProb))
                    TSuitableNb <- NeighboursIndex[SuitableNb]
                    if(length(TSuitableNb) > 1)  Destination <- sample(TSuitableNb, 1, prob = HabitatProb[NbHabCat[SuitableNb]])
                    if(length(TSuitableNb) == 1) Destination <- TSuitableNb
                  }
                  
                  ## Here we add the list with the moved pixels
                  MovedPixMale[[jj]] <- c(MovedPixMale[[jj]], Destination)
                  ## Is this only a suitable habitat pixel or a suitable habitat and has a female(s) with no or 1 male
                  
                  TMPIndPix <- (WBMat[Destination, 5] == 1 & 
                                sum(PopMatWB[PopMatWB[ ,7] == Destination & 
                                PopMatWB[ ,9] != 3, 3] == 1) > 0 & 
                                sum(PopMatWB[PopMatWB[ ,7] %in% Destination, 9] != 3) > 0 &
                                sum(PopMatWB[PopMatWB[ ,7] == Destination, 3] == 0) < 2) #& unique(PopMatWB[PopMatWB[,7]%in%tmp5,2])<2
                  # the first part of the c() represent a good habitat cell and not occupied
                  Freepixel      <- which(c(WBMat[Destination, 5] == 1 & !(Destination %in% PopMatWB[ ,7]), TMPIndPix))
                  
                  ## We decide whether they will settle in this pixel or not based on a random process
                  ## Notice here that the animals that die during splitting do not affect the splitting. 
                  ## We checked above that all groups in splitting matrix have animals in the Population matrix
                  if(length(Freepixel) > 0){
                    Settle <- rbinom(1, 1, prob = SettlingProb[max(Freepixel)]) == 1
                    if(Settle&!TMPIndPix){
                      indexspPigs             <- which(PopMatWB[ ,2] == GroupsToSplit[jj, 1] & PopMatWB[,12]==1 & PopMatWB[,9]!=3)
                      #newGroupIDM             <- rep(max(PopMatWB[,2])+1,length(indexspPigs))
                      #PopMatWB[indexspPigs,2] <- newGroupIDM
                      PopMatWB[indexspPigs,7] <- Destination
                      PopMatWB[indexspPigs,8] <- Destination
                      PopMatWB[indexspPigs,12]<- 0
                      ToRemovejj              <- c(ToRemovejj,jj)
                    }
                    
                    ## If the group will settle with a pre-existing female group, then their group number will be the same as the females
                    if(Settle&TMPIndPix){
                      indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1 & PopMatWB[,9]!=3)
                      IndexSelPigs            <- which(PopMatWB[,7]%in%Destination & PopMatWB[,9]!=3 & !PopMatWB[,2]%in%GroupsToSplit[jj,1] & !PopMatWB[,2]%in%GroupsToSplit[GroupsToSplit[,4]%in%PopMatWB[,7],1])
                      newGroupIDM             <- rep(unique(PopMatWB[IndexSelPigs,2]),length(indexspPigs)) #this would allow settling in an infected pixel with carcass 
                      PopMatWB[indexspPigs,2] <- newGroupIDM
                      PopMatWB[indexspPigs,7] <- Destination
                      PopMatWB[indexspPigs,8] <- Destination
                      PopMatWB[indexspPigs,12]<- 0
                      ToRemovejj              <- c(ToRemovejj,jj)
                    }
                    if(!Settle){
                      indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1 & PopMatWB[,9]!=3)
                      PopMatWB[indexspPigs,8] <- Destination
                      GroupsToSplit[jj,4]     <- Destination
                    }
                  }
                  
                  if(length(Freepixel) == 0) {
                    indexspPigs             <- which(PopMatWB[ ,2]  == GroupsToSplit[jj, 1] & 
                                                     PopMatWB[ ,12] == 1 & 
                                                     PopMatWB[ ,9]  != 3)
                    PopMatWB[indexspPigs, 7] <- Destination
                    PopMatWB[indexspPigs, 8] <- Destination ## we set the home pixel as the current pixel for splitting groups as they do not have home yet
                    GroupsToSplit[jj, 4]     <- Destination
                  }
                  
                }
              }
            }
          }
          names(MovedPixMale) <- GroupsToSplit[,1]
          if(length(ToRemovejj)>0) GroupsToSplit  <- GroupsToSplit[-ToRemovejj, ,drop = FALSE]
        }
      }
       
        
      ### SEED INFECTION
        

       if(TimeSeedInf == gTime){
       
       ## Here I seed randomly in 1 group and infect 1 animal (make it directly infectious). 
       ## We make sure that the groups have at least 2 animals.
               
       NumPGroup  <- tapply(PopMatWB[ ,2], PopMatWB[ ,2], length)
       temp1      <- NumPGroup[NumPGroup > 2]
       temp2      <- as.numeric(names(temp1))
       temp3      <- temp2[temp2 > 0]
       RandSeedG  <- sample(temp3, NumbGSeedInf)
       SampSeedAn <- sapply(RandSeedG, function(x) sample(which(PopMatWB[ ,2] == x), 1))
       PopMatWB[SampSeedAn, 9]  <- 2
       PopMatWB[SampSeedAn, 14] <- gTime   # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
       PopMatWB[SampSeedAn, 15] <- PopMatWB[SampSeedAn,14] + round(rpert(length(SampSeedAn), 1, 5, 7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
       
          }
  

# Carcass Persistence -----------------------------------------------------
  
  ## Winter
  if(gTime %in% ((49*7) + ((Year - 1)*365)):((9*7 + ((Year)*365)))){
    CarcDaySur <- 90
  } else
  ## Summer
  if(gTime %in% ((22*7 + (Year*365))):((35*7 + (Year*365)))){
    CarcDaySur <- 5
    }
  

# ASF Dynamics ------------------------------------------------------------
      
      ## A - Here we model within group infection and infection from connecting edges via carcasses
   
      if(sum(PopMatWB[ ,9]) > 0){
        infGroups <- unique(PopMatWB[PopMatWB[ ,9] %in% 2:3, 2])
        for(i in infGroups){
          NumInfPG  <- sum(PopMatWB[PopMatWB[,2] == i ,9] == 2)
          TotDead   <- sum(PopMatWB[PopMatWB[,2] == i, 9] == 3 & 
                             (gTime - PopMatWB[PopMatWB[ ,2] == i, 15]) <= CarcDaySur)
          ## Prob from infectious individuals
          probWGDCInf <- 1 - (1 - WGProbInf)^NumInfPG
          ## Prob from infectious carcasses within the group
          ProbWGCInf  <- 1 - (1 - CarcProbInf)^TotDead
          ## Prob from infectious carcasses in connecting edges (based on current pixel NOT Home pixel as they maybe splitting)
          tmpConEdgs  <- unlist(WBMat[unique(PopMatWB[PopMatWB[ ,2] == i, 8]), 6:13])
          tmpConEdgs  <- tmpConEdgs[tmpConEdgs > 0]
          tmpDInCE    <- sum(PopMatWB[ ,8] %in% tmpConEdgs & PopMatWB[ ,9] == 3)
          ProbBGCInf  <- 1 - ( 1 - CarcProbInf)^tmpDInCE
          ## Estimate total probability of infection from contact
          ProbWGInf   <- 1 - (( 1- probWGDCInf)*(1 - ProbWGCInf)*(1 - ProbBGCInf))
          
          ## Susceptible animals
          tmpSus <- which(PopMatWB[ ,9] == 0 & PopMatWB[ ,2] == i)
          if(length(tmpSus) > 0){
            NewInfPG <- tmpSus[rbinom(length(tmpSus), 1, prob = ProbWGInf) == 1]
            if(length(NewInfPG) > 0){
              PopMatWB[NewInfPG,  9] <- 1 ## exposed individuals
              PopMatWB[NewInfPG, 14] <- gTime + round(rpert(length(NewInfPG), 1, 5, 9))  # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
              PopMatWB[NewInfPG, 15] <- PopMatWB[NewInfPG, 14] + round(rpert(length(NewInfPG), 1, 5, 7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
            }
          }
      
          ## B - Infection of neighboring susceptible groups. 
          ## Here we model the susceptible groups that are around infectious herd i, but we also take into
          ## account that the susceptible herd could have risk from other infectious herds 
          ## Or carcasses around it or a carcass in its own pixel
          
          GroupsAtRisk <- unique(PopMatWB[(PopMatWB[ ,8] %in% tmpConEdgs & !PopMatWB[ ,2] %in% infGroups), 2])
          if(length(GroupsAtRisk) > 0){
            GroupsAtRisk <- GroupsAtRisk[GroupsAtRisk > 0]
            PixelsAtRisk <- unique(PopMatWB[PopMatWB[ ,2] %in% GroupsAtRisk, 8])
            NumInfAtRisk <- sapply(PixelsAtRisk, function(x) {tmp1 <- c(unlist(WBMat[x, 6:13]), x) # include the possibility of an infected carcass in the pixel from another group
            tmp2 <- PopMatWB[ ,8] %in% tmp1
            sum(tmp2 & PopMatWB[ ,9] %in% 2:3)}) 
            if(length(NumInfAtRisk) > 0){
              ProbInfGAR   <- 1 - (1 - CarcProbInf)^NumInfAtRisk
              StatusGAR    <- GroupsAtRisk[rbinom(length(GroupsAtRisk), 1, ProbInfGAR) == 1]
              if(length(StatusGAR) > 0){
                NumAnimPGAR <- sapply(StatusGAR, function(x) sum(PopMatWB[ ,2] == x))
                NumInfAnim  <- round(runif(length(StatusGAR), 1, 3))
                NumInfAnim[NumInfAnim > NumAnimPGAR] <- 1
                for(ss in 1:length(StatusGAR)){
                  IndexToSelect1 <- which(PopMatWB[ ,2] == StatusGAR[ss])
                  if(NumAnimPGAR[ss] > NumInfAnim[ss]) IndexToSelect2 <- sample(IndexToSelect1, NumInfAnim[ss]) else IndexToSelect2 <- IndexToSelect1
                  PopMatWB[IndexToSelect2,  9] <- 1
                  PopMatWB[IndexToSelect2, 14] <- gTime + round(rpert(length(IndexToSelect2), 1, 5, 9))  # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
                  PopMatWB[IndexToSelect2, 15] <- PopMatWB[IndexToSelect2, 14] + round(rpert(length(IndexToSelect2), 1, 5, 7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
                  
                }
              }
            }
          }
        }
            
            ## C - Groups that can be infected via carcass, while moving.
            ## C.1 - First male groups can be infected while moving (if susceptible) or infect others (if infectious)
        
            GroupsSplit    <- as.numeric(names(MovedPixMale))
            if(length(GroupsSplit) > 0){
              MovedGroups    <- unname(unlist(lapply(MovedPixMale, function(x)length(x) > 1)))
              MovedGroupsTod <- GroupsSplit[MovedGroups]
              if(length(MovedGroupsTod) > 0){
                SusMovMalG     <- sapply(MovedGroupsTod,function(x) all(PopMatWB[ ,9] == 0 & PopMatWB[,2] == x & PopMatWB[ ,12] == 1))
                SusMaleGroups  <- MovedGroupsTod[SusMovMalG]
                if(length(SusMaleGroups)>0){
                  for(i in SusMaleGroups){
                    tmp <- unname(unlist(MovedPixMale[names(MovedPixMale)==i]))
                    tmp <- tmp[tmp>0]
                    ## pixels that have carcasses
                    AnyInfPopInPix <- tmp[sapply(tmp,function(x) any(PopMatWB[,8]%in%x & PopMatWB[,9]%in%3))] 
                    if(length(AnyInfPopInPix)>0){
                      ## number of carcasses per pixel
                      NumCarPP     <- sapply(AnyInfPopInPix,function(x) sum(PopMatWB[,8]%in%x & PopMatWB[,9]==3) )
                      probInfSplit <- 1-(1-CarcProbInf)^sum(NumCarPP)
                      NewInfSG     <- rbinom(1,1,probInfSplit)
                      if(NewInfSG==1){
                        SusAnimal   <- which(PopMatWB[,2]==i & PopMatWB[,12]==1)
                        NumInfAnim  <- round(runif(1,1,3))
                        NumInfAnim[NumInfAnim>length(SusAnimal)] <- SusAnimal
                        if(SusAnimal > NumInfAnim) IndexToSelect3 <- sample(SusAnimal,NumInfAnim) 
                        else IndexToSelect3 <- SusAnimal
                        PopMatWB[IndexToSelect3,9]  <- 1
                        PopMatWB[IndexToSelect3,14] <- gTime + round(rpert(length(IndexToSelect3),1,5,9))  # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
                        PopMatWB[IndexToSelect3,15] <- PopMatWB[IndexToSelect3,14] + round(rpert(length(IndexToSelect3),1,5,7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
                      }
                    }
                  }
                }
              }
            }
       
            ## C.2 - Female Groups that can be infected while splitting 
            GroupsSplit <- as.numeric(names(PixelsMoved))
            if(length(GroupsSplit) > 0){
              MovedGroups    <- unname(unlist(lapply(PixelsMoved, function(x) length(x) > 1)))
              MovedGroupsTod <- GroupsSplit[MovedGroups]
              if(length(MovedGroupsTod) > 0){
                SusMovFMalG    <- sapply(MovedGroupsTod,function(x) all(PopMatWB[ ,9] == 0 &
                                                                        PopMatWB[ ,2] == x & 
                                                                        PopMatWB[ ,10] == 1))
                SusFemaleGroups<- MovedGroupsTod[SusMovFMalG]
                if(length(SusFemaleGroups) > 0){
                  for(i in SusFemaleGroups){
                    tmp <- unname(unlist(PixelsMoved[names(PixelsMoved) == i]))
                    tmp <- tmp[tmp > 0]
                    ## pixels that have carcasses
                    AnyInfPopInPix <- tmp[sapply(tmp, function(x) any(PopMatWB[ ,8] %in% x & PopMatWB[ ,9] %in% 3))] 
                    if(length(AnyInfPopInPix) > 0){
                      ## number of carcasses per pixel
                      NumCarPP     <- sapply(AnyInfPopInPix, function(x) sum(PopMatWB[ ,8] %in% x & PopMatWB[,9] == 3) )
                      probInfSplit <- 1 - (1 - CarcProbInf)^sum(NumCarPP)
                      NewInfSG     <- rbinom(1, 1, probInfSplit)
                      if(NewInfSG == 1){
                        SusAnimal   <- which(PopMatWB[ ,2] == i & PopMatWB[ ,10] == 1)
                        NumInfAnim  <- round(runif(1, 1, 3))
                        NumInfAnim[NumInfAnim > length(SusAnimal)] <- SusAnimal
                        if(SusAnimal > NumInfAnim) IndexToSelect3 <- sample(SusAnimal, NumInfAnim) 
                        else  IndexToSelect3 <- SusAnimal
                        PopMatWB[IndexToSelect3, 9]  <- 1
                        PopMatWB[IndexToSelect3, 14] <- gTime + round(rpert(length(IndexToSelect3), 1, 5, 9))  # the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5.3,1.3)
                        PopMatWB[IndexToSelect3, 15] <- PopMatWB[IndexToSelect3, 14] + round(rpert(length(IndexToSelect3), 1, 5, 7))# the rpert is estimated from Olsen et al. (2017) using rnorm(1000,5,0.3)
                      }
                    }
                  } 
                }
              }
            }
            
            ## D - Groups that moved due to hunting that can be infected
            
            GroupsMoved <- as.numeric(names(PixelsMovedHS))
            if(length(GroupsMoved) > 0) {
              MovedGroups    <- unname(unlist(lapply(PixelsMovedHS, function(x) length(x) > 1)))
              MovedGroupsTod <- GroupsSplit[MovedGroups]
              if(length(MovedGroupsTod) > 0){
                SusMovG    <- sapply(MovedGroupsTod, function(x) all(PopMatWB[ ,9] == 0 &
                                                                     PopMatWB[ ,2] == x & 
                                                                     PopMatWB[ ,16] == 1))
                SusceptibleGroups <- MovedGroupsTod[SusMovG]
                if(length(SusceptibleGroups) > 0){
                  for(i in SusceptibleGroups){
                    tmp <- unname(unlist(PixelsMovedHS[names(PixelsMovedHS) == i]))
                    tmp <- tmp[tmp > 0]
                    ## pixels that have carcasses
                    AnyInfPopInPix <- tmp[sapply(tmp, function(x) any(PopMatWB[ ,8] %in% x & PopMatWB[ ,9] %in% 3))]
                    if(length(AnyInfPopInPix) > 0){
                      ## number of carcasses per pixel
                      NumCarPP     <- sapply(AnyInfPopInPix, function(x) sum(PopMatWB[ ,8] %in% x & PopMatWB[ ,9] == 3))
                      probInfSplit <- 1 - (1 - CarcProbInf)^sum(NumCarPP)
                      NewInfSG     <- rbinom(1, 1, probInfSplit)
                      if(NewInfSG == 1){
                        SusAnimal   <- which(PopMatWB[ ,2] == i & PopMatWB[ ,16] == 1)
                        NumInfAnim  <- round(runif(1, 1, 3))
                        NumInfAnim[NumInfAnim > length(SusAnimal)] <- SusAnimal
                        if(SusAnimal > NumInfAnim) IndexToSelect3 <- sample(SusAnimal, NumInfAnim) 
                        else  IndexToSelect3 <- SusAnimal
                        PopMatWB[IndexToSelect3, 9]  <- 1
                        PopMatWB[IndexToSelect3, 14] <- gTime + round(rpert(length(IndexToSelect3), 1, 5, 9))  
                        PopMatWB[IndexToSelect3, 15] <- PopMatWB[IndexToSelect3, 14] + round(rpert(length(IndexToSelect3), 1, 5, 7))
                      }
                    }
                  } 
                }
              }
            }
            
            ### UPDATE
            
            
            # Update animals to Infectious status
            IndexToInf <- PopMatWB[ ,9] == 1 & PopMatWB[ ,14] == gTime
            if(sum(IndexToInf) > 0) PopMatWB[IndexToInf, 9] <- 2
            
            # Update animals to Dead/Immune state
            IndexToDOrI <- which(PopMatWB[ ,9] == 2 & PopMatWB[ ,15] == gTime)
            if(length(IndexToDOrI) > 0){
              tmpToDie <- rbinom(length(IndexToDOrI), 1, ProbDeathInf)
              if(any(tmpToDie == 1)) {
                PopMatWB[IndexToDOrI[tmpToDie == 1], 9]  <- 3
                ## Dead animals may die within their own home range 20% or any of the connecting edges 80% (splitting male groups are excluded).
                tmpDieOSG   <- IndexToDOrI[tmpToDie == 1 & (!PopMatWB[IndexToDOrI, 2] %in% GroupsToSplit[ ,1])]
                DieOutGroup <- tmpDieOSG[rbinom(length(tmpDieOSG), 1, DieInConEdge)]
                if(length(DieOutGroup) > 0){
                  EdgesDie <- sapply(DieOutGroup, function(x) {
                    tmp1 <- unlist(WBMat[PopMatWB[x, 8], 6:13])
                    tmp2 <- WBMat[tmp1, 5]
                    tmp3 <- tmp1[tmp2 %in% 1:2]
                    tmp3 <- tmp3[tmp3 > 0]
                    if(length(tmp3) > 1) sample(tmp3, 1) else(PopMatWB[x, 8])
                  })
                  ### If an animal will die in a connecting edge, its pixel home and current will change but not the group number as the probability
                  ### of infection from Carcasses within a home-range and the connecting edges is the same according to Lange et al. (2015,EFSA) .
                  PopMatWB[DieOutGroup, 7] <- EdgesDie
                  PopMatWB[DieOutGroup, 8] <- EdgesDie
                }
              }
              if(any(tmpToDie == 0)) PopMatWB[IndexToDOrI[tmpToDie == 0], 9]  <- 4
            }
            
            ## Remove dead carcasses over 28 days from death from the matrix
            indexDeadOld <- PopMatWB[ ,9] == 3 & ((gTime - PopMatWB[ ,15]) > CarcDaySur)
            
            if(sum(indexDeadOld) > 0){
              PopMatWB[indexDeadOld, ] <- 0
              PopMatWB        <- PopMatWB[order(PopMatWB[ ,2], PopMatWB[ ,3], PopMatWB[ ,5], decreasing = T), ]
            } 
            
            ## Dead animals must be removed from the splitting mechanism
            # Identify group number of animals that died today
            IndexDiedTod  <- sort(unique(PopMatWB[PopMatWB[ ,9] == 3 & PopMatWB[ ,15] == gTime, 2]))
            MIndexDiedTod <- IndexDiedTod %in% GroupsToSplit[ ,1]
            if(sum(MIndexDiedTod) > 0){
              IndexToRemSM <- IndexDiedTod[MIndexDiedTod]
              NumRemTodFSM <- tapply((PopMatWB[PopMatWB[ ,2] %in% IndexToRemSM, 9] == 3 & 
                                      PopMatWB[PopMatWB[ ,2] %in% IndexToRemSM, 15] == gTime),
                                      PopMatWB[PopMatWB[ ,2] %in% IndexToRemSM, 2], sum)
              GroupsToSplit[match(IndexToRemSM, GroupsToSplit[ ,1]), 2] <- GroupsToSplit[match(IndexToRemSM, GroupsToSplit[ ,1]), 2] - NumRemTodFSM
              if(any(GroupsToSplit[ ,2] < 1)){
                IndexDelFSM   <- which(GroupsToSplit[ ,2] < 1)
                GroupsToSplit <- GroupsToSplit[-IndexDelFSM, , drop = FALSE]
              }
            }
      } #closes ASF Dynamics
      
      
      ### REINITIATE & SUMMARISE
      
      
      ## Here re-initiate the daily matrix for males and females movements to make sure that there is no bleed from previous days
      MovedPixMale  <- as.list(matrix(0, ncol = 1))
      PixelsMoved   <- as.list(matrix(0, ncol = 1))
      PixelsMovedHS <- as.list(matrix(0, ncol = 1))
      
      ## Make the daily summary
      DayOutInfMat[gTime, Iter] <- length(unique(PopMatWB[PopMatWB[ ,9] %in% c(2:4), 2]))
      DayOutPopMat[gTime, Iter] <- length(unique(PopMatWB[ ,2]))
      DayOutAniMat[gTime, Iter] <- sum(PopMatWB[ ,2] > 0)
      InfGroupsD         <- unique(PopMatWB[PopMatWB[ ,9] %in% 2:4, 2])
      InfGroupsD         <- InfGroupsD[!InfGroupsD %in% NewInfGroups[NewInfGroups[ ,1] == Iter, 3]]
      InfCells           <- unique(PopMatWB[InfGroupsD, 8])
      InfAnimals         <- which(PopMatWB[  ,9] == 2)
      InfAnimals         <- InfAnimals[!InfAnimals %in% NewInfAnimals[NewInfAnimals[ ,1] == Iter, 3]]
      InfCarcass         <- which(PopMatWB[ ,9] == 3)
      InfCarcass         <- InfCarcass[!InfCarcass %in% NewInfCarcass[NewInfCarcass[ ,1] == Iter, 3]]
      
      cumDeathInf <- cumDeathInf + length(InfCarcass)
      
      
      if(length(InfGroupsD) > 0){
        NewInfGroups    <- rbind(NewInfGroups, cbind(Iter, gTime, InfGroupsD, InfCells))
      }
      if(length(InfAnimals) > 0){
        NewInfAnimals   <- rbind(NewInfAnimals, cbind(Iter, gTime, InfAnimals))
      }
      if(length(InfCarcass) > 0){
        NewInfCarcass <- rbind(NewInfCarcass, cbind(Iter, gTime, InfCarcass))
      }
      
      ## A variable to count whether the disease faded out this iteration and started up again. 
      ## This is counted once/iteration despite in while loop
      
      TMPOutToday <- TMPOutYesterday & sum(PopMatWB[ ,9] %in% 1:2) > 0
      if(TMPOutToday & OnlyOnce){
        FreqRelapse[Iter] <- 1
        OnlyOnce <- FALSE
      }
      
      TMPOutYesterday <- gTime > TimeSeedInf & sum(PopMatWB[ ,9] %in% 1:2) == 0 & sum(PopMatWB[ ,9] == 3) > 0 
     
     print(c(runID, Iter, gTime))
      
     # } #While (Criteria & (gTime ...
     
     
     if(Iter == MaxIterations){  
       
       PopMatWB <- as.data.frame(PopMatWB)
       
       # #Daily Population Count
       #Popcount <- table(PopMatWB[,"Current_pixel"])
       Popcount <- PopMatWB %>% group_by(Current_pixel) %>%
         dplyr::summarise(n())
       Popcount <- t(Popcount)
       #Final population count (number of animals) each day (1825 rows) and in each cell (1 383 columns)
       #Current pixel that matches Popcount in Fpopcount (column) will have the sum of individuals 
       Fpopcount[gTime, Popcount[1, -1]] = Popcount[2, -1]
       #Fpopcount[gTime, as.numeric(names(Popcount))[-1]] = Popcount[-1]
       
       # if(Criteria == F | gTime==MaxDays){
       #   NAMEPOPC <- paste(runID, 'PopCount.txt', sep ='-')
       #   write.table(Fpopcount, NAMEPOPC, sep=' ', col.names=T, row.names = T)
       # }
     
       
       if(TimeSeedInf > 0){
         
         #Daily Infected Count
         Infcount <- PopMatWB %>%
           group_by(Current_pixel, Infect_status) %>%
           filter(Infect_status == 2) %>%
           dplyr::summarise(n())
         Infcount <- t(Infcount)
         # if(gTime==TimeSeedInf){
         # print(which(Infcount>0))
         
         #Final Infected count (number of infected animals) each day in each cell
         Finfcount[gTime, Infcount[1, ]] = Infcount[3, ]
         
         #Carcass Daily Count
         CarcassCount <- PopMatWB %>%
           group_by(Current_pixel, Infect_status) %>%
           filter(Infect_status == 3) %>%
           dplyr::summarise(n())
         CarcassCount <- t(CarcassCount)
         Fcarcount[gTime, CarcassCount[1, ]] = CarcassCount[3 ,]
         
         #Immune Count
         Immunecount <- PopMatWB %>%
           group_by(Current_pixel, Infect_status) %>%
           filter(Infect_status == 4) %>%
           dplyr::summarise(n())
         Immunecount <- t(Immunecount)
         Fimmunecount[gTime, Immunecount[1, ]] = Immunecount[3, ]
         
         # if(Criteria == F | gTime==MaxDays){
         #   
         #   NAMEINFC <- paste(runID, 'InfectedCount.txt', sep='-')
         #   NAMECARC <- paste(runID, 'CarcassCount.txt', sep='-')
         #   NAMEIMC <- paste(runID, 'ImmuneCount.txt', sep='-')
         #   
         #   write.table(Finfcount, NAMEINFC, sep = ' ', col.names = T, row.names = T)
         #   write.table(Fcarcount, NAMECARC, sep = ' ', col.names = T, row.names = T)
         #   write.table(Fimmunecount, NAMEIMC, sep = ' ', col.names = T, row.names = T)
         #   
         # }
         
       }
       
     }
     
  } #While(criteria...)
    
    
    #Save Popcount file at the end of the big loop. If there is infection, Maxdays could be less than 5 years...
      NAMEPOPC <- paste(runID, 'PopCount.txt', sep = '-')
      write.table(Fpopcount, NAMEPOPC, sep=' ', col.names = T, row.names = T)
    
    
    #Save Infected count, carcass and immune only if there is infection
    if(TimeSeedInf>0){
      
      NAMEINFC <- paste(runID, 'InfectedCount.txt', sep = '-')
      NAMECARC <- paste(runID, 'CarcassCount.txt', sep = '-')
      NAMEIMC  <- paste(runID, 'ImmuneCount.txt', sep = '-')
      
      write.table(Finfcount, NAMEINFC, sep = ' ', col.names = T, row.names = T)
      write.table(Fcarcount, NAMECARC, sep = ' ', col.names = T, row.names = T)
      write.table(Fimmunecount, NAMEIMC, sep = ' ', col.names = T, row.names = T)
      
    }
    
    
    ## Make the summaries per iteration
    EpDuration[Iter]            <- gTime - TimeSeedInf
    cumDeath[Iter]              <- cumDeathPar + sum(PopMatWB[ ,2] > 0)
    InfCarcass[Iter]            <- cumDeathInf 
    
    #print(Iter)
    
  } #Close for(Iter in...)
  
  NAMENA    <- paste(runID, "FNewInfAnimals.txt", sep = "-")
  NAMECA    <- paste(runID, "FNewInfCarcass.txt", sep = "-")
  NAMENI    <- paste(runID, "FNewInfGroups.txt" , sep = "-")
  NAMEIG    <- paste(runID, "FDayOutInfMat.txt" , sep = "-")
  NAMETG    <- paste(runID, "FDayOutPopMat.txt" , sep = "-")
  NAMETA    <- paste(runID, "FDayOutAniMat.txt" , sep = "-")
  
  #(NewInfAnimals, 'NAMENA', append = T, sep =" ", col.names=F, row.names=F)
  write.table(NewInfAnimals, NAMENA, sep = " ", col.names = T, row.names = T)
  write.table(NewInfCarcass, NAMECA, sep = " ", col.names = T, row.names = T)
  write.table(NewInfGroups,  NAMENI, sep = " ", col.names = T, row.names = T)
  write.table(DayOutInfMat,  NAMEIG, sep = " ", col.names = T, row.names = T)
  write.table(DayOutPopMat,  NAMETG, sep = " ", col.names = T, row.names = T)
  write.table(DayOutAniMat,  NAMETA, sep = " ", col.names = T, row.names = T)
  
  OutPutList <- cbind(EpDuration, FreqRelapse, cumDeath, InfCarcass)
  NAMEOL     <- paste(runID, "FOutPutList.txt", sep = "-")
  write.table(OutPutList, NAMEOL, sep = " ", col.names = T, row.names = T)
  
  #Save data frame with the main characteristics of the model, for statistic analysis/comparison between models
  #List with the Output list result as first element and a data frame with the model main characteristics as second element
    
    OutPutParamsList <- list(OutPutList,
                             NewInfCarcass,
                             Params = data.frame(runID = runID,
                                                 Hunt = paste(ProbHarvest*100,'% +',ProbHarvestAM*100, '%'),
                                                 MovHunt = paste(ProbMovHunt*100, "%"),
                                                 Cells = "9 km2",
                                                 InitGroup = InitialOccCells,
                                                 HabitatProb = paste0('HC',length(HabitatProb)),
                                                 HuntPiglet = "No",
                                                 SurvivalP = SurvivalProbAdM,
                                                 MaleReset = 'No',
                                                 FemaleProp = '3-6',
                                                 GroupMoveHunt = 'Yes',
                                                 SeasonalCarcass = 'Yes'
                                                 ))
  save(OutPutParamsList, file = paste0(runID, "-FOutPutList.Rda"))
  
}



# 
# library(doParallel)
# 
# cl <- parallel::makeCluster(length(Scenarios$Sc))
# doParallel::registerDoParallel(cl)

#For Scenario 1 - 4 km

HC3 = c(0.75, 0.20, 0.05)
HC2 = c(0.75,0.25)

Scenarios <- data.frame(Sc = c(1:42),
                        ProbHarvest   = c(rep(0, 2), rep(0.40, 2), 0, rep(0.40, 20), 0, 0.70, 0.90, rep(0.40, 13), 0),
                        ProbHarvestAM = c(rep(0, 2), 0.40, 0.70, 0, rep(0.70, 16), rep(0.40, 4), 0, 0.70, 0.90, rep(0.70, 6), 0.60, 0.50, rep(0.60, 5), 0),
                        ProbMovHunt   = c(rep(0, 5), 0.30, 0.70, rep(0.30, 5), rep(0.70, 5), rep(0.30, 2), rep(0.70, 4), rep(0.30, 2), 0, rep(0.70, 12), 0, rep(0.30, 2), 0),
                        Cells         = c(4, rep(9, 27), 4, rep(9, 13)),
                        InitGroup     = c(rep(50, 7), 100, 200, 300, 400, 500, 100, 200, 300, 400, 500, 400, 500, 400, 500, 400, 500, 400, rep(500, 18)),
                        HabitatProb   = c(rep("HC2", 4), rep("HC3", 13), rep("HC2", 8), rep('HC3', 17)),
                        HuntPigglet   = c(rep('NA', 2), rep('Yes', 2), 'NA', rep('Yes', 24), rep('No', 12), 'NA'),
                        SurvivalP     = c(rep(0.60, 30), 0.70, 0.80, rep(0.90, 10)),
                        MaleReset     = c(rep('Yes', 33), rep('No', 9)),
                        FemaleProp    = c(rep("3-4", 36), "3-5", rep("3-6", 5)),
                        GroupMoveHunt = c(rep("No", 38), rep("Yes", 3), "No"),
                        CarcSeason    = c(rep("No", 39), rep("Yes",3))
                        )




# foreach (i = c(26:28), .export = c('HC2','HC3'), .packages = "dplyr")  %dopar% {
# 
# HabitatProb = get(as.character(Scenarios[Scenarios$Sc == i,"HabitatProb"]))
# 
#   WBModel(ProbHarvest = Scenarios[Scenarios$Sc == i,"ProbHarvest"],
#           ProbHarvestAM= Scenarios[Scenarios$Sc == i,"ProbHarvestAM"],
#           ProbMovHunt = Scenarios[Scenarios$Sc == i,"ProbMovHunt"],
#           InitialOccCells = Scenarios[Scenarios$Sc == i,"InitGroup"],
#           HabitatProb = HabitatProb,
#           runID = paste0("WB_MovH_ASF_Scenario_", i)
#          )
# }
# 
#  parallel::stopCluster(cl)

 # Scenarios$Sc
 
for (i in 41:42) {
  
  HabitatProb = get(as.character(Scenarios[Scenarios$Sc == i,"HabitatProb"]))

  WBModel(ProbHarvest      = Scenarios[Scenarios$Sc == i, "ProbHarvest"],
          ProbHarvestAM    = Scenarios[Scenarios$Sc == i, "ProbHarvestAM"],
          ProbMovHunt      = Scenarios[Scenarios$Sc == i, "ProbMovHunt"],
          InitialOccCells  = Scenarios[Scenarios$Sc == i, "InitGroup"],
          HabitatProb      = HabitatProb,
          SurvivalProbAdF  = Scenarios[Scenarios$Sc == i, "SurvivalP"],
          SurvivalProbAdM  = Scenarios[Scenarios$Sc == i, "SurvivalP"],
          SurvivalProbSAdF = Scenarios[Scenarios$Sc == i, "SurvivalP"],
          SurvivalProbSAdM = Scenarios[Scenarios$Sc == i, "SurvivalP"],
          runID            = paste0("WB_Model_Dep64_Scenario_", i))
}


# for (InitialOccCells in seq(100, 500, 100)) {
#   for(ProbMovHunt in c(0.3, 0.7)){
#     WBModel(ProbMovHunt = ProbMovHunt, InitialOccCells = InitialOccCells)
#   }
# }
 

