library(reshape2)
library(gtools) 
library(lme4)
library(plyr)
library(dplyr)
library(lattice)
library(car)
library(HH)
library(agricolae) 
library(multcomp) 
library(MASS)
library(heplots) 
library (coin) 
library(ggplot2)
library(Hmisc)
library(lmerTest)
library(LMERConvenienceFunctions)
library(languageR)
library(doBy)
library(RColorBrewer)
library("eyetrackingR")
library(pbapply)
install.packages('doParallel')
library(doParallel)
#amusing sounds you can use to let you know when your code is finished, if you don't want to use just get rid of anything with BRRR
if(!require(devtools)) {install.packages("devtools")}
devtools::install_github("brooke-watson/BRRR")

library(doMC)
library(doSNOW)

registerDoSNOW(makeCluster(4, type = "SOCK"))

doMC::registerDoMC(cores=4) # or however many cores you have access to
cores <- detectCores()
cores
getwd()

########################################
###
###   Load E-Prime Data for Trials   ###
###
########################################
#
#   Read in Samples Master File
#   Made by MVOR_MakeMaster_Samples.r
#   S. Olsen
##  Modified PG - June 29, 2018
##  Modified K. Stockwell June 2019
#   Meshes SMI with EMerge files
#   Builds AOI data for EyeTrackingR
##
########################################

#This prevents Scientific Notation from being invoked
#  Should not be an issue - SMI merging script recalculates times
options(scipen=999)

# STOP!! Open and remove random "x" column here! (don't have to, x column is just annoying)
#  Delete "UseSample column as well unless done in prev. script (WHY? Kayden)
Samples <- read.csv(file = "ItemTEST_MV2_Full_Samples_Kids.csv", header = TRUE)
#Remove unneeded row with Mask.png message
Samples <- subset(Samples, L.Raw.X..px. !="# Message: Mask.png")

#unique(Samples$ExperimentName)
unique(Samples$Trial)
unique(Samples$Type)
unique(Samples$Subject)
nrow(Samples)

### Eprime File ###
### Read in the csv file 
EPrime <- read.csv(file = "../MV2_Kids_EMerge.csv", header=TRUE)
nrow(EPrime)
colnames(EPrime)

## For KIDS ONLY- AGE RANGE ANALYSIS Load age qualifiers and merge into SMI file
AgeList <- read.csv(file = "../AgeList.csv", header = TRUE)
unique(AgeList$Subject)
Samples <- merge(AgeList, Samples, by=c("Subject"), all=T)
colnames(Samples)

###Reduce the number of variables
EPrime <- EPrime[ , c("ExperimentName", "Subject", "Trial", "ArrayFile", "ComplementFile", "ComplementLocation",
                      "Condition", "DiffExemplarFile", "DiffExemplarLocation", "DistractorFile", "DistractorLocation", "FvPMatch",
                      "ItemNumber", "ListNumber", "Q1File", "Q1Item", "Q2File", "Q2Item", "Q3File", "Q3Item", "Q4File", "Q4Item", 
                      "TargetFile", "TargetLocation", "TrialType")]


#Remove practice trials
EPrime <- EPrime[EPrime$Trial > 4, ]

unique(EPrime$ExperimentName)
unique(EPrime$Subject)
unique(Samples$Subject)
unique(EPrime$Trial)
nrow(EPrime)

###Merge the eye tracking and EPrime files###
MVOR2_MERG <- merge(EPrime, Samples, by=c("Subject", "Trial"), all=T)
nrow(MVOR2_MERG)
unique(MVOR2_MERG$Subject)
unique(MVOR2_MERG$ExperimentName)
nrow(MVOR2_MERG)

# Sort MVOR2 (using 'arrange' from ddply)
newData <- arrange(MVOR2_MERG, MVOR2_MERG$Subject, MVOR2_MERG$Trial, MVOR2_MERG$Time)
MVOR2_MERG <- newData
write.csv(MVOR2_MERG, file = "MV2_Full_EPRIME-SMI.csv")

###########################################

##########Experimental Trials##############
########AOIS and Timestamp Columns#########

###########################################

##Read in the file if needed
#MVOR2_MERG <- read.csv(file = "MV2_Full_EPRIME-SMI.csv", header=TRUE)
#colnames(MVOR2_MERG)
#head(MVOR2_MERG, n = 5)

unique(MVOR2_MERG$Subject)
unique(MVOR2_MERG$ExperimentName)
unique(MVOR2_MERG$TrialType)


#################################
### STOP AND SORT!! #############
#################################
# 'order' is the quick-sort function
#MV2_Prime <- MV2_Prime[
#  order( MV2_Prime[,"Subject"], MV2_Prime[,"Trial"], MV2_Prime[, "PrimeOnset"] ),
#  ]

ldf <- MVOR2_MERG
ldf$PrimeOrArray <- 0
#Label rows as either Prime or Array
cat ("Number of Files Opened = ", nrow(ldf))
for (r in 1:nrow(ldf)) {
  #  Add flag for "prime" and in-stimulus-array" marker
    k = r - 1
    if ((grepl('# Message: p', ldf [r, ("L.Raw.X..px.")]))
        | (grepl('# Message: f', ldf [r, ("L.Raw.X..px.")]))) {
        ldf [r, ("PrimeOrArray")] <- "Prime"
    } else if ((ldf [r, ("L.Raw.X..px.")] == '# Message: ArrayOnset')) {
              ldf [r, ("PrimeOrArray")] <- "Array"
    } else {
      ldf [r, ("PrimeOrArray")] <- ldf [k, ("PrimeOrArray")]
    }
    if ((r %% 5000) == 0) { 
      cat ("Iteration = ", r)
    }
  }


write.csv(ldf, file = "MV2_Full_PrimeorArray_Sorted.csv")
#Read in if needed
#MV2_FULL <- read.csv(file = "MV2_Full_PrimeorArray_Sorted.csv", header = TRUE)

unique (MV2_FULL$ExperimentName)
unique (MV2_FULL$Subject)
unique (MV2_FULL$TrialType)

###Seperate into prime and an array data frames
MV2_Prime <- subset(MV2_FULL, PrimeOrArray != "Array")
MV2_Array <- subset(MV2_FULL, PrimeOrArray != "Prime")


###remove messages from data set
MV2_Prime_SMP <- MV2_Prime[MV2_Prime$Type != 'MSG', ]

##################################

##  ADD Prime AOI               ##
##  Might be able to            ##
##  Use Eyetracking R function  ##

##################################
#
#  FIRST - Create a 'SampleLost' column for EyeTrackingR
# Creating a column that specifies whether the eye tracker lost the eye for a given sample
#
MV2_Prime_SMP$SampleLost <- 999

for (i in 1:nrow(MV2_Prime_SMP)) {
  #Target
  if ((MV2_Prime_SMP[i, c("L.POR.X..px.")] == 0) & (MV2_Prime_SMP[i, c("L.POR.Y..px.")] == 0)) {
    MV2_Prime_SMP[i, c("SampleLost")] <- 1
  } else {
    MV2_Prime_SMP[i, c("SampleLost")] <- 0
  }
  if ((i %% 1000) == 0) {
    print(i)
  }
}

#write.csv(MV2_Prime_SMP, file = "MV2_Prime_SMP_Only.csv")
#MV2_Prime_SMP <- read.csv(file = "MV2_Prime_SMP_Only.csv", header = TRUE)
MV2_Prime_SMP_EXP <- MV2_Prime_SMP[MV2_Prime_SMP$TrialType == "Experimental", ]

#######For Loop marking 1 if the sample is in the Prime AOI, and 0 if it is not. Must do on Prime data, not Array!######
#Creates column that will denote if sample is in Prime AOI
MV2_Prime_SMP_EXP$PrimeAOI <- 0

for (i in 1:nrow(MV2_Prime_SMP_EXP)) {
  x = MV2_Prime_SMP_EXP[i, ("L.POR.X..px.")]
  y = MV2_Prime_SMP_EXP[i, ("L.POR.Y..px.")]
  prime = 0
  #coordinates from Illustrator. Confirm this is correct 6/18/19 KAYDEN
  if (((x >= 686) & (x <= 1234)) & ((y >= 310) & (y <= 770))) { prime = 1 }
  # Writes value back into dataframe
  MV2_Prime_SMP_EXP[i,"PrimeAOI"] <- prime
  
  if ((i %% 1000) == 0) {print(i)}
}

#######################################
#####Make time stamp number column#####
#######################################

###make column with time of first sample
MV2_Prime_SMP_EXP_T <- ddply(MV2_Prime_SMP_EXP, c("Subject", "Trial"), transform, TimeOfFirstSample=min(Time))

###make column with sample time starting from the time of the first sample recorded after appearance of the array image
MV2_Prime_SMP_EXP_T$SampleTime <- (MV2_Prime_SMP_EXP_T$Time - MV2_Prime_SMP_EXP_T$TimeOfFirstSample) 

###Create a column (TimeStamp) which labels each incrementing sample time as a time stamp number starting with one, and increasing sequentially by one 
MV2_Prime_SMP_EXP_T <- ddply(MV2_Prime_SMP_EXP_T, c("Subject", "Trial"), transform, TimeStamp=seq(1, (length(SampleTime)), by=1))

###Write file with timestamp number###
##DO NOT REWRITE FILE!!!##
write.csv(MV2_Prime_SMP_EXP_T, file = "MV2_SMP_Prime_EXP_ONLY_TSTAMP_Kids.csv")
#Read in if needed
#MV2_Prime_SMP_EXP_T <- read.csv(file = "MV2_SMP_Prime_EXP_ONLY_TSTAMP_Kids.csv", header = TRUE)

#######################################
#########Prime AOI Examination#########
#######################################

#Creates dataframe showing percent of time slices subject is looking in the Prime AOI for each trial
PercentInAOIBySubj <- ddply(MV2_Prime_SMP_EXP_T, c("Subject", "Trial"), summarise, PctInAOI=sum(PrimeAOI)/length(TimeStamp))

#Creates dataframe showing percent of time slices subject is looking in the Prime AOI for each Prime stimulus
PercentInAOIByImage <- ddply(MV2_Prime_SMP_EXP_T, c("TargetFile"), summarise, PctInAOI=sum(PrimeAOI)/length(TimeStamp))

#What percent of time slices subject is looking in the Prime AOI for trial and for Prime stimulus excluding SampleLost samples
MV2_Prime_SMP_EXP_T_NoSL <- subset(MV2_Prime_SMP_EXP_T, SampleLost != 1)
PercentInAOIBySubj_NoSL <- ddply(MV2_Prime_SMP_EXP_T_NoSL, c("Subject", "Trial"), summarise,
                                         PctInAOI=sum(PrimeAOI)/length(TimeStamp))
PercentInAOIByImage_NoSL <- ddply(MV2_Prime_SMP_EXP_T_NoSL, c("TargetFile"), summarise, PctInAOI=sum(PrimeAOI)/length(TimeStamp))

#What percentage of out of Prime AOI looks are off screen for by subject each trial (SampleLost)
PercentOutAOIBySubj_SL <- ddply(MV2_Prime_SMP_EXP_T, c("Subject", "Trial"), summarise,
                                  PctOutAOI_SL=sum(PrimeAOI == 0 & SampleLost == 1)/sum(PrimeAOI == 0))
#What percentage of out of Prime AOI looks are off screen for by Prime stimulus (SampleLost)
#Not sure I believe this output... 6/19/19 KAYDEN
PercentOutAOIByImage_SL <- ddply(MV2_Prime_SMP_EXP_T, c("TargetFile"), summarise,
                                  PctOutAOI_SL=sum(PrimeAOI == 0 & SampleLost == 1)/sum(PrimeAOI == 0))

#Remove trials for each subject that are less than 50% of samples in Prime AOI
###Merge files
MV2_Prime_SMP_EXP_T_Merged <- merge(MV2_Prime_SMP_EXP_T, PercentInAOIBySubj, by=c("Subject", "Trial"), all=T)
MV2_Prime_SMP_EXP_T_Cleaned <- subset(MV2_Prime_SMP_EXP_T_Cleaned, PctInAOI > 0.5)

# 'order' is the quick-sort function
MV2_Prime_SMP_EXP_T_Cleaned <- MV2_Prime_SMP_EXP_T_Cleaned[
  order(MV2_Prime_SMP_EXP_T_Cleaned[ , "Subject"], MV2_Prime_SMP_EXP_T_Cleaned[ , "Trial"]), ]

#Find out what trials remain for each subject (neater way to do this would be nice)
Retained_SubjTrials <- MV2_Prime_SMP_EXP_T_Cleaned %>% distinct(Subject, Trial)

###Write file with timestamp number###
write.csv(MV2_Prime_SMP_EXP_T_Cleaned, file = "MV2_SMP_Prime_EXP_Cleaned_TSTAMP_Kids.csv")
#Read in if needed
#MV2_Prime_SMP_EXP_T_Cleaned <- read.csv(file = "MV2_SMP_Prime_EXP_Cleaned_TSTAMP_Kids.csv", header = TRUE)

#######################################
##########Quad AOI Assignment##########
#######################################

#######For Loop marking 1 if the sample is in one of the 4 array AOIs, and 0 if it is not. Must do on Array data, not Prime!######

###remove messages from Array data set
MV2_Array_SMP <- MV2_Array[MV2_Array$Type != 'MSG', ]

#  FIRST - Create a 'SampleLost' column for EyeTrackingR
# Creating a column that specifies whether the eye tracker lost the eye for a given sample
#
MV2_Array_SMP$SampleLost <- 999

for (i in 1:nrow(MV2_Array_SMP)) {
  #Target
  if ((MV2_Array_SMP[i, c("L.POR.X..px.")] == 0) & (MV2_Array_SMP[i, c("L.POR.Y..px.")] == 0)) {
    MV2_Array_SMP[i, c("SampleLost")] <- 1
  } else {
    MV2_Array_SMP[i, c("SampleLost")] <- 0
  }
  if ((i %% 1000) == 0) {
    print(i)
  }
  #Code makes sound when done, last line is hardcoded for now
  if (i == 445337) {
    BRRR::skrrrahh(2)
  }
}

#Create file with ONLY Array experimental trials
MV2_Array_SMP_EXP <- MV2_Array_SMP[MV2_Array_SMP$TrialType == "Experimental", ]

##
##  Quandrant ID = QuadID
##  Sarah used binary flags - use 1-4 instead?  Faster...
##  Quad indicates location of the look
##  Set to 999 first to indicate NO AOI
##
MV2_Array_SMP_EXP$Target <- 0   #  These first 4 are the AOI columns - set to 1 if present, 0 otherwise
MV2_Array_SMP_EXP$Complement <- 0
MV2_Array_SMP_EXP$DiffExemplar <- 0
MV2_Array_SMP_EXP$Distractor <- 0
MV2_Array_SMP_EXP$InQuad <- 0

#  When doing a data check:
#  DO NOT Check with all rows - change for a test to about 1000!  (Takes a long time)
#  Future - convert to lapply
for (i in 1:nrow(MV2_Array_SMP_EXP)) { 
  x = MV2_Array_SMP_EXP[i, c("L.POR.X..px.")]
  y = MV2_Array_SMP_EXP[i, c("L.POR.Y..px.")]
  quad = 999
  
  if ( ((x >= 1103) & (x <= 1528)) & ((y >= 0) & (y <= 425)) ) { quad = 1 }
  if ( ((x >= 421) & (x <= 846)) & ((y >= 0) & (y <= 425)) ) { quad = 2 }
  if ( ((x >= 421) & (x <= 846)) & ((y >= 655) & (y <= 1080)) ) { quad = 3 }
  if ( ((x >= 1103) & (x <= 1528)) & ((y >= 655) & (y <= 1080)) ) { quad = 4 }
  
  MV2_Array_SMP_EXP[i,"InQuad"] <- quad  # Writes value back into dataframe!
  
  if (quad == 1) {
    if (MV2_Array_SMP_EXP[i,"TargetLocation"] == 'Q1') { MV2_Array_SMP_EXP[i,"Target"] = 1}
    if (MV2_Array_SMP_EXP[i,"ComplementLocation"] == 'Q1') { MV2_Array_SMP_EXP[i,"Complement"] = 1}
    if (MV2_Array_SMP_EXP[i,"DiffExemplarLocation"] == 'Q1') { MV2_Array_SMP_EXP[i,"DiffExemplar"] = 1}
    if (MV2_Array_SMP_EXP[i,"DistractorLocation"] == 'Q1') { MV2_Array_SMP_EXP[i,"Distractor"] = 1}
  }
  
  if (quad == 2) {
    if (MV2_Array_SMP_EXP[i,"TargetLocation"] == 'Q2') { MV2_Array_SMP_EXP[i,"Target"] = 1}
    if (MV2_Array_SMP_EXP[i,"ComplementLocation"] == 'Q2') { MV2_Array_SMP_EXP[i,"Complement"] = 1}
    if (MV2_Array_SMP_EXP[i,"DiffExemplarLocation"] == 'Q2') { MV2_Array_SMP_EXP[i,"DiffExemplar"] = 1}
    if (MV2_Array_SMP_EXP[i,"DistractorLocation"] == 'Q2') { MV2_Array_SMP_EXP[i,"Distractor"] = 1}
  }
  
  if (quad == 3) {
    if (MV2_Array_SMP_EXP[i,"TargetLocation"] == 'Q3') { MV2_Array_SMP_EXP[i,"Target"] = 1}
    if (MV2_Array_SMP_EXP[i,"ComplementLocation"] == 'Q3') { MV2_Array_SMP_EXP[i,"Complement"] = 1}
    if (MV2_Array_SMP_EXP[i,"DiffExemplarLocation"] == 'Q3') { MV2_Array_SMP_EXP[i,"DiffExemplar"] = 1}
    if (MV2_Array_SMP_EXP[i,"DistractorLocation"] == 'Q3') { MV2_Array_SMP_EXP[i,"Distractor"] = 1}
  }
  
  if (quad == 4) {
    if (MV2_Array_SMP_EXP[i,"TargetLocation"] == 'Q4') { MV2_Array_SMP_EXP[i,"Target"] = 1}
    if (MV2_Array_SMP_EXP[i,"ComplementLocation"] == 'Q4') { MV2_Array_SMP_EXP[i,"Complement"] = 1}
    if (MV2_Array_SMP_EXP[i,"DiffExemplarLocation"] == 'Q4') { MV2_Array_SMP_EXP[i,"DiffExemplar"] = 1}
    if (MV2_Array_SMP_EXP[i,"DistractorLocation"] == 'Q4') { MV2_Array_SMP_EXP[i,"Distractor"] = 1}
  }
  
  if ((i %% 1000) == 0) {  
    print(i)
  }
    #Code makes sound when done
    if (i == nrow(MV2_Array_SMP_EXP)) {
      BRRR::skrrrahh(3)
  }
}

unique(MV2_Array_SMP_EXP$InQuad)
unique(MV2_Array_SMP_EXP$TrialType)

#######################################
#####Make time stamp number column#####
#######################################

###make column with time of first sample
MV2_Array_SMP_EXP <- ddply(MV2_Array_SMP_EXP, c("Subject", "Trial"), transform, TimeOfFirstSample=min(Time))

###make column with sample time starting from the time of the first sample recorded after appearance of the array image
MV2_Array_SMP_EXP$SampleTime <- (MV2_Array_SMP_EXP$Time - MV2_Array_SMP_EXP$TimeOfFirstSample) 

###Create a column (TimeStamp) which labels each incrementing sample time as a time stamp number starting with one, and increasing sequentially by one 
MV2_Array_SMP_EXP <- ddply(MV2_Array_SMP_EXP, c("Subject", "Trial"), transform, TimeStamp=seq(1, (length(SampleTime)), by=1))

###Write file with timestamp number###
##DO NOT REWRITE FILE!!!##
write.csv(MV2_Array_SMP_EXP, file = "MV2_Array_SMP_EXP_ONLY_TSTAMP_Kids.csv")

#Create file with only trials where they looked at Prime more than 50% of time
MV2_Array_SMP_EXP_Trimmed <- MV2_Array_SMP_EXP %>% semi_join(Retained_SubjTrials)
#Match to Retained_SubjTrials
MV2_Array_SMP_EXP_Trimmed %>% distinct(Subject, Trial)

###Write trimmed file with timestamp number###
##DO NOT REWRITE FILE!!!##
write.csv(MV2_Array_SMP_EXP_Trimmed, file = "MV2_Array_SMP_EXP_ONLY_TSTAMP_Kids_Trimmed.csv")

###########################################
##
####        Filler Trials           #######
####  AOIS and Timestamp Columns    #######
##
###########################################
MV2_SMP <- read.csv(file = "MV2_SMP_EXP_ONLY_TSTAMP_Kids.csv", header=TRUE)

###remove messages from data set
MV2_SMP <- MV2_FULL[MV2_FULL$Type != 'MSG', ]

#Create file with only filler trials
MV2_SMP_FILLER <- MV2_SMP[MV2_SMP$TrialType == 'Filler', ]
unique (MV2_SMP_FILLER$TrialType)

MV2_SMP_FILLER$SampleLost <- 999

for (i in 1:nrow(MV2_SMP_FILLER)) {
  #Target
  if ((MV2_SMP_FILLER[i, c("L.POR.X..px.")] == 0) & (MV2_SMP_FILLER[i, c("L.POR.Y..px.")] == 0)) {
    MV2_SMP_FILLER[i, c("SampleLost")] <- 1
  } else {
    MV2_SMP_FILLER[i, c("SampleLost")] <- 0
  }
  if ((i %% 1000) == 0) {
    print(i)
  }
}

##
##  Quandrant ID = QuadID
##  Sarah used binary flags - use 1-4 instead?  Faster...
##  Quad indicates location of the look
##  Set to 999 first to indicate NO AOI
##
MV2_SMP_FILLER$Target <- 0   #  These first 4 are the AOI columns - set to 1 if present, 0 otherwise
MV2_SMP_FILLER$Complement <- 0
MV2_SMP_FILLER$DiffExemplar <- 0
MV2_SMP_FILLER$Distractor <- 0
MV2_SMP_FILLER$InQuad <- 0

for (i in 1:nrow(MV2_SMP_FILLER)) { #  DO NOT Check with all rows - change for a test to about 1000!
  x = MV2_SMP_FILLER[i, c("L.POR.X..px.")]
  y = MV2_SMP_FILLER[i, c("L.POR.Y..px.")]
  quad = 999
  
  if ( ((x >= 1103) & (x <= 1528)) & ((y >= 0) & (y <= 425)) ) { quad = 1 }
  if ( ((x >= 421) & (x <= 846)) & ((y >= 0) & (y <= 425)) ) {  quad = 2 }
  if ( ((x >= 421) & (x <= 846)) & ((y >= 655) & (y <= 1080)) ) { quad = 3 }
  if ( ((x >= 1103) & (x <= 1528)) & ((y >= 655) & (y <= 1080)) ) { quad = 4 }
  
  MV2_SMP_FILLER[i,"InQuad"] <- quad  # Write value back into dataframe!
  
  if (quad == 1) {
    if (MV2_SMP_FILLER[i,"TargetLocation"] == 'Q1') { MV2_SMP_FILLER[i,"Target"] = 1}
    if (MV2_SMP_FILLER[i,"ComplementLocation"] == 'Q1') { MV2_SMP_FILLER[i,"Complement"] = 1}
    if (MV2_SMP_FILLER[i,"DiffExemplarLocation"] == 'Q1') { MV2_SMP_FILLER[i,"DiffExemplar"] = 1}
    if (MV2_SMP_FILLER[i,"DistractorLocation"] == 'Q1') { MV2_SMP_FILLER[i,"Distractor"] = 1}
  }
  
  if (quad == 2) {
    if (MV2_SMP_FILLER[i,"TargetLocation"] == 'Q2') { MV2_SMP_FILLER[i,"Target"] = 1}
    if (MV2_SMP_FILLER[i,"ComplementLocation"] == 'Q2') { MV2_SMP_FILLER[i,"Complement"] = 1}
    if (MV2_SMP_FILLER[i,"DiffExemplarLocation"] == 'Q2') { MV2_SMP_FILLER[i,"DiffExemplar"] = 1}
    if (MV2_SMP_FILLER[i,"DistractorLocation"] == 'Q2') { MV2_SMP_FILLER[i,"Distractor"] = 1}
  }
  
  if (quad == 3) {
    if (MV2_SMP_FILLER[i,"TargetLocation"] == 'Q3') { MV2_SMP_FILLER[i,"Target"] = 1}
    if (MV2_SMP_FILLER[i,"ComplementLocation"] == 'Q3') { MV2_SMP_FILLER[i,"Complement"] = 1}
    if (MV2_SMP_FILLER[i,"DiffExemplarLocation"] == 'Q3') { MV2_SMP_FILLER[i,"DiffExemplar"] = 1}
    if (MV2_SMP_FILLER[i,"DistractorLocation"] == 'Q3') { MV2_SMP_FILLER[i,"Distractor"] = 1}
  }
  
  if (quad == 4) {
    if (MV2_SMP_FILLER[i,"TargetLocation"] == 'Q4') { MV2_SMP_FILLER[i,"Target"] = 1}
    if (MV2_SMP_FILLER[i,"ComplementLocation"] == 'Q4') { MV2_SMP_FILLER[i,"Complement"] = 1}
    if (MV2_SMP_FILLER[i,"DiffExemplarLocation"] == 'Q4') { MV2_SMP_FILLER[i,"DiffExemplar"] = 1}
    if (MV2_SMP_FILLER[i,"DistractorLocation"] == 'Q4') { MV2_SMP_FILLER[i,"Distractor"] = 1}
  }
  
  if ((i %% 1000) == 0) {  print(i) }
}

unique(MV2_SMP_FILLER$InQuad)
unique(MV2_SMP_FILLER$TrialType)

#######################################
#####Make time stamp number column#####
#######################################

###make column with time of first sample
system.time(MV2_SMP_FILLER <- ddply(MV2_SMP_FILLER, c("Subject", "Trial"), transform, TimeOfFirstSample=min(Time)))

###make column with sample time starting from the time of the first sample recorded after appearance of the array image
MV2_SMP_FILLER$SampleTime <- (MV2_SMP_FILLER$Time - MV2_SMP_FILLER$TimeOfFirstSample) 

###Create a column (TimeStamp) which labels each incrementing sample time as a time stamp number starting with one, and increasing sequentially by one 
MV2_SMP_FILLER <- ddply(MV2_SMP_FILLER, c("Subject", "Trial"), transform, TimeStamp=seq(1, (length(SampleTime)), by=1))

###Write file with timestamp number###
##DO NOT REWRITE FILE!!!##
write.csv(MV2_SMP_FILLER, file = "MV2_SMP_FILLER_TSTAMP_Kids.csv")

MV2K_SMP = rbind(MV2_SMP_EXP, MV2_SMP_FILLER)
colnames(MV2_SMP_EXP)
colnames(MV2_SMP_FILLER)
write.csv(MV2K_SMP, file = "MV2K_EyetrackingR_Kids.csv")

MV2K_SMP <- read.csv(file = "MV2K_EyetrackingR_Kids.csv", header = TRUE)

#  Next File is MV2Adults_eyetrackingR_PGv2.r