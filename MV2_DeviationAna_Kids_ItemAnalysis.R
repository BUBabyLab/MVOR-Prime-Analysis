##  KAYDEN EDITS ITEM ANALYSIS Divergence analysis from EyeTrackingR
#  PG - Oct, 2018
#

library(gtools) 
library(lme4)
library(plyr)
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
library(reshape2)
getwd()

#source ("eyeRTrackLossCreateColumn.r")
source("../MVOR-Scripts-fun.r") # PG-defined functions for script

#InReadData <- read.csv(file = "MV2K_EyetrackingR_Kids.csv", header = TRUE)
#Read in files for trimmed (based on looks to prime) vs not analysis June 2019
InReadData <- read.csv(file = "MV2_Array_SMP_EXP_ONLY_TSTAMP_Kids.csv", header = TRUE)
InReadData_Trimmed <- read.csv(file = "MV2_Array_SMP_EXP_ONLY_TSTAMP_Kids_Trimmed.csv", header = TRUE)
unique(InReadData_Trimmed$Subject)

unique(InReadData$Subject)
nrow(InReadData)
#Remove subjects with overall too much data loss (from prior analysis)
# See file "MV2A_Kids_ALL_FALSE_TracklossNoNonAOI" for explanation
InReadData_Trimmed <- subset(InReadData_Trimmed, Subject != 1629059)
InReadData_Trimmed <- subset(InReadData_Trimmed, Subject != 1629061)
InReadData_Trimmed <- subset(InReadData_Trimmed, Subject != 1629062)
unique(InReadData_Trimmed$Subject)
nrow(InReadData)

MV2K <- InReadData_Trimmed

#Eliminate 2 objects:  Saw2 and Camera2 (Kids did not recognize well)
MV2K <- subset(MV2K, TargetFile != "fsaw_2_1.png")
MV2K <- subset(MV2K, TargetFile != "fsaw_2_2.png") 
MV2K <- subset(MV2K, TargetFile != "psaw_2_1.png")
MV2K <- subset(MV2K, TargetFile != "psaw_2_2.png") 

MV2K <- subset(MV2K, TargetFile != "fcamera_2_1.png")
MV2K <- subset(MV2K, TargetFile != "fcamera_2_2.png")
MV2K <- subset(MV2K, TargetFile != "pcamera_2_1.png") 
MV2K <- subset(MV2K, TargetFile != "pcamera_2_2.png")
nrow(MV2K)

#Filter to only have the image of interest in dataframe
MV2K_fairplane <- subset(MV2K, TargetFile == "fairplane_1_1.png")
MV2K_pairplane <- subset(MV2K, TargetFile == "pairplane_1_1.png")

#  Uncomment another age group to test different ages
ageLabel <- "5 Yr-Olds"
yngKids = filter (MV2K, AgeYear == 3)
middleKids = filter (MV2K, AgeYear == 4)
oldKids_p = filter (MV2K_pairplane, AgeYear == 5)
oldKids_f = filter (MV2K_fairplane, AgeYear == 5)

younger = yngKids = filter (MV2K, AgeYear != 5)
older = filter (MV2K, AgeYear != 3)

####create file for eyetrackingR
MV2K_ETR_KIDS_FALSE <- make_eyetrackingr_data(oldKids_p, 
                                              participant_column = "Subject",
                                              trial_column = "Trial",
                                              time_column = "SampleTime",
                                              trackloss_column = "SampleLost",
                                              aoi_columns = c('Target','Complement', 'DiffExemplar', 'Distractor'),
                                              treat_non_aoi_looks_as_missing = FALSE)

###########################################################################################
#######Look at track loss, and remove participants and trials with too much track loss#####
###########################################################################################
#Create a file removing any time points after 1600ms 
MV2K_ETR_1600 <- subset_by_window(MV2K_ETR_KIDS_FALSE,
                                  window_start_time = 0, 
                                  window_end_time = 1900, 
                                  rezero = TRUE, 
                                  remove = TRUE)

#How many subjects
unique(MV2K_ETR_1600[, c('Subject')])

###Look at track loss data (data where eyetracker lost the eye)
trackloss <- trackloss_analysis(data = MV2K_ETR_1600)
trackloss
write.csv(trackloss, file = "MV2A_Kids_ALL_FALSE_FULL_CHECK_TracklossNoNonAOI.csv")
#trackloss <- read.csv("MV2A_Adults_TracklossNoNonAOI.csv")

#Graph trackloss by subject
meanTLoss_by_Subj <- trackloss %>% 
  group_by(Subject) %>% 
  summarize(averaged.TLoss = mean(TracklossForParticipant))

#Graph samples remaining/present per subject
MeanSamplesPresent  <- trackloss %>% 
  group_by(Subject) %>% 
  summarize(averaged.Samples = mean(Samples))

describe(MeanSamplesPresent)
describe(meanTLoss_by_Subj)

nrow(MV2K_ETR_1600)
###Remove subjects and trials with over a proportion certain proportion of track loss
# NOTE:  .95 because initial run with .25 resulted in removal of three subjects;
#  (see top of file)
unique(MV2K_ETR_1600$Subject)
MV2K_Clean <- clean_by_trackloss(MV2K_ETR_1600, participant_prop_thresh = 0.25, 
                                 trial_prop_thresh = 0.25)
nrow(MV2K_Clean)
unique(MV2K_Clean$Subject)
#  YES!  You are doing this intentionally!
#  Filter with treat_non_aoi_looks_as_missing = FALSE first,
#  Remove the 198 trials that do not meet 25% criterion (step just above)
#  Then refilter with treat_non_aoi_looks_as_missing = TRUE to remove non-AOI looks
#  Note!  You are using MV2K_Clean, not MV2K as input to the next call!
MV2K_ETR_KIDS <- make_eyetrackingr_data(MV2K_Clean, 
                                        participant_column = "Subject",
                                        trial_column = "Trial",
                                        time_column = "SampleTime",
                                        trackloss_column = "SampleLost",
                                        aoi_columns = c('Target','Complement', 'DiffExemplar', 'Distractor'),
                                        treat_non_aoi_looks_as_missing = TRUE)

MV2K_Clean <- clean_by_trackloss(MV2K_ETR_KIDS, participant_prop_thresh = 0.95, 
                                 trial_prop_thresh = .65)
nrow(MV2K_Clean)
########################################################
###     Descriptives and figure time!            #######
########################################################
#
# Create proportion of looks across two AOIs for each of the features, parts and 
#   target-only ("fillers") conditions, in time bins of 25 msec
#
TimeBins <- make_time_sequence_data(MV2K_Clean, 
                                    time_bin_size = 25, 
                                    predictor_columns = c("Condition", "TrialType"),
                                    aois = c("Target", "Complement", "DiffExemplar", "Distractor"),
                                    summarize_by = "Subject"
)
colnames(TimeBins)
TimeBins$AOI<-as.factor(TimeBins$AOI)

####Remove NaN's from proportion
TB_NoNA <- TimeBins[TimeBins$Prop != 'NaN', ]
#TB_NoNA

# How many NAs?
TB_NAs <- TimeBins[TimeBins$Prop == 'NaN', ]
#TB_NAs

###Create df with only Experimental Trials
TB_Exp <- TB_NoNA[TB_NoNA$TrialType == 'Experimental', ]

# Within that df, restric to Parts only
TB_Exp_Parts <- TB_Exp[TB_Exp$Condition == 'Parts', ]

# Within that df, restric to Parts only
TB_Exp_Feat <- TB_Exp[TB_Exp$Condition == 'Features', ]

###Create file with only Filler Trials
# And divide into parts & features
TB_Filler <- TB_NoNA[TB_NoNA$TrialType == 'Filler', ]
TB_Filler_Parts <- TB_Filler[TB_Filler$Condition == "Parts", ]
TB_Filler_Feat <- TB_Filler[TB_Filler$Condition == "Features", ]

######################
##  PARTS:
#####################

############################################
######## Target vs, complement analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Parts_TvsC <- TB_Exp_Parts[(TB_Exp_Parts$AOI == 'Target' | TB_Exp_Parts$AOI == 'Complement') , ]
myClusterTTest_Compare (LookProp_Exp_Parts_TvsC)

############################################
######## Target vs, DiffExempar analysis####
############################################

###Create file with only target and Diff Exemp
LookProp_Exp_Parts_TvsDE <- TB_Exp_Parts[(TB_Exp_Parts$AOI == 'Target' | TB_Exp_Parts$AOI == 'DiffExemplar') , ]

myClusterTTest_Compare (LookProp_Exp_Parts_TvsDE)
myCluster_LMER_Compare(LookProp_Exp_Parts_TvsDE)

############################################
######## Target vs. Distractor analysis#####
############################################
LookProp_Exp_Parts_TvsDist <- TB_Exp_Parts[(TB_Exp_Parts$AOI == 'Target' | TB_Exp_Parts$AOI == 'Distractor') , ]
myClusterTTest_Compare (LookProp_Exp_Parts_TvsDist)
#myClusterCompare (LookProp_Exp_Parts_TvsDist)

############################################
######## Complement vs, DiffExempar analysis####
############################################

###Create file with only target and Diff Exemp
LookProp_Exp_Parts_CompvsDE <- TB_Exp_Parts[(TB_Exp_Parts$AOI == 'Complement' | TB_Exp_Parts$AOI == 'DiffExemplar') , ]

myClusterTTest_Compare (LookProp_Exp_Parts_CompvsDE)
myCluster_LMER_Compare(LookProp_Exp_Parts_CompvsDE)

#########################
###  FEATURES
########################

############################################
######## Target vs. complement analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Feat_TvsC <- TB_Exp_Feat[(TB_Exp_Feat$AOI == 'Target' | TB_Exp_Feat$AOI == 'Complement') , ]
myClusterTTest_Compare (LookProp_Exp_Feat_TvsC)
myCluster_LMER_Compare(LookProp_Exp_Feat_TvsC)

#LookProp_Exp_Feat_TvsC$TargetRef<- relevel(LookProp_Exp_Feat_TvsC$AOI, ref = "Target")

############################################
######## Target vs, DiffExempar analysis#####
############################################

###Create file with only target and complement
LookProp_Exp_Feat_TvsDE <- TB_Exp_Feat[(TB_Exp_Feat$AOI == 'Target' | TB_Exp_Feat$AOI == 'DiffExemplar') , ]
myClusterTTest_Compare (LookProp_Exp_Feat_TvsDE)

############################################
######## Target vs. Distractor analysis#####
############################################

LookProp_Exp_Feat_TvsDist <- TB_Exp_Feat[(TB_Exp_Feat$AOI == 'Target' | TB_Exp_Feat$AOI == 'Distractor') , ]
myClusterTTest_Compare (LookProp_Exp_Feat_TvsDist)
myCluster_LMER_Compare(LookProp_Exp_Feat_TvsDist)

############################################
######## Complement vs. Diff Exemplar  #####
############################################

LookProp_Exp_CompvsDiffExem <- TB_Exp_Feat[(TB_Exp_Feat$AOI == 'Complement' | TB_Exp_Feat$AOI == 'DiffExemplar') , ]
myClusterTTest_Compare (LookProp_Exp_CompvsDiffExem)
#####################
#FixProp <- ddply(TB_Exp, c("Condition", "AOI", "TimeBin"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
#FixProp
nrow(FixProp)

###################################
###Parts x AOI Figure##############
###################################

source("MVOR-Scripts-fun.r") # PG-defined functions for script
FixProp <- ddply(TB_Exp, c("Condition", "AOI", "Time"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
FixProp_Parts <- FixProp[FixProp$Condition == 'Parts', ]

###Reorder bars###
#  Set for 1900 msec
FixProp_Parts <- FixProp_Parts[c(115:152, 1:38, 39:76, 77:114), ]
FixProp_Parts$AOI <- factor(FixProp_Parts$AOI, levels = c("Target", "Complement", "DiffExemplar", "Distractor"))
myGGplot(FixProp_Parts, "5-yr-olds Parts Priming pAirplane 1_1 Trimmed")

###################################
###  Features x AOI Figure  #######
###################################

source("MVOR-Scripts-fun.r") # PG-defined functions for script
FixProp <- ddply(TB_Exp, c("Condition", "AOI", "Time"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
FixProp_Feat <- FixProp[FixProp$Condition == 'Features', ]

###Reorder bars###
# 1600 msec cutoff
#FixProp_Feat <- FixProp_Feat[c(193:256, 1:64, 65:128, 129:192), ]
# 1900 msec cutoff
FixProp_Feat <- FixProp_Feat[c(229:304, 1:76, 77:152, 153:228), ]
FixProp_Feat$AOI <- factor(FixProp_Feat$AOI, levels = c("Target", "Complement", "DiffExemplar", "Distractor"))
myGGplot(FixProp_Feat, "5-Yr-Olds Feature Priming fAirplane 1_1")

#######################################################
##  FILLER (TARGET ONLY) TRIALS:
## This is for Parts and Features! Separate trials
######################################################
Test_Filler_Parts <- TB_Filler_Parts[(TB_Filler_Parts$AOI == 'Target' | TB_Filler_Parts$AOI == 'DiffExemplar') , ]
Test_Filler_Feat <- TB_Filler_Feat[(TB_Filler_Feat$AOI == 'Target' | TB_Filler_Feat$AOI == 'DiffExemplar') , ]

# This step averages across the 3 distractors so that the T-test can look at only one distractor level
# Does not quite work!  PG 12/14/18 - need to actually transform to average this...
TB_Filler$AOI <- ifelse((TB_Filler$AOI == 'Target'), 'Target', 'Distractor')
# This should combine the 3 distractors into one but it does not
TB_Filler_Ave <- ddply(TB_Filler, c("Condition", "AOI", "TimeBin"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))

myClusterTTest_Compare (Test_Filler_Parts)
myClusterTTest_Compare (Test_Filler_Feat)

#Re-code to compare Parts to Features across trials (Target-Only trials)
TB_Target_Only <- TB_Filler[TB_Filler$AOI == "Target", ]
TB_Filler$AOI_Comp <- ifelse((TB_Filler$Condition == 'Parts'), 'Target', 'Distractor')

###Make Variable so I can average across all the distractors
TB_Filler$AOI_Filler <- ifelse((TB_Filler$AOI == 'Target'), 'Target', 'Distractor')

###Fixation proportion###
#FixProp_Filler

############################################
# Compare Parts and Features on fix prop
############################################
Test_Filler_AllTarg <- TB_Filler[(TB_Filler$AOI == 'Target'), ]
myClusterTTest_CompareBtwnConditions (Test_Filler_AllTarg)

###Reorder bars###
FixProp_Filler <- ddply(TB_Filler, c("Condition", "AOI_Filler", "Time"), summarise, MeanProp=mean(Prop), SEProp=sd(Prop)/sqrt(length(Prop)))
FixProp_Filler <- FixProp_Filler[c(193:256, 65:128, 129:192, 1:64), ]

FixProp_Filler$ArrayObject[1:46] <- "PartsTarget"
FixProp_Filler$ArrayObject[65:128] <- "FeatTarget"
FixProp_Filler$ArrayObject[129:192] <- "PartsDistr"
FixProp_Filler$ArrayObject[193:256] <- "FeatDistr"

source("MVOR-Scripts-fun.r") # PG-defined functions for script
FixProp_Filler$AOI_Filler <- factor(FixProp_Filler$ArrayObject, levels = c("PartsTarget", "FeatTarget", "PartsDistr", "FeatDistr"))
myFillerGGplot(FixProp_Filler, "Kids Target Only")