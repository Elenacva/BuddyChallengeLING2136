###LCBS LAB### Descriptives and LME analysis of behavioural data 

##(C) Daniel Kleinman and Matteo Mascelloni, 2020


######################################################################

##loading packages####

library("lme4",) 
library("lmerTest")
library("plyr")
library("MASS")
library("trimr")
library('emmeans')
library('car')
library('multcomp')

emm_options(lmer.df="satterthwaite")
emm_options(lmerTest.limit = 100000)

####FUNCTIONS#####

# define function for subsequent use
find.random.slopes.to.remove <- function(fitted.lmer, propVarThreshold, override=0) {
  # this function takes a fitted LMER and, for each random factor, identifies the
  # number of random slopes n that account for less variance than a specified threshold
  # (typically .01) using the rePCA() function. then, it identifies the n random slopes
  # that allegedly account for the lowest variance. it outputs this information so it can
  # be easily copied/pasted into a subsequent LMER call.
  # 
  # NOTE: **This function does not work as intended if the model includes random correlations!!**
  #       It has also not been tested for optimizers other than 'bobyqa'.
  
  # check to make sure there are no random correlations
  variance.groups <- as.data.frame(summary(fitted.lmer)[["varcor"]])[["grp"]]
  
  rand.corrs.absent = length(variance.groups) == length(unique(variance.groups))
  
  if( rand.corrs.absent | override ) {
    
    # is.na() selects vars from the table as opposed to correlations
    fitted.lmer.varCor <- subset(as.data.frame(summary(fitted.lmer)[["varcor"]]), is.na(var2))
    
    rePCA_summary <- summary(rePCA(fitted.lmer))
    factorNames <- names(rePCA_summary)
    
    allToPrint <- c()
    
    for (i in 1:length(factorNames)) {
      nextFactor <- factorNames[i]
      nextPropVars <- rePCA_summary[[nextFactor]][["importance"]][2,]
      numSlopesToRemove <- sum(nextPropVars < propVarThreshold)
      if( numSlopesToRemove > 0 ) {
        rowNumsOfFactor <- which(grepl(paste(nextFactor, "\\.*", sep=""), fitted.lmer.varCor[["grp"]]))
        factorVarCors <- fitted.lmer.varCor[rowNumsOfFactor,]
        rowNumsOfSlopesToRemove <- order(factorVarCors[,"vcov"])[seq(numSlopesToRemove)]
        nameSlopesToRemove <- factorVarCors[rowNumsOfSlopesToRemove,"var1"]
        nameSlopesToRemove[nameSlopesToRemove == "(Intercept)"] <- "1"   # replace "(Intercept)" if present
        factorsToSubtract <- paste(nameSlopesToRemove, collapse=' + ')
      } else {
        factorsToSubtract <- "remove nothing"
      }
      
      nextToPrint <- paste(nextFactor, ": ", numSlopesToRemove, " / ", length(nextPropVars), " : (", factorsToSubtract, ")", sep="")
      allToPrint <- c(allToPrint, nextToPrint)
    }
    
    if( !rand.corrs.absent ) {
      cat(paste("\n*** This function does not work properly (and thus does not return standard output)",
                "\n*** when random correlations are present, as they are here. However, the override",
                "\n*** parameter was manually set, so it is returning output anyway",
                sep="")
      )
    }
    cat(paste("\n", paste(allToPrint, collapse="\n\n"), "\n\n", sep=""))
  } else {
    
    # identify correlated random slopes
    rand.comps.with.corr.slopes <- names(summary(as.factor(variance.groups))[summary(as.factor(variance.groups)) != 1])
    corr.slopes <- subset(as.data.frame(summary(fitted.lmer)[["varcor"]])[variance.groups %in% rand.comps.with.corr.slopes,], is.na(var2))
    rand.comps.with.corr.slopes <- unique(corr.slopes[["grp"]])   # same as above, but now in the same order as the function call
    
    # construct output to identify random slopes to user (and align spaces for easy reading)
    corr.slopes.output <- vector(mode="character", length=length(rand.comps.with.corr.slopes))
    max.uniqueComp.name.length <- max(sapply(rand.comps.with.corr.slopes, nchar))
    for( nextUniqueCompIndex in 1:length(rand.comps.with.corr.slopes) ) {
      nextUniqueComp <- rand.comps.with.corr.slopes[nextUniqueCompIndex]
      extra.spaces <- paste(rep(" ",max.uniqueComp.name.length - nchar(nextUniqueComp)), collapse="")
      corr.slopes.output[nextUniqueCompIndex] <- paste("   ", nextUniqueComp, ":   ", extra.spaces, paste(subset(corr.slopes, grp==nextUniqueComp)[["var1"]], collapse=", "), sep="")
    }
    
    output.text <- paste(
      "",
      "It appears this (g)lmer object has correlations between random slopes.",
      "The method used by this function to identify the slopes accounting for",
      "the least variance does not work properly when those correlations are",
      "present. (Note that this includes correlations between different levels",
      "of a single variable when it is specified as a nominal variable in the",
      "random effects; such variables must be explicitly recoded as numeric",
      "variables.) If you believe this is a mistake, rerun the function while",
      "setting optional parameter override=1...",
      "",
      "The correlated random slopes identified are shown below:",
      paste(corr.slopes.output, collapse="\n"),
      sep="\n")
    
    cat(output.text)
  }
}
#########################################################################################################

setwd("./Data_Analysis/")

##### Open the data files #######
##open csv file

exp1_data = data = read.csv('data_auditory.csv', stringsAsFactors = TRUE)


## set subj number

exp1_n=20

## remove errors and empty cell
exp1_data_clean = subset(exp1_data, exp1_data$accuracy=="1" )
exp1_data_clean = subset(exp1_data_clean, exp1_data_clean$rt!='-99')


## trim data 3 SD by subject, min RT 250ms

##Raw values
exp1_trimmedData <- sdTrim(data = exp1_data_clean, minRT = 250, sd = 3, 
                      perCondition = FALSE, perParticipant = TRUE, 
                      returnType = "raw", digits = 2)



exp2_data = data = read.csv('data_written.csv', stringsAsFactors = TRUE)

## set subj number

exp2_n=20

## remove errors and empty cell
exp2_data_clean = subset(exp2_data, exp2_data$accuracy=="1" )
exp2_data_clean = subset(exp2_data_clean, exp2_data_clean$rt!='NaN')


## trim data 3 SD by subject, min RT 250ms

##Raw values
exp2_trimmedData <- sdTrim(data = exp2_data_clean, minRT = 250, sd = 3, 
                      perCondition = FALSE, perParticipant = TRUE, 
                      returnType = "raw", digits = 2)



###### merge datasets #####

# prep Exp. 1 dataset for merging
exp1_trimmedData$experiment = 'exp1'
exp1_col_num_to_rename = which(names(exp1_trimmedData) == "sound")
names(exp1_trimmedData) = c(names(exp1_trimmedData)[1:(exp1_col_num_to_rename-1)], "distractor_word", names(exp1_trimmedData)[(exp1_col_num_to_rename+1):length(names(exp1_trimmedData))])
exp1_trimmedData$distractor_modality = 'auditory'
exp1_trimmedData$participant = paste(exp1_trimmedData$experiment, '_s', exp1_trimmedData$participant, sep="")

# prep Exp. 2 dataset for merging
exp2_trimmedData$experiment = 'exp2'
exp2_col_num_to_rename = which(names(exp2_trimmedData) == "word")
names(exp2_trimmedData) = c(names(exp2_trimmedData)[1:(exp2_col_num_to_rename-1)], "distractor_word", names(exp2_trimmedData)[(exp2_col_num_to_rename+1):length(names(exp2_trimmedData))])
exp2_trimmedData$distractor_modality = 'visual'
exp2_trimmedData$participant = paste(exp2_trimmedData$experiment, '_s', exp2_trimmedData$participant, sep="")
exp2_trimmedData$SOA = 0

# merge datasets and factorize nominal variables
allExp_trimmedData = rbind(  
    exp1_trimmedData[,c("experiment","participant","distractor_modality","condition","SOA","pic","distractor_word","rt")]
  , exp2_trimmedData[,c("experiment","participant","distractor_modality","condition","SOA","pic","distractor_word","rt")]
)
cols_to_factorize = c("experiment", "participant", "distractor_modality", "condition", "SOA", "pic", "distractor_word")
allExp_trimmedData[cols_to_factorize] = lapply(allExp_trimmedData[cols_to_factorize], factor)

## prepare model variables: distractor_modality, distractor_type, and distractor_relatedness. for each one:
#     (1) create the variable 
#     (2) order factor levels
#     (3) set contrasts
#     (4) create an explicitly numeric version of the variable to use in random effects

# distractor_modality 
allExp_trimmedData[["distractor_modality"]]               <- factor(allExp_trimmedData[["distractor_modality"]], c("auditory", "visual"))
contrasts(allExp_trimmedData[["distractor_modality"]])    <- -contr.sum(2)/2
allExp_trimmedData[["distractor_modality.num"]]           <- model.matrix(~ distractor_modality, data=allExp_trimmedData)[,2]

# distractor_type
allExp_trimmedData[["distractor_type"]]                   <- as.factor(as.character(ifelse(allExp_trimmedData[["condition"]] %in% c("category", "unrelated_cat"), "category", ifelse(allExp_trimmedData[["condition"]] %in% c("mediated", "unrelated_med"), "mediated", "???"))))
allExp_trimmedData[["distractor_type"]]                   <- factor(allExp_trimmedData[["distractor_type"]], c("category", "mediated"))
contrasts(allExp_trimmedData[["distractor_type"]])        <- -contr.sum(2)/2
allExp_trimmedData[["distractor_type.num"]]               <- model.matrix(~ distractor_type, data=allExp_trimmedData)[,2]

# distractor_relatedness
allExp_trimmedData[["distractor_relatedness"]]            <- as.factor(as.character(ifelse(allExp_trimmedData[["condition"]] %in% c("category", "mediated"), "related", ifelse(allExp_trimmedData[["condition"]] %in% c("unrelated_cat", "unrelated_med"), "unrelated", "???"))))
allExp_trimmedData[["distractor_relatedness"]]            <- factor(allExp_trimmedData[["distractor_relatedness"]], c("unrelated", "related"))
contrasts(allExp_trimmedData[["distractor_relatedness"]]) <- -contr.sum(2)/2
allExp_trimmedData[["distractor_relatedness.num"]]        <- model.matrix(~ distractor_relatedness, data=allExp_trimmedData)[,2]


# datasets for analysis of individual experiments
exp1_trimmedData <- subset(allExp_trimmedData, experiment=='exp1')
exp2_trimmedData <- subset(allExp_trimmedData, experiment=='exp2')



#################################### EXPERIMENT 1 (auditory distractors) ANALYSIS  ####################################
# Maximal model -> Do NOT converge
lmer.exp1 <- lmer(rt ~ distractor_type * distractor_relatedness + (1 + distractor_type.num * distractor_relatedness.num | participant) + (1 + distractor_type.num * distractor_relatedness.num | pic) + (1 + distractor_relatedness.num | distractor_word), data=exp1_trimmedData, control=lmerControl(optimizer="bobyqa"), verbose=2)

# Maximal model - No random correlation -> Do NOT converge
lmer.exp1.noRanCorrs <- lmer(rt ~ distractor_type * distractor_relatedness + (1 + distractor_type.num * distractor_relatedness.num || participant) + (1 + distractor_type.num * distractor_relatedness.num || pic) + (1 + distractor_relatedness.num || distractor_word), data=exp1_trimmedData, control=lmerControl(optimizer="bobyqa"), verbose=2)

# remove all random slopes except for one
find.random.slopes.to.remove(lmer.exp1.noRanCorrs, .01)


# Reduced model -> Does Converge
lmer.exp1.noRanCorrs.reduced <- lmer(rt ~ distractor_type * distractor_relatedness + (1 + distractor_relatedness.num | participant) + (1 | pic) + (1 | distractor_word), exp1_trimmedData, control=lmerControl(optimizer="bobyqa"), verbose=2)

lmer.exp1.noRanCorrs.reduced.summary <- summary(lmer.exp1.noRanCorrs.reduced)
#                                         Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                               754.064     23.167   26.418  32.548  < 2e-16 ***
# distractor_type1                           10.716      9.277   56.105   1.155  0.25296    
# distractor_relatedness1                    20.345      5.297   18.784   3.841  0.00112 ** 
# distractor_type1:distractor_relatedness1   -7.491      9.445 2623.554  -0.793  0.42778    

#  emmeans contrasts
lmer.exp1.noRanCorrs.reduced.emmeans <- emmeans(lmer.exp1.noRanCorrs.reduced, ~ distractor_relatedness + distractor_type)

lmer.exp1.noRanCorrs.reduced.emmeans.contrasts <- contrast(lmer.exp1.noRanCorrs.reduced.emmeans,
    as.data.frame(cbind(
      # Relatedness (U=Unrelated, R=Related ):      U   R   U   R
      # Type        (C=Category,  M=Mediated):      C   C   M   M
      category_interference                    = c(-1,  1,  0,  0)
    , mediated_interference                    = c( 0,  0, -1,  1)
  )))
# contrast              estimate   SE   df t.ratio p.value
# category_interference     24.1 7.07 59.4 3.406   0.0012 
# mediated_interference     16.6 7.12 61.0 2.331   0.0230 

#################################### EXPERIMENT 2 (written distractors) ANALYSIS  ####################################

# Maximal model -> Do NOT converge
lmer.exp2 <- lmer(rt ~ distractor_type * distractor_relatedness + (1 + distractor_type.num * distractor_relatedness.num | participant) + (1 + distractor_type.num * distractor_relatedness.num | pic) + (1 + distractor_relatedness.num | distractor_word), data=exp2_trimmedData, control=lmerControl(optimizer="bobyqa"), verbose=2)

#  Maximal model - No random correlation -> Do NOT converge
lmer.exp2.noRanCorrs <- lmer(rt ~ distractor_type * distractor_relatedness + (1 + distractor_type.num * distractor_relatedness.num || participant) + (1 + distractor_type.num * distractor_relatedness.num || pic) + (1 + distractor_relatedness.num || distractor_word), data=exp2_trimmedData, control=lmerControl(optimizer="bobyqa"), verbose=2)

# remove most random slopes and a random intercept 
find.random.slopes.to.remove(lmer.exp2.noRanCorrs, .01)

# note that the function above (find.random.slopes.to.remove) identifies which random (intercepts and) slopes account for < 1% of the variance
# *of their random factor*. 
summary(lmer.exp2.noRanCorrs)$varcor

# inspection of the summary showed that for distractor_word, both the intercept and the slope account for virtually no variance; therefore, both have been removed.

# Reduced model -> converged after 95 iterations
lmer.exp2.noRanCorrs.reduced <- lmer(rt ~ distractor_type * distractor_relatedness + (1 + distractor_type.num || participant) + (1 | pic), data=exp2_trimmedData, control=lmerControl(optimizer="bobyqa"), verbose=2)

lmer.exp2.noRanCorrs.reduced.summary <- summary(lmer.exp2.noRanCorrs.reduced)
#                                         Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                               654.220     21.367   24.388  30.618  < 2e-16 ***
# distractor_type1                           -6.892      5.379   18.968  -1.281  0.21557    
# distractor_relatedness1                     6.009      4.501 2649.752   1.335  0.18204    
# distractor_type1:distractor_relatedness1  -24.970      9.001 2648.959  -2.774  0.00557 ** 

# emmeans contrasts
lmer.exp2.noRanCorrs.reduced.emmeans <- emmeans(lmer.exp2.noRanCorrs.reduced, ~ distractor_relatedness + distractor_type)

lmer.exp2.noRanCorrs.reduced.emmeans.contrasts <- contrast(lmer.exp2.noRanCorrs.reduced.emmeans,
    as.data.frame(cbind(
      # Relatedness (U=Unrelated, R=Related ):      U   R   U   R
      # Type        (C=Category,  M=Mediated):      C   C   M   M
      category_interference                    = c(-1,  1,  0,  0)
    , mediated_interference                    = c( 0,  0, -1,  1)
  )))
# contrast              estimate   SE   df t.ratio p.value
# category_interference    18.49 6.39 2650  2.895  0.0038 
# mediated_interference    -6.48 6.34 2649 -1.021  0.3072 


#################################### COMBINED ANALYSIS  ####################################



#    -maximal random effects:
#        participant:            type * relatedness   [not modality, since each participant only did PWI with distractors of a single modality]
#        pic:         modality * type * relatedness   [pictures were paired with all kinds of distractors in both experiments]
#        word:        modality        * relatedness   [not type, since an individual word was either part of the category-related OR mediated-related set but not both]
#    -bobyqa optimizer 
#    -using nominal variables (with contrasts set appropriately) for fixed effects, and numeric variables for random effects


# Maximal Model -> do NOT coverge
lmer.crossExp            <- lmer(rt ~ distractor_modality * distractor_type * distractor_relatedness + (1 + distractor_type.num * distractor_relatedness.num |  participant) + (1 + distractor_modality.num * distractor_type.num * distractor_relatedness.num |  pic) + (1 + distractor_modality.num * distractor_relatedness.num |  distractor_word), data=allExp_trimmedData, control=lmerControl(optimizer="bobyqa"), verbose=2)

# Maximal Model -  No Random correlations -> do NOT coverge
lmer.crossExp.noRanCorrs <- lmer(rt ~ distractor_modality * distractor_type * distractor_relatedness + (1 + distractor_type.num * distractor_relatedness.num || participant) + (1 + distractor_modality.num * distractor_type.num * distractor_relatedness.num || pic) + (1 + distractor_modality.num * distractor_relatedness.num || distractor_word), data=allExp_trimmedData, control=lmerControl(optimizer="bobyqa"), verbose=2)

# remove 10/16 random slopes
find.random.slopes.to.remove(lmer.crossExp.noRanCorrs, .01)


# Reduced model -> converged after 284 iterations
lmer.crossExp.noRanCorrs.reduced <- lmer(rt ~ distractor_modality * distractor_type * distractor_relatedness + (1 | participant) + (1 + distractor_modality.num + distractor_relatedness.num || pic) + (1 + distractor_modality.num || distractor_word), data=allExp_trimmedData, control=lmerControl(optimizer="bobyqa"), verbose=2)

lmer.crossExp.noRanCorrs.reduced.summary <- summary(lmer.crossExp.noRanCorrs.reduced)
#                                                              Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                                                    704.150     16.650   58.409  42.292  < 2e-16 ***
# distractor_modality1                                           -99.797     29.551   39.972  -3.377  0.00164 ** 
# distractor_type1                                                 1.901      4.713   37.381   0.403  0.68895    
# distractor_relatedness1                                         13.010      3.850   15.616   3.379  0.00394 ** 
# distractor_modality1:distractor_type1                          -17.524     11.072   57.271  -1.583  0.11899    
# distractor_modality1:distractor_relatedness1                   -14.352      6.538 5268.056  -2.195  0.02819 *  
# distractor_type1:distractor_relatedness1                       -16.404      6.538 5268.729  -2.509  0.01213 *  
# distractor_modality1:distractor_type1:distractor_relatedness1  -17.165     13.076 5266.933  -1.313  0.18936    

# emmeans contrasts
lmer.crossExp.noRanCorrs.reduced.emmeans <- emmeans(lmer.crossExp.noRanCorrs.reduced, ~  distractor_relatedness + distractor_type + distractor_modality)
lmer.crossExp.noRanCorrs.reduced.contrastsDf <- as.data.frame(cbind(
    # Relatedness (U=Unrelated, R=Related ):     U     R     U     R     U     R     U     R
    # Type        (C=Category,  M=Mediated):     C     C     M     M     C     C     M     M
    # Modality    (A=Auditory,  W=Visual  ):     A     A     A     A     V     V     V     V
    crossModality_category_interference    = c(1/2, -1/2,    0,    0, -1/2,  1/2,    0,    0)
  , crossModality_mediated_interference    = c(  0,    0,  1/2, -1/2,    0,    0, -1/2,  1/2)
))


lmer.crossExp.noRanCorrs.reduced.emmeans.contrasts <- contrast(lmer.crossExp.noRanCorrs.reduced.emmeans, lmer.crossExp.noRanCorrs.reduced.contrastsDf)
# contrast                            estimate   SE   df t.ratio p.value
# crossModality_category_interference    -2.89 4.62 5268 -0.624  0.5327 
# crossModality_mediated_interference   -11.47 4.62 5267 -2.481  0.0131 
