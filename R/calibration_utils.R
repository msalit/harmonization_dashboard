# This function will analyze the data
#  • it will separate the long tibble (r data frame) into a grouped tibble
#  • it will calibrate the relationship between the IU and the signal
#  • differently for dPCR and qPCR
calAndPredict <- function(studyData, makePlot = TRUE) {
    # make an analyzable nested dataframe for each target within each lab -- for calibration
    myGroupedData <- group_by(studyData, Lab, Target)

    # fit and predict all the nested chunks
    calAndPredictModel <- myGroupedData %>%
        group_modify(~ fitGroup(myLabTarget = .x, gr = .y))

    # here modify the originally passed in data frame
    studyData <- merge(studyData, calAndPredictModel)

    return(calAndPredictModel)
}

# this is the function to fit a Lab:Target calibration curve
# written in a single spot at a point when it was in 3 places
# maintaining it as a dynamic function instead of a "calibrate 1x"
# so we can perhaps create dynamic calibrant selection checkboxes in Shiny
fitCal <- function(myLabTarget) {
    if (myLabTarget$dPCR[1]) { # fit dPCR model

        # test for more than 1 calibrant
        if ((length(levels(fct_drop(factor(myLabTarget$SamName[which(myLabTarget$isCalib & !is.na(myLabTarget$Sig))]))))) > 1) {
            # fit the model log sig ~ logConc
            myFit <- lm(logSig ~ logConc,
                data = myLabTarget[which(myLabTarget$isCalib), ],
                weights = myLabTarget[which(myLabTarget$isCalib), ]$wts,
                na.action = "na.omit"
            )
        } else { # only 1 calibrant

            # force calib through 0, no weights
            myFit <- lm(logSig ~ logConc - 1,
                data = myLabTarget[which(myLabTarget$isCalib), ],
                na.action = "na.omit"
            )
        }
    } else { # fit qPCR model

        # test for more than 1 calibrant
        if ((length(levels(fct_drop(factor(myLabTarget$SamName[which(myLabTarget$isCalib & !is.na(myLabTarget$Sig))]))))) > 1) {
            # fit
            myFit <- lm(Sig ~ logConc,
                data = myLabTarget[which(myLabTarget$isCalib), ],
                weights = myLabTarget[which(myLabTarget$isCalib), ]$wts,
                na.action = "na.omit"
            )
        } else # force calib through 0
        {
            # fit forced through 0, no weights
            myFit <- lm(Sig ~ logConc - 1,
                data = myLabTarget[which(myLabTarget$isCalib), ],
                na.action = "na.omit"
            )
        }
    }

    return(myFit)
}

# This function will get a grouped data element
#  that represents data from a particular lab and target
#  and calibrate, fit, and predict
fitGroup <- function(myLabTarget, gr = .y) {
    # do the fit
    myFit <- fitCal(myLabTarget)

    # now do predictions
    if (myLabTarget$dPCR[1]) {
        # calculate dPCR predictions

        # test for more than 1 calibrant
        if ((length(levels(fct_drop(
            factor(myLabTarget$SamName[which(myLabTarget$isCalib &
                !is.na(myLabTarget$Sig))])
        )))) > 1) {
            # do the prediction of logConc from Sig
            myLabTarget$logConcPredicted <- (myLabTarget$logSig - myFit$coefficients[1]) /
                myFit$coefficients[2]
        } else {
            # only 1 calibrant

            # predict from straight proportion, no intercept
            # this should be for NMIs that did only 1 level
            # disregard fit
            straightProportion <- mean(
                myLabTarget[which(myLabTarget$SamName == "WHO-IS_dilution_0"), ]$conc /
                    myLabTarget[which(myLabTarget$SamName == "WHO-IS_dilution_0"), ]$Sig,
                na.rm = TRUE
            )

            myLabTarget$logConcPredicted <- log(myLabTarget$Sig * straightProportion, 10)
        }
    } else {
        # fit qPCR model

        # test for more than 1 calibrant
        if ((length(levels(fct_drop(
            factor(myLabTarget$SamName[which(myLabTarget$isCalib &
                !is.na(myLabTarget$Sig))])
        )))) > 1) {
            # do the prediction of logConc from Sig
            myLabTarget$logConcPredicted <- (myLabTarget$Sig - myFit$coefficients[1]) /
                myFit$coefficients[2]
        } else {
            # predict from straight proportion, no intercept
            myLabTarget$logConcPredicted <- myLabTarget$Sig / myFit$coefficients[1]
        }
    }

    # correct for dilutions, multiply predicted values by Dilution column
    myLabTarget$logConcPredicted <- log10((10^myLabTarget$logConcPredicted) / myLabTarget$Dilution)

    return(myLabTarget)
}



# calculate but don't plot
calcResults <- function(studyData) {
    # make an analyzeable nested dataframe for each target within each lab -- for calibration
    myGroupedData <- group_by(studyData, Lab, Target)

    # and fit all the nested chunks
    myResults <- myGroupedData %>%
        group_modify(~ fitAndPredict(data = .x, gr = .y))

    return(myResults)
}
