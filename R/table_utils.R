# table_utils.R contains functions which create tables from the study data


# this function will return a tibble of a results table
materialTable <- function(myResults) {
    # make a grouped tibble of the results
    myGrpResults <- group_by(myResults[which(!myResults$isCalib), ], Lab, Target, SamName)
    samVals <- myGrpResults %>%
        summarise(meds = round(median(logConcPredicted, na.rm = TRUE), 2)) %>%
        arrange(Lab)
    myTable <- pivot_wider(ungroup(samVals), names_from = c(SamName), values_from = meds)
    return(myTable)
}

# This function will return a table of medians and 2x mads (95% CI) for each material
materialSummary <- function(myResults) {
    # make a grouped tibble of the results
    myGrpResults <- group_by(myResults[which(!myResults$isCalib), ], SamName, Lab, Target)

    samVals <- myGrpResults %>%
        summarise(meds = median(logConcPredicted, na.rm = TRUE)) %>%
        arrange(Lab)

    samVals <- group_by(ungroup(samVals), SamName)

    samMedian <- samVals %>%
        summarise(
            gmeds = median(meds, na.rm = TRUE),
            mads = mad(meds, na.rm = TRUE) * 2
        ) %>%
        arrange(SamName)

    samMedian <- ungroup(samMedian)
    samMedian$gmeds <- signif(samMedian$gmeds, 4)
    samMedian$mads <- signif(samMedian$mads, 3)
    return(samMedian)
}

# this function will return a tibble of a full Lab:Target table
labTargTable <- function(myLT) {
    # reorder the table columns for presentation
    myRawBrowser <- data.frame(
        Sample = myLT$SamName,
        Signal = myLT$Sig,
        Wts = myLT$wts,
        Dilution = myLT$Dilution,
        log10_IU_per_mL_pred = myLT$logConcPredicted,
        Replicate = myLT$Rep,
        Batch = myLT$Batch,
        isCalib = myLT$isCalib,
        isDpcr = myLT$dPCR
    )
    # return( na.omit( myRawBrowser ) )
    return(myRawBrowser)
}
