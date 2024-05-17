# data_intake.R
#
# Utility functions to read in study data.

getPoC <- function() {
    PoC <- read_csv("..\\data\\20230223 PoC.csv",
        col_types = cols(
            Replicate = col_integer(),
            Experiment = col_integer()
        )
    )
    return(PoC)
}


getStudyData <- function() {
    # read in the "standard" data -- massaged in Excel to be well-posed
    easyData <- read_csv("..\\data\\20210810_normal_labs.csv")

    # read in the MassCPR data -- this is LONG already, as replicates were done in different batches
    MassCPR_l <- read_csv("..\\data\\20210603 MassCPR.csv")

    # read in the NIST and NML data --
    NMIdata <- read_csv("..\\data\\20210603 NIST-NML-NIB.csv")

    # make the 9-lab data frame tidy/long for analysis and plotting
    easyData_l <- pivot_longer(easyData, cols = c(9:12), names_to = "Rep", values_to = "Sig")

    # make the NMI data frame long
    NMIdata_l <- pivot_longer(NMIdata, cols = c(9:12), names_to = "Rep", values_to = "Sig")

    # bind in the rows of the all the data
    allData_l <- rbind(easyData_l, MassCPR_l, NMIdata_l)

    # order the rows of the data frame by the Lab, Target, and Sample Name
    allData_l <- arrange(allData_l, Lab, Target, SamName)

    # make factors out of Lab, SamName, Batch
    allData_l$Lab <- as.factor(allData_l$Lab)
    allData_l$SamName <- as.factor(allData_l$SamName)
    allData_l$Batch <- as.factor(allData_l$Batch)

    # create a concentration (in IU) column for the WHO IS dilutions
    allData_l$conc <- -1.11

    # To prepare the bulk material, quantification of the SARS-CoV-2 genome copies
    # within the inactivated material was determined relative to a plasmid standard
    # curve by in-house real-time RT-PCR using primer/probe targeting the E-gene [9].
    # The material was prepared to contain 1x10^8 genomes per mL and as with the
    # chimeric LVP, was formulated in universal buffer containing a background of
    # 1x10^5 copies/mL of human genomic DNA.
    #     from WHO/BS/2020.2403 Page 5

    # 3. UNITAGE
    # The assigned potency of the WHO International Standard for SARS-CoV- 2 RNA
    # for NAT-based assays is 7.40 Log10 IU/ampoule. After reconstitution in 0.5mL
    # of molecular grade water or PBS,
    # the final concentration of the preparation is 7.70 Log10 IU/mL.
    #       from
    # WHO International Standard
    # First WHO International Standard for SARS-CoV-2 RNA NIBSC code: 20/146
    # Instructions for use
    # (Version 2.0, Dated 05/01/2021)


    # fill in the concentrations (in IU) for the WHO IS dilutions
    allData_l[which(allData_l$SamName == "WHO-IS_dilution_0"), ]$conc <- 10^7.7
    allData_l[which(allData_l$SamName != "WHO-IS_dilution_0"), ]$conc <- NA
    allData_l[which(allData_l$SamName == "WHO-IS_dilution_1"), ]$conc <- 10^6.7
    allData_l[which(allData_l$SamName == "WHO-IS_dilution_2"), ]$conc <- 10^5.7
    allData_l[which(allData_l$SamName == "WHO-IS_dilution_3"), ]$conc <- 10^4.7
    allData_l[which(allData_l$SamName == "WHO-IS_dilution_4"), ]$conc <- 10^3.7
    allData_l[which(allData_l$SamName == "WHO-IS_dilution_5"), ]$conc <- 10^2.7
    allData_l[which(allData_l$SamName == "WHO-IS_dilution_6"), ]$conc <- 10^1.7

    # create a logical column to identify the calibrants
    allData_l$isCalib <- str_detect(as.character(allData_l$SamName), "WHO-IS", negate = FALSE)

    # special case MUSC data -- they rehydrated the WHO IS in 1.0 mL not 0.5
    allData_l[which(allData_l$Lab == "MUSC" & allData_l$isCalib), ]$conc <-
        allData_l[which(allData_l$Lab == "MUSC" & allData_l$isCalib), ]$conc / 2

    # create a log Conc column which we'll fit against as the independent variable
    allData_l$logConc <- log(allData_l$conc, 10)

    # create column for log Signal (needed for dPCR)
    allData_l$logSig <- -99

    # for dPCR, calculate and store log signal
    allData_l$logSig[which(allData_l$dPCR)] <- log(allData_l$Sig[which(allData_l$dPCR)], 2)

    # for anything that's NOT dPCR, make value NA
    allData_l$logSig[which(!allData_l$dPCR)] <- NA

    # for any -Inf make it NA
    allData_l$logSig[is.infinite(allData_l$logSig)] <- NA

    # create a logConcPredicted column for predictions
    allData_l$logConcPredicted <- log(-1, 10)
    allData_l$logConcPredicted[is.nan(allData_l$logConcPredicted)] <- NA

    # here remove control targets from logical column isCtrl
    allData_l <- allData_l[which(!allData_l$isCtrl), ]

    # let's relevel the SamName factor so it's got the Int Std dilution series first
    # drop the unused sample names that were controls
    allData_l$SamName <- fct_drop(factor(allData_l$SamName, levels = levels(allData_l$SamName)[c(9:15, 1:8, 16)]))

    # calculate weights for a robust fit
    # group by Lab:Target:SamName
    myGroupedAll <- group_by(allData_l, Lab, Target, SamName, Batch)

    # calculate the weights for dPCR, where the dependent variable is log Signal
    myGroupedAllWtsQpcr <- myGroupedAll[which(!myGroupedAll$dPCR), ] %>% mutate(wts = 1 / sd(Sig, na.rm = TRUE)^2)

    # calculate the weights for qPCR, where the dependent variable is Cq, which is already log space
    myGroupedAllWtsDpcr <- myGroupedAll[which(myGroupedAll$dPCR), ] %>% mutate(wts = 1 / sd(logSig, na.rm = TRUE)^2)

    # ungroup the data for further use...
    myGroupedAllWtsDpcr <- ungroup(myGroupedAllWtsDpcr)
    myGroupedAllWtsQpcr <- ungroup(myGroupedAllWtsQpcr)

    # make a composite data frame
    myDataWts <- rbind(myGroupedAllWtsDpcr, myGroupedAllWtsQpcr)

    # set the weights to non-calibrant samples to NA
    myDataWts[which(!myDataWts$isCalib), ]$wts <- NA

    # return this munged tibble
    return(myDataWts)
}
