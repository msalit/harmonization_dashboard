# plotting_utils.R contains functions which plot the study data in various methods for visualization.

# plotCal plots a calibration curve and the data around it,
# given data for a specific lab and target
plotCal <- function(myData) {
    # get the plot title
    myTitle <- paste(myData$Lab, myData$Target)

    # different plots for dPCR and qPCR
    if (myData$dPCR[1]) {
        # graph the dPCR fit and the predictions
        myCalPlot <- ggplot(
            myData[which(myData$isCalib), ],
            aes(x = 10^logConc, y = Sig)
        )
        p <- myCalPlot + geom_point() +
            geom_smooth(method = "lm", formula = "y ~ x", mapping = aes(weight = wts)) +
            labs(title = myTitle, subtitle = "weighted fit") +
            scale_x_log10(
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))
            ) +
            scale_y_continuous(
                trans = "log2",
                breaks = scales::trans_breaks("log2", function(x) 2^x),
                labels = scales::trans_format("log2", scales::math_format(2^.x))
            ) +
            ylab("dPCR Count") +
            xlab("Concentration in IU/mL") +
            theme_bw() +
            theme(
                panel.grid.minor = element_blank()
            ) +
            labs(color = "Sample") +
            geom_rug(data = myData[which(!myData$isCalib), ], aes(x = 10^logConcPredicted, y = Sig, color = SamName))

        suppressWarnings(print(p))
    } else { # plot the qPCR calibration

        # graph the calibration curve and predictions
        myCalPlot <- ggplot(
            myData[which(myData$isCalib), ],
            aes(x = 10^logConc, y = Sig)
        )
        p <- myCalPlot +
            geom_point() +
            geom_smooth(method = "lm", formula = "y ~ x", mapping = aes(weight = wts)) +
            labs(title = myTitle, subtitle = "weighted fit") +
            ylab("Cq") + xlab("Concentration in IU/mL") +
            scale_y_reverse() +
            scale_x_log10(
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))
            ) +
            theme_bw() +
            theme(
                panel.grid.minor = element_blank()
            ) +
            labs(color = "Sample") +
            geom_rug(data = myData[which(!myData$isCalib), ], aes(x = 10^logConcPredicted, color = SamName))

        suppressWarnings(print(p))
    }
}

# this function will plot data for a single sample passed in
sampleByRank <- function(aSample) {
    # make a composite Lab:Target variable, this will be X
    aSample$LabTarg <- factor(paste(aSample$Lab, ":", aSample$Target, sep = ""))

    # order this factor by the median values of 10^logConcPredicted
    aSample$LabTarg <- fct_reorder(aSample$LabTarg, 10^aSample$logConcPredicted, .fun = median, na.rm = TRUE)

    # now calculate a naive robust estimate of the median of the medians
    aSampleGroup <- group_by(aSample, LabTarg)

    # get the medians for each LabTarg
    aSampleMeds <- aSampleGroup %>%
        mutate(meds = median(logConcPredicted, na.rm = TRUE)) %>%
        summarise(smed = median(meds))

    # get the median of the medians
    aSampleMedian <- median(aSampleMeds$smed, na.rm = TRUE)

    # get the mad of the medians
    aSampleMad <- mad(aSampleMeds$smed, na.rm = TRUE)

    # create a ggplot object of the sample values
    myValues <- ggplot(aSample, aes(x = LabTarg, y = 10^logConcPredicted))

    # make a boxplot
    p <- myValues +
        geom_boxplot() +
        geom_jitter(aes(color = LabTarg, alpha = 0.4, shape = Batch)) +
        geom_hline(yintercept = 10^aSampleMedian, colour = "orange", linetype = 2, size = 1.2) +
        geom_hline(yintercept = 10^(aSampleMedian - aSampleMad), colour = "grey70", linetype = 2, size = 0.8) +
        geom_hline(yintercept = 10^(aSampleMedian + aSampleMad), colour = "grey70", linetype = 2, size = 0.8) +
        theme_bw() +
        scale_y_log10(
            limits = c(10^(aSampleMedian - 1), 10^(aSampleMedian + 1)),
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x)),
        ) +
        annotation_logticks(sides = "l") +
        xlab("Lab:Target") + ylab("Concentration in IU/mL") +
        labs(title = paste(aSample$SamName), subtitle = "Predictions from calibration to WHO IS -- Orange line is median of all results") +
        theme(
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(colour = "black", size = rel(1.05), face = "plain"),
            axis.title = (element_text(colour = "black", size = rel(1.5), face = "bold")),
            strip.text = (element_text(size = rel(1.4))),
            legend.position = "none"
        ) +
        guides(col = guide_legend(ncol = 1, title = NULL), alpha = FALSE)

    suppressWarnings(print(p))
}

# this function will plot data for a single sample passed in,
# by lab, to try to see lab effects
sampleByLabTarg <- function(aSample) {
    # make a composite Lab:Target variable, this will be X
    aSample$LabTarg <- factor(paste(aSample$Lab, ":", aSample$Target, sep = ""))

    # now calculate a naive robust estimate of the median of the medians
    aSampleGroup <- group_by(aSample, LabTarg)

    # get the medians for each LabTarg
    aSampleMeds <- aSampleGroup %>%
        mutate(meds = median(logConcPredicted, na.rm = TRUE)) %>%
        summarise(smed = median(meds))

    # get the median of the medians
    aSampleMedian <- median(aSampleMeds$smed, na.rm = TRUE)

    # get the mad of the medians
    aSampleMad <- mad(aSampleMeds$smed, na.rm = TRUE)

    # create a ggplot object of the sample values
    myValues <- ggplot(aSample, aes(x = LabTarg, y = 10^logConcPredicted))

    # make a boxplot
    p <- myValues +
        geom_boxplot() +
        geom_jitter(aes(color = LabTarg, alpha = 0.4, shape = Batch)) +
        # geom_rug(sides="r", aes(color = LabTarg)) +
        geom_hline(yintercept = 10^aSampleMedian, colour = "orange", linetype = 2, size = 1.2) +
        geom_hline(yintercept = 10^(aSampleMedian - aSampleMad), colour = "grey70", linetype = 2, size = 0.8) +
        geom_hline(yintercept = 10^(aSampleMedian + aSampleMad), colour = "grey70", linetype = 2, size = 0.8) +
        theme_bw() +
        scale_y_log10(
            limits = c(10^(aSampleMedian - 1), 10^(aSampleMedian + 1)),
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x)),
        ) +
        annotation_logticks(sides = "l") +
        xlab("Lab:Target") + ylab("Concentration in IU/mL") +
        labs(title = paste(aSample$SamName), subtitle = "Predictions from calibration to WHO IS -- Orange line is median of all results") +
        theme(
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text = element_text(colour = "black", size = rel(1.05), face = "plain"),
            axis.title = (element_text(colour = "black", size = rel(1.5), face = "bold")),
            strip.text = (element_text(size = rel(1.4))),
            legend.position = "none"
        ) +
        guides(col = guide_legend(ncol = 1, title = NULL), alpha = FALSE)

    suppressWarnings(print(p))
}

# build and plot the deviations for each lab:target from global medians
plotDeviations <- function(myResults) {
    # first get a local table of lab:target medians from materialSummary fn
    ltMeds <- materialTable(myResults)

    # get a column of lab:target
    ltMeds$LabTarg <- factor(paste(ltMeds$Lab, ":", ltMeds$Target, sep = ""))

    ltDevs <- ltMeds # use this to hold the deviations

    # now get a local table of material medians
    matMeds <- materialSummary(myResults)


    # clumsy, but loop over the lab:targets and materials
    for (matl in 1:dim(matMeds)[1]) { # loop over materials
        for (lt in 1:dim(ltMeds)[1]) { # loop over lab:targ
            ltDevs[lt, 2 + matl] <- ltMeds[lt, 2 + matl] - matMeds[matl, 2]
        }
    }

    # pivot to get deviations
    allDevs <- pivot_longer(ltDevs, cols = 3:10, names_to = "Material", values_to = "Deviation")

    pDevs <- ggplot(allDevs, aes(y = Deviation, x = Material))

    p <- pDevs +
        geom_point() +
        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            # axis.text = element_text(colour = "black", size = rel(1.05), face = "plain"),
            # axis.title = (element_text(colour = "black", size = rel(1.5), face = "bold")),
            # strip.text = (element_text(size=rel(1.4))),
            legend.position = "none"
        ) +
        ylim(c(-1, 1)) +
        ylab("Deviation from Median, log10 IU/mL") +
        facet_wrap(vars(LabTarg), nrow = 5)

    print(p)

    return(ltDevs)
}



# function to plot material median results against the nominal values
comparePlot <- function(samMedian) {
    # here are the nominal values in alphabetical order, in log10 copies/mL
    myNominals <- c(10.30103, 4, 5.195899652, 3.698970004, 4.505149978, 6.73, 3.698970004, 4.698970004)

    # here's a nice plottable data frame
    myCompar <- cbind(samMedian, myNominals)
    names(myCompar) <- c("Material", "studyVal", "studyValCI", "nominal")

    # make the ggplot object
    myCompPlot <- ggplot(myCompar, aes(x = nominal, y = studyVal, color = Material))

    # and a pointrange (error bar) plot
    myCompPlot +
        geom_pointrange(aes(ymin = studyVal - studyValCI, ymax = studyVal + studyValCI), alpha = 0.7) +
        theme_bw() +
        ylab("Study Results (log10 IU/mL)") +
        xlab("Nominal Value (log10 genome copies/mL)") +
        theme(legend.position = "bottom") +
        geom_abline(intercept = 0, slope = 7.7 / 8, color = "grey")
}
