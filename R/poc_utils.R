# poc_utils.R contains all functions required for plotting and calibrating the Proof of Concept data.

plot_poc_results <- function(shared_axis) {
  PoC <- getPoC()
  # remove the Asuragen data -- the sample shared was NOT CALIBRATED, we
  # don't have a clue about its IU concentration
  PoC <- subset(PoC, SamName != "Asuragen")

  # make a column of consistent log-scale signals
  # this is Cq and log(dPCR copies)
  PoC$logVal <- PoC$Raw.Value
  PoC$logVal[PoC$dPCR == TRUE] <- log(PoC$logVal[PoC$dPCR == TRUE], 10)

  # make a column of LabTargExpt
  PoC$LabTargExpt <- paste(PoC$Lab, PoC$Target, PoC$Experiment, sep = ":")

  # make a dataframe for only the calibration data
  PoC_rawcal <- subset(PoC, PoC$isCalib)

  # make a dataframe for only the sample data
  PoC_rawsample <- subset(PoC, !PoC$isCtrl)

  # group calibration data by Lab:Target
  calGroup <- PoC_rawcal %>% group_by(LabTargExpt, SamName)

  # get the medians for each LabTarg
  PoCCalVal <- calGroup %>%
    mutate(avg = mean(logVal, na.rm = T)) %>%
    summarize(mmean = mean(avg))

  # Here are the values in IU for each calibrant
  IUvals <- data.frame(
    SamName = c("AcroMetrix", "LGC"),
    log10IU = c(4.033, 3.911),
    CI = c(0.46, 0.315)
  )

  # get the calibrant values in the grouped data
  PoCCalVal <- PoCCalVal %>% left_join(IUvals)

  # calculate the calibrated results

  # combine calibration data with raw data
  PoC_to_calib <- left_join(PoC_rawsample, PoCCalVal, by = "LabTargExpt")

  # change calibrant column name
  PoC_to_calib <- rename(PoC_to_calib, Calibrant = SamName.y)

  # calculate the qPCR calibrated values
  PoC_calibrated_qPCR <- PoC_to_calib %>%
    filter(!dPCR) %>%
    mutate(log10IU_Sample = log(10^log10IU * 2^(mmean - logVal), 10))

  # calculate the dPCR calibrated values
  PoC_calibrated_dPCR <- PoC_to_calib %>%
    filter(dPCR) %>%
    mutate(log10IU_Sample = log(10^logVal / (10^mmean / 10^log10IU), 10))

  # merge qPCR and d PCR results
  PoC_results <- bind_rows(PoC_calibrated_qPCR, PoC_calibrated_dPCR)

  # get rid of no-signal rows -- this will clean up the faceting of graphs by sample
  PoC_results <- subset(PoC_results, !is.na(PoC_results$logVal))

  # make a column of LabTargExptCalib
  PoC_results$LabTargExptCalib <- paste(PoC_results$LabTargExpt, PoC_results$Calibrant, sep = ":")

  myPoC_bySample <- ggplot(PoC_results, aes(x = LabTargExpt, y = log10IU_Sample))

  if (shared_axis) {
    axis_scale_setting <- "fixed"
  } else {
    axis_scale_setting <- "free_y"
  }

  myPoC_bySample + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(legend.position = "bottom") +
    geom_jitter(aes(color = Calibrant)) +
    facet_wrap(~SamName.x, scales = axis_scale_setting) +
    scale_colour_manual(values = c("#3182bd", "#fdae6b")) +
    ggtitle("Calibrated results on Clinical Samples") +
    ylab("log10 IU") + xlab("Lab:Target:Run:Calibrant")
}


calibPoC <- function() {
  # first read the data
  PoC <- read_csv("..\\data\\20230223 PoC.csv",
    col_types = cols(
      Replicate = col_integer(),
      Experiment = col_integer()
    )
  )

  # remove the Asuragen data
  PoC <- subset(PoC, SamName != "Asuragen")

  # make a column of consistent log-scale signals
  # this is Cq and log(dPCR copies)
  PoC$logVal <- PoC$Raw.Value
  PoC$logVal[PoC$dPCR == TRUE] <- log(PoC$logVal[PoC$dPCR == TRUE], 10)

  # make a column of LabTargExpt
  PoC$LabTargExpt <- paste(PoC$Lab, PoC$Target, PoC$Experiment, sep = ":")

  # make a dataframe for only the calibration data
  PoC_rawcal <- subset(PoC, PoC$isCalib)

  # make a dataframe for only the sample data
  PoC_rawsample <- subset(PoC, !PoC$isCtrl)

  # group calibration data by Lab:Target
  calGroup <- PoC_rawcal %>% group_by(LabTargExpt, SamName)

  # get the medians for each LabTarg
  PoCCalVal <- calGroup %>%
    mutate(avg = mean(logVal, na.rm = T)) %>%
    summarize(mmean = mean(avg))

  # Here are the values in IU for each calibrant
  IUvals <- data.frame(
    SamName = c("AcroMetrix", "LGC"),
    log10IU = c(4.033, 3.911),
    CI = c(0.46, 0.315)
  )

  # get the calibrant values in the grouped data
  PoCCalVal <- PoCCalVal %>% left_join(IUvals)

  # calculate the calibrated results

  # combine calibration data with raw data
  PoC_to_calib <- left_join(PoC_rawsample, PoCCalVal, by = "LabTargExpt")

  # change calibrant column name
  PoC_to_calib <- rename(PoC_to_calib, Calibrant = SamName.y)

  # calculate the qPCR calibrated values
  PoC_calibrated_qPCR <- PoC_to_calib %>%
    filter(!dPCR) %>%
    mutate(log10IU_Sample = log(10^log10IU * 2^(mmean - logVal), 10))

  # calculate the dPCR calibrated values
  PoC_calibrated_dPCR <- PoC_to_calib %>%
    filter(dPCR) %>%
    mutate(log10IU_Sample = log(10^logVal / (10^mmean / 10^log10IU), 10))

  # merge qPCR and d PCR results
  PoC_results <- bind_rows(PoC_calibrated_qPCR, PoC_calibrated_dPCR)

  # Mon Aug 29 18:47:54 2022 ------------------------------
  # get rid of no-signal rows -- this will clean up the faceting of graphs by sample
  PoC_results <- subset(PoC_results, !is.na(PoC_results$logVal))

  # Thu Oct  6 16:59:48 2022 ------------------------------
  # make a column of LabTargExptCalib
  PoC_results$LabTargExptCalib <- paste(PoC_results$LabTargExpt, PoC_results$Calibrant, sep = ":")

  # Mon Aug 29 18:48:04 2022 ------------------------------

  # Mon Aug 29 20:18:36 2022 ------------------------------
  # Thu Sep  1 12:34:32 2022 ------------------------------
  # DIDN'T WORK, normalised signal didn't make sense --
  # make a tibble with a column of normalized log signal -- can we plot dPCR and Cq on the same graph?
  #  PoC_normSig <- PoC_results %>% group_by(LabTargExpt)
  #  PoC_normSig <- PoC_normSig %>% mutate( normSig = (log10( (10^logVal - 10^min(logVal)) / 10^max(logVal))))
  # Mon Aug 29 20:18:36 2022 ------------------------------


  # make graphs!
  # Figure 3

  # Mon Aug 29 16:51:32 2022 ------------------------------
  # Thu Sep  1 12:07:26 2022 ------------------------------
  # Thu Sep  1 12:06:52 2022 ------------------------------
  # Left side of possible Fig 3 -- uncalibrated Cq only
  myPoC_Cq_bySample <- ggplot(subset(PoC_results, (Calibrant == "LGC" & !dPCR)), aes(x = LabTargExpt, y = logVal))

  myPoC_Cq_bySample + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter((aes(colour = Lab))) +
    scale_colour_manual(values = c("#31a354", "#c51b8a")) +
    facet_wrap(~SamName.x) +
    ggtitle("Cq on Clinical Samples") +
    ylab("Cq") + xlab("Lab:Target:Run") +
    scale_y_reverse()

  # for ddPCR
  myPoC_copy_bySample <- ggplot(subset(PoC_results, (Calibrant == "LGC" & dPCR == TRUE)), aes(x = LabTargExpt, y = logVal))

  myPoC_copy_bySample + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter(aes(color = as.factor(Experiment))) +
    # scale_colour_manual(values=c("#31a354", "#c51b8a")) +
    facet_wrap(~SamName.x) +
    ggtitle("ddPCR on Clinical Samples") +
    ylab("log Copy") + xlab("Lab:Target:Run")
  # scale_y_reverse()

  # Right side of possible Fig 3

  myPoC_bySample <- ggplot(PoC_results, aes(x = LabTargExpt, y = log10IU_Sample))

  myPoC_bySample + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter(aes(color = Calibrant)) +
    scale_colour_manual(values = c("#3182bd", "#fdae6b")) +
    facet_wrap(~SamName.x) +
    ggtitle("Calibrated results on Clinical Samples") +
    ylab("log10 IU") + xlab("Lab:Target:Run")

  ggsave("allCalibCommonScale.pdf", width = 8, height = 10, units = "in", device = "pdf")

  # Fri Feb 24 10:33:47 2023 ------------------------------
  # plot a graph by sample, faceted by calibrant to look for bias

  myPoC_bySampleAll <- ggplot(PoC_results, aes(x = SamName.x, y = log10IU_Sample))

  # only points, not box plot
  myPoC_bySampleAll + theme_bw() +
    # geom_boxplot( aes(colour = LabTarg) ) +
    geom_jitter(aes(colour = LabTargExpt, alpha = 0.2), width = 0.25) +
    facet_grid(Calibrant ~ Lab) +
    ggtitle("All Results, Calibrated to IU")

  ggsave("all Labs v Calibrant.pdf", width = 8, height = 10, units = "in", device = "pdf")


  # Thu Oct  6 14:40:59 2022 ------------------------------
  # new analysis -- before and after cal for LANL
  # LANL only before and after

  myLANLb4 <- ggplot(subset(PoC_results, (Lab == "LANL" & Calibrant == "LGC")), aes(x = LabTargExpt, y = logVal))

  myLANLb4 + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter((aes(colour = Lab))) +
    scale_colour_manual(values = c("#31a354", "#c51b8a")) +
    facet_wrap(~SamName.x) +
    ggtitle("Cq on Clinical Samples") +
    ylab("Cq") + xlab("Lab:Target:Run") +
    scale_y_reverse()

  ggsave("LANLbefore.pdf", width = 8, height = 10, units = "in", device = "pdf")

  myLANLafter <- ggplot(subset(PoC_results, (Lab == "LANL")), aes(x = LabTargExptCalib, y = log10IU_Sample))

  myLANLafter + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter((aes(colour = Calibrant))) +
    scale_colour_manual(values = c("#31a354", "#c51b8a")) +
    facet_wrap(~SamName.x) +
    ggtitle("Calibrated Results on Clinical Samples") +
    ylab("log10 IU") + xlab("Lab:Target:Run:Calibrant")

  ggsave("LANLafter.pdf", width = 8, height = 10, units = "in", device = "pdf")

  # Thu Feb 23 13:00:57 2023 ------------------------------
  # new analysis -- before and after cal for Biogazelle
  # Biogazelle only before and after

  myGazb4 <- ggplot(subset(PoC_results, (Lab == "Biogazelle" & Calibrant == "LGC")), aes(x = LabTargExpt, y = logVal))

  myGazb4 + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter((aes(colour = Lab))) +
    scale_colour_manual(values = c("#31a354", "#c51b8a")) +
    facet_wrap(~SamName.x) +
    ggtitle("Cq on Clinical Samples") +
    ylab("Cq") + xlab("Lab:Target:Run") +
    scale_y_reverse()

  ggsave("GazBefore.pdf", width = 8, height = 10, units = "in", device = "pdf")

  myGazAfter <- ggplot(subset(PoC_results, (Lab == "Biogazelle")), aes(x = LabTargExptCalib, y = log10IU_Sample))

  myGazAfter + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter((aes(colour = Calibrant))) +
    scale_colour_manual(values = c("#31a354", "#c51b8a")) +
    facet_wrap(~SamName.x) +
    ggtitle("Calibrated Results on Clinical Samples") +
    ylab("log10 IU") + xlab("Lab:Target:Run:Calibrant")

  ggsave("GazAfter.pdf", width = 8, height = 10, units = "in", device = "pdf")

  # Thu Oct  6 16:50:37 2022 ------------------------------
  # Biodesix b4 and after

  myBioDb4 <- ggplot(subset(PoC_results, (Lab == "Biodesix" & Calibrant == "LGC")), aes(x = LabTargExpt, y = logVal))

  myBioDb4 + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter((aes(colour = Lab))) +
    scale_colour_manual(values = c("#31a354", "#c51b8a")) +
    facet_wrap(~SamName.x) +
    ggtitle("ddPCR Counts on Clinical Samples") +
    ylab("ddPCR Copies") + xlab("Lab:Target:Run")

  ggsave("BiodBefore.pdf", width = 8, height = 10, units = "in", device = "pdf")

  myBioDafter <- ggplot(subset(PoC_results, (Lab == "Biodesix")), aes(x = LabTargExptCalib, y = log10IU_Sample))

  myBioDafter + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter((aes(colour = Calibrant))) +
    scale_colour_manual(values = c("#31a354", "#c51b8a")) +
    facet_wrap(~SamName.x) +
    ggtitle("Calibrated Results on Clinical Samples") +
    ylab("log10 IU") + xlab("Lab:Target:Run:Calibrant")

  ggsave("BioDafter.pdf", width = 8, height = 10, units = "in", device = "pdf")

  # Fri Oct  7 10:33:46 2022 ------------------------------

  # let's do an after plot segregated by calibrant

  myAllafter <- ggplot(PoC_results, aes(x = LabTargExptCalib, y = log10IU_Sample))

  myAllafter + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(legend.position = "bottom") +
    geom_jitter((aes(colour = Calibrant, shape = Lab))) +
    scale_colour_manual(values = c("#31a354", "#c51b8a")) +
    facet_wrap(~SamName.x, scales = "free_y") +
    ggtitle("Calibrated Results on Clinical Samples") +
    ylab("log10 IU") + xlab("Lab:Target:Run:Calibrant")

  ggsave("AllafterFreeScale.pdf", width = 8, height = 10, units = "in", device = "pdf")

  return(myPoC_bySample)


  # Tue Dec 13 16:45:14 2022 ------------------------------
  # redo allafter with a range of 10^4
  # this graph is an unforgiving view of the biases
  # Fri Feb 24 10:31:00 2023 ------------------------------
  # this seems broken today -- doesn't find "facet_grid_sc"

  myAllafterSameRange <-
    myAllafter + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(legend.position = "bottom") +
    geom_jitter((aes(colour = Calibrant, shape = Lab))) +
    scale_colour_manual(values = c("#31a354", "#c51b8a")) +
    facet_wrap(~SamName.x, scales = "free_y") +
    ggtitle("Calibrated Results on Clinical Samples") +
    ylab("log10 IU") + xlab("Lab:Target:Run:Calibrant")

  scales_our_y <- list(
    "s01" = scale_y_continuous(limits = c(7, 11)),
    "s02" = scale_y_continuous(limits = c(7, 11)),
    "s03" = scale_y_continuous(limits = c(7, 11)),
    "s04" = scale_y_continuous(limits = c(7, 11)),
    "s05" = scale_y_continuous(limits = c(7, 11)),
    "s06" = scale_y_continuous(limits = c(7, 11)),
    "s07" = scale_y_continuous(limits = c(5, 9)),
    "s08" = scale_y_continuous(limits = c(4, 8)),
    "s09" = scale_y_continuous(limits = c(4, 8)),
    "s10" = scale_y_continuous(limits = c(3, 7)),
    "s11" = scale_y_continuous(limits = c(0, 4)),
    "s12" = scale_y_continuous(limits = c(0, 4))
  )
  myAllafter + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(legend.position = "bottom") +
    geom_jitter((aes(colour = Calibrant, shape = Lab))) +
    scale_colour_manual(values = c("#31a354", "#c51b8a")) +
    facet_grid_sc(rows = SamName.x, scales = list(y = scales_our_y)) +
    ggtitle("Calibrated Results on Clinical Samples") +
    ylab("log10 IU") + xlab("Lab:Target:Run:Calibrant")



  # ***************
  # these are mostly obsolete graphs and analyses
  # ***************

  # Fri Oct  7 11:43:33 2022 ------------------------------
  # try anova

  # LANL
  # design is different levels of target and different calibrants
  # levels of Lab, Target, Calibrant, SamName.x -- dependent var == log10IU

  # create a new dataframe
  myLANLresults <- data.frame()

  # Possible Supplement Graphs

  # make smoothed density plots
  PoC_distrib <- ggplot(PoC_results, aes(x = log10IU_Sample))

  PoC_distrib + geom_density(na.rm = TRUE, aes(fill = Calibrant, alpha = 0.2), color = NA) + facet_grid(LabTargExpt ~ .) +
    theme_bw() + xlab("Log10 IU") + geom_rug(aes(colour = Calibrant)) +
    ggtitle("PoC Results: Clinical Samples") +
    theme(aspect.ratio = 0.5)

  # make a correlation matrix
  # need wide results
  PoC_resultsW <- pivot_wider(PoC_results[c(1:6, 9, 7, 11:12, 16)], names_from = "Calibrant", values_from = "log10IU_Sample")

  # make the correlation matrix WITHOUT border columns
  myCorr <- ggpairs(PoC_resultsW,
    columns = c(10:11),
    aes(colour = LabTargExpt, alpha = 0.1),
    axisLabels = "show"
  )

  # plot it
  myCorr + theme_bw() + theme(aspect.ratio = 1) +
    ggtitle("PoC Results: Clinical Samples")

  # make the correlation matrix with border columns
  myCorr <- ggpairs(PoC_resultsW,
    columns = c(9, 10:11), aes(colour = LabTargExpt, alpha = 0.2),
    lower = list(continuous = "density", combo = "dot_no_facet"),
    upper = list(continuous = "points", combo = "box_no_facet"),
    axisLabels = "internal"
  )

  # plot it
  myCorr + theme_bw() + theme(aspect.ratio = 1) +
    ggtitle("PoC Results: Clinical Samples") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  # Wed Aug 17 16:36:58 2022 ------------------------------
  # Simplify our story -- box plot per sample, per lab:targ

  # data come from PoC_results
  # myPoC_boxes = ggplot(PoC_results, aes( x = SamName.x, y = log10IU_Sample))
  #
  # myPoC_boxes + theme_bw() +
  #   geom_boxplot(aes(colour = LabTargExpt)) +
  #   facet_wrap(.~Calibrant)
  #
  # myPoC_boxes + theme_bw() +
  #   geom_boxplot(aes(colour = Calibrant)) +
  #   facet_grid(LabTargExpt ~ . )
  #
  # myPoC_boxes + theme_bw() +
  #   geom_boxplot(aes(colour = Calibrant)) +
  #   facet_grid(.~LabTargExpt ) +
  #   geom_jitter(aes(colour=Calibrant))


  myPoC_allUncal_bySample <-
    ggplot(subset(PoC_normSig, Calibrant == "LGC"), aes(x = LabTargExpt, y = normSig))

  myPoC_allUncal_bySample + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter((aes(colour = Lab))) +
    facet_wrap(~SamName.x) +
    ggtitle("Normalized Signal on Clinical Samples") +
    ylab("log Signal") + xlab("Lab:Target:Run") +
    scale_y_reverse()

  # let's look only at calibrated qPCR v. Cq
  myPoC_cq <- ggplot(subset(PoC_results, (Calibrant == "LGC" & !dPCR)), aes(x = SamName.x, y = logVal))

  myPoC_cq + theme_bw() +
    geom_boxplot(aes(colour = LabTargExpt)) +
    geom_jitter(aes(alpha = 0.2))

  # Mon Aug 29 15:45:30 2022 ------------------------------
  # jitterred points only...
  # let's look only at calibrated qPCR v. Cq
  myPoC_cq + theme_bw() +
    # geom_boxplot( aes(colour = LabTargExpt) ) +
    geom_jitter(aes(colour = LabTargExpt, alpha = 0.2))

  myPoC_PCR_calib <- ggplot(subset(PoC_results, (!dPCR)), aes(x = SamName.x, y = log10IU_Sample))

  myPoC_cq + theme_bw() +
    # geom_boxplot( aes(colour = LabTarg) ) +
    geom_jitter(aes(colour = LabTargExpt, alpha = 0.2)) +
    facet_grid(Calibrant ~ .) +
    ggtitle("qPCR only Results, uncalibrated, Cq") +
    ylab("Cq") + xlab("Sample")

  myPoC_qPCR <- ggplot(subset(PoC_results, (Calibrant == "LGC" & !dPCR)), aes(y = logVal))

  myPoC_qPCR + theme_bw() +
    geom_histogram(aes(fill = LabTargExpt, alpha = 0.1), position = "dodge") +
    ggtitle("qPCR only Results, uncalibrated, Cq") +
    xlab("Cq")

  myPoC_qPCR + theme_bw() + scale_x_reverse() +
    geom_density(aes(fill = LabTargExpt, alpha = 0.1)) +
    ggtitle("qPCR only Results, uncalibrated, Cq") +
    xlab("Cq") +
    facet_grid(LabTargExpt ~ .)

  # only for qPCR, density plot
  myPoC_qPCR_Calib <- ggplot(subset(PoC_results, !dPCR), aes(x = log10IU_Sample))

  myPoC_qPCR_Calib + theme_bw() +
    geom_density(aes(fill = Calibrant), alpha = 0.3, color = NA) +
    ggtitle("qPCR only Results, CALIBRATED, log10 IU") +
    xlab("log10 IU") +
    facet_grid(LabTargExpt ~ .)

  # include dPCR, density plot
  myPoC_all_Calib <- ggplot(PoC_results, aes(x = log10IU_Sample))

  myPoC_all_Calib + theme_bw() +
    geom_density(aes(fill = Calibrant), alpha = 0.3, color = NA) +
    ggtitle("All Results, CALIBRATED, log10 IU") +
    xlab("log10 IU") +
    facet_grid(LabTargExpt ~ .)

  # only points, not box plot
  myPoC_PCR_calib + theme_bw() +
    # geom_boxplot( aes(colour = LabTarg) ) +
    geom_jitter(aes(colour = LabTargExpt, alpha = 0.2), width = 0.25) +
    facet_grid(Calibrant ~ .) +
    ggtitle("qPCR only Results, calibrated to IU")

  # only points, not box plot
  myPoC_PCR_calib + theme_bw() +
    # geom_boxplot( aes(colour = LabTarg) ) +
    geom_jitter(aes(colour = LabTargExpt, alpha = 0.2), width = 0.25) +
    facet_grid(Calibrant ~ .) +
    ggtitle("qPCR only Results, calibrated to IU")

  # show dPCR and qPCR in facets with free_y scaling
  PoC_all_graph <- ggplot(PoC_results, aes(x = SamName.x, y = log10IU_Sample))

  PoC_all_graph + theme_bw() +
    geom_jitter(aes(alpha = 0.15, shape = Calibrant, colour = as.factor(Experiment)), size = 2) +
    facet_grid(~LabTargExpt)
}

# Thu Aug 18 14:23:25 2022 ------------------------------
# Publication-ready graphs

# Figure 1 needs some rank-ordered histograms
# maybe all 8 faceted in a single figure?

# Figure 2 needs the "results" figure of IU v. Nominal Copy

# Figure 3 needs demonstration of the harmonization hypothesis
# Density figure, each assay (lab:target:experiment)
# Plot of harmonized data points -- faceted by assay

pubPlots <- function() {
  # get the data in scope
  myResults <- getStudyData()

  # Now compute the results
  i <- 0
  myResults <- study
  for (Lab in levels(study$Lab)) {
    newLab <- calAndPredict(study[which(study$Lab == Lab), ])
    if (i == 0) {
      myResults <- newLab
    } else {
      myResults <- rbind(myResults, newLab)
    }
    i <- i + 1
  }

  # get a factor list of sample names (not WHO_IS)
  myMaterials <- levels(fct_drop(factor(myResults$SamName[which(!myResults$isCalib & !myResults$isCtrl)])))

  # iterate the samples to get rank-ordered plots for Fig 1 panel B
  for (Matl in myMaterials) {
    # make the rank-ordered box plot for each material
    sampleByRank(myResults[which(myResults$SamName == Matl), ])

    # save a PDF of the graph, with the material as the base filename
    ggsave(paste(Matl, ".pdf", sep = ""), width = 9.60, height = 5.94, units = "in", device = "pdf")
  }
}
