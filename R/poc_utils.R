# poc_utils.R contains all functions required for computing and plotting the Proof of Concept data.

plot_poc_results <- function(shared_axis) {
  
  # read and calibrate the PoC results
  PoC_results = computePoC()
  
  # build the plot -- we can plot it with free_y or fixed scales
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


computePoC <- function( ) {
  
  # first read the data
  PoC <- getPoC()

  # remove the Asuragen data, it wasn't calibrated
  PoC <- subset(PoC, SamName != "Asuragen")
  
  # make a column of consistent log-scale signals
  # this is Cq and log(dPCR copies)
  PoC$logVal = PoC$Raw.Value
  PoC$logVal[PoC$dPCR == TRUE] = log(PoC$logVal[PoC$dPCR == TRUE], 10)
  
  # make a column of LabTargExpt
  PoC$LabTargExpt = paste(PoC$Lab, PoC$Target, PoC$Experiment, sep=":")
  
  # make a dataframe for only the calibration data
  PoC_rawcal = subset(PoC, PoC$isCalib)
  
  # make a dataframe for only the sample data
  PoC_rawsample = subset(PoC, !PoC$isCtrl)
  
  # group calibration data by Lab:Target
  calGroup <- PoC_rawcal %>% group_by(LabTargExpt, SamName)
  
  # get the medians for each LabTarg
  PoCCalVal <- calGroup %>% mutate( avg = mean( logVal, na.rm=T)) %>% summarize( mmean = mean(avg))
  
  # Here are the values in IU for each calibrant 
    # Asuragen dropped -- the sample shared was NOT CALIBRATED, we 
    # don't have a clue about its IU concentration  
  IUvals = data.frame(
    SamName = c( "AcroMetrix", "LGC"),
    log10IU = c(4.033 , 3.911),
    CI = c(0.46, 0.315)
  )
  
  # get the calibrant values in the grouped data
  PoCCalVal <- PoCCalVal %>% left_join( IUvals )
  
  # calculate the calibrated results
  
  # combine calibration data with raw data
  PoC_to_calib <- left_join( PoC_rawsample, PoCCalVal, by = "LabTargExpt" )
  
  # change calibrant column name
  PoC_to_calib <- rename(PoC_to_calib, Calibrant = SamName.y )
  
  # calculate the qPCR calibrated values
  PoC_calibrated_qPCR <- PoC_to_calib %>% filter( !dPCR ) %>% mutate( log10IU_Sample = log( 10^log10IU * 2^( mmean - logVal ) , 10))
  
  # calculate the dPCR calibrated values
  PoC_calibrated_dPCR <- PoC_to_calib %>% filter( dPCR ) %>% 
    mutate( log10IU_Sample = log(10^logVal / (10^mmean/10^log10IU), 10))
  
  # merge qPCR and d PCR results
  PoC_results <- bind_rows( PoC_calibrated_qPCR, PoC_calibrated_dPCR)
  
  # get rid of no-signal rows -- this will clean up the faceting of graphs by sample
  PoC_results = subset( PoC_results, !is.na(PoC_results$logVal))
  
  # make a column of LabTargExptCalib
  PoC_results$LabTargExptCalib = paste(PoC_results$LabTargExpt, PoC_results$Calibrant, sep=":")
  
  return( PoC_results )
}  


