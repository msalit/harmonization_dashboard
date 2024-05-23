# RNA Harmonization Study Analysis -- Publication-ready Graphs
# Marc Salit
# 22 May 2024

# Publication-ready graphs

# Figure 1 needs some rank-ordered histograms
# Made these as 8 separate .tiff files

# Figure 2a needs the "results" figure of IU v. Nominal Copy as a nice display item

# Figure 2b is a demonstration of the harmonization hypothesis with PoC results

# make and save a ranked-by-Lab:Target box plot for each material
pubBoxPlots1b <- function() {
  
  # register system fonts for the publication-ready box plots
  #library(extrafont)
  #library(extrafontdb)
  #font_import()
  
  # get the data in scope
  #myResults = getStudyData()
  study = getStudyData()
  
  # Now compute the results
  i = 0
  myResults = study
  for( Lab in levels(study$Lab)) {
    
    newLab = calAndPredict(study[which(study$Lab == Lab),], FALSE)
    if (i == 0)
      myResults = newLab
    else
      myResults = rbind(myResults, newLab)
    i = i+1
  }
  
  # get a factor list of sample names (not WHO_IS)
  myMaterials = levels(fct_drop(factor(myResults$SamName[which(!myResults$isCalib & !myResults$isCtrl)])))
  
  # iterate the samples to get rank-ordered plots for Fig 1 panel B
  for( Matl in myMaterials ) {
    
    # make the rank-ordered box plot for each material
    sampleByRankForPub(myResults[which(myResults$SamName == Matl),])
    
    # save a PDF of the graph, with the material as the base filename 
    ggsave(paste(Matl, ".tiff", sep=""), width=180, height = 180, units ="mm", dpi=600, device="tiff" ) 
    #embed_fonts(paste(Matl, ".pdf", sep=""))
    
  }
}


# Tue Apr 30 14:39:04 2024 ------------------------------
# this function will plot data for a single sample passed in
sampleByRankForPub <- function( aSample ) {
  
  # make a composite Lab:Target variable, this will be X
  aSample$LabTarg = factor(paste(aSample$Lab, ':', aSample$Target, sep=''))
  
  # order this factor by the median values of 10^logConcPredicted
  aSample$LabTarg = fct_reorder( aSample$LabTarg, 10^aSample$logConcPredicted, .fun = median, na.rm=TRUE )
  
  # now calculate a naive robust estimate of the median of the medians
  aSampleGroup <- group_by(aSample, LabTarg)
  
  # get the medians for each LabTarg
  aSampleMeds = aSampleGroup %>% mutate(meds = median(logConcPredicted, na.rm = TRUE)) %>% summarise(smed = median(meds))
  
  # get the median of the medians
  aSampleMedian = median(aSampleMeds$smed, na.rm = TRUE)
  
  # get the mad of the medians
  aSampleMad = mad(aSampleMeds$smed, na.rm = TRUE)
  
  # make a plot title dummy variable
  aSample$tempTitle = paste(aSample$SamName)
  
  # create a ggplot object of the sample values
  myValues = ggplot( aSample, aes(x=LabTarg, y=10^logConcPredicted))
  
  # make a boxplot
  p = myValues + 
    geom_boxplot() + 
    geom_jitter(aes(color = LabTarg, alpha =0.4, shape = Batch)) +
    #geom_rug(sides="r", aes(color = LabTarg)) +
    geom_hline(yintercept = 10^aSampleMedian, colour = "orange", linetype = 2, size = 1.2) +
    #geom_ribbon( aes(x=LabTarg, ymin = 10^(aSampleMedian - aSampleMad),  ymax = 10^(aSampleMedian + aSampleMad)), fill = "grey30" ) +
    #geom_ribbon( aes(x=LabTarg, ymin = 10^(aSampleMedian - 0.8),  ymax = 10^(aSampleMedian + 0.8)), fill = "grey70" ) +
    geom_hline(yintercept = 10^(aSampleMedian - aSampleMad), colour = "grey70", linetype = 2, size = 0.8) +
    geom_hline(yintercept = 10^(aSampleMedian + aSampleMad), colour = "grey70", linetype = 2, size = 0.8) +
    theme_bw() + 
    scale_y_log10(
      limits = c( 10^(aSampleMedian -1 ), 10^(aSampleMedian + 1) ),
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
    )  +   
    annotation_logticks(sides = "l") +
    #xlab("Lab:Target") + ylab("Concentration in IU/mL") +
    xlab(NULL) + ylab(NULL) +
    #labs(title = paste(aSample$SamName)) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
      #axis.text = element_text(colour = "black", size = rel(1.05), face = "plain"),
      #axis.title = (element_text(colour = "black", size = rel(1.5), face = "bold")),
      strip.text = (element_text(size=rel(1.4))),
      legend.position="none"
    ) +
    theme(text=element_text(size=16,  family="Trade Gothic LT Std")) +
    #theme(text=element_text(size=16,  family="Avenir Next Condensed")) +
    #theme(plot.title = element_text(vjust = -7, hjust=0.07))+
    facet_grid(. ~ tempTitle ) +
    theme( strip.text.x = element_text( hjust = 0, family="Trade Gothic LT Std", face="bold" )) +
    guides(col = guide_legend(ncol = 1, title=NULL), alpha = FALSE)
  
  # add 10^4 to min y value to set upper range
  # ggplot_build(p)$layout$panel_params[[1]]$y.range[2] = 
  #   ggplot_build(p)$layout$panel_params[[1]]$y.range[1] + 10^4
  
  suppressWarnings( print(p) )
}



# function to plot material median results against the nominal values
pubResults2a <- function( samMedian ) {
  
  # calculate the results
  # get the data in scope
  myResults = getStudyData()
  
  # Now compute the results
  i = 0
  myResults = study
  for( Lab in levels(study$Lab)) {
    
    newLab = calAndPredict(study[which(study$Lab == Lab),], FALSE)
    if (i == 0)
      myResults = newLab
    else
      myResults = rbind(myResults, newLab)
    i = i+1
  }
  
  # get the sample medians
  samMedian = materialSummary( myResults )
  
  # here are the nominal values in alphabetical order, in log10 copies/mL
  myNominals = c(10.30103, 4, 5.195899652, 3.698970004, 4.505149978, 6.73, 3.698970004, 4.698970004)
  
  # here's a nice plottable data frame
  myCompar = cbind( samMedian, myNominals)
  names(myCompar) = c( "Material", "studyVal", "studyValCI", "nominal")
  
  # make the ggplot object
  myCompPlot = ggplot( myCompar, aes( x=10^nominal, y = 10^studyVal, color = Material))
  
  # and a pointrange (error bar) plot
  myCompPlot + 
    geom_pointrange(aes( ymin = 10^(studyVal - studyValCI), ymax = 10^(studyVal + studyValCI)), alpha =0.7) + 
    theme_bw() + 
    ylab("Study Results (IU/mL)") +
    xlab( "Nominal Value (genome copies/mL)") +
    scale_y_log10(label = label_log(), breaks=breaks_log()) +
    scale_x_log10(label = label_log(), breaks=breaks_log()) +
    theme( legend.position = "bottom") +
    theme(legend.text=element_text(size=10), legend.title=element_text(size=12)) +
    geom_abline(intercept = 0, slope = 7.7/8, color = "grey") +
    theme(text=element_text(size=16,  family="Trade Gothic LT Std"))
  
  ggsave(paste("fig2a", ".tiff", sep=""), width=180, height = 180, units ="mm", dpi=600, device="tiff" ) 
}


pubPoC2b <- function( samMedian ) {
  
  # read the data and compute the PoC experiment results
  myPoC = computePoC()
  
  # make Figure 2b for publication
  myPoC_bySample = ggplot(myPoC, aes(x=LabTargExpt, y=10^log10IU_Sample))
  
  myPoC_bySample + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_jitter(aes(color=Calibrant)) +
    scale_y_log10(label = label_log(), breaks=breaks_log()) +
    scale_colour_manual(values=c("#3182bd", "#fdae6b")) +
    facet_wrap(~SamName.x) +
    #ggtitle( "Calibrated results on Clinical Samples") +
    ylab( "IU" ) + xlab ("Lab:Target:Run") +
    theme(text=element_text(size=14,  family="Trade Gothic LT Std")) +
    theme(axis.text.x = element_text( size = 9 )) +
    theme(legend.text=element_text(size=10), legend.title=element_text(size=12)) +
    theme( strip.text.x = element_text(size=16, hjust = 0, family="Trade Gothic LT Std", face="bold" ))
  
  #ggsave( "fig2b", width=8, height = 10, units ="in", device="pdf" )
  ggsave(paste("fig2b", ".tiff", sep=""), width=180, height = 180, units ="mm", dpi=600, device="tiff" ) 
  
}
