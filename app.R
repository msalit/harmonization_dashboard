# RNA Harmonization Study Analysis App
# Single entry point at the project root. Shiny's R/ auto-loader is disabled
# (see R/_disable_autoload.R) so the source order below is deterministic.

# required packages (all CRAN)
library(shiny)
library(shinythemes)
library(DT)
library(tidyverse)   # ggplot2, dplyr, forcats, readr, tidyr, stringr
library(scales)
library(jtools)
library(huxtable)

# utilities
source(local = TRUE, "R/data_intake.R")
source(local = TRUE, "R/calibration_utils.R")
source(local = TRUE, "R/plotting_utils.R")
source(local = TRUE, "R/table_utils.R")
source(local = TRUE, "R/poc_utils.R")
source(local = TRUE, "R/makePubPlots.R")

# precompute the calibrated results once at app start
study <- getStudyData()
myResults <- do.call(rbind, lapply(levels(study$Lab), function(Lab) {
    calAndPredict(study[study$Lab == Lab, ])
}))

# dropdown options for the UI selectors
labDropdownOptions <- levels(myResults$Lab)
materialDropdownOptions <- fct_drop(factor(myResults$SamName[which(!myResults$isCalib & !myResults$isCtrl)]))
targetDropdownOptions <- levels(fct_drop(factor(myResults[which(myResults$Lab == labDropdownOptions[1]), ]$Target)))

# ui and server
source(local = TRUE, "R/app_ui.R")
source(local = TRUE, "R/app_server.R")

shinyApp(ui = ui, server = server)
