# RNA Harmonization Study Analysis App

# source the libraries
require(tidyverse)
require(ggplot2)
require(plotly)
require(dplyr)
require(broom)
require(forcats)
require(rmarkdown)
require(sjPlot)
require(knitr)
require(shiny)
require(robustbase)
require(GGally)

# app components
source("app_ui.R")
source("app_server.R")

# utilities
source("data_intake.R")
source("calibration_utils.R")
source("plotting_utils.R")
source("table_utils.R")
source("poc_utils.R")
source("makePubPlots.R")


# Run the shiny app
shinyApp(ui = ui, server = server)
