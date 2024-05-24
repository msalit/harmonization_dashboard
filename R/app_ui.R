# RNA Harmonization Analysis App

library(shiny)
library(shinythemes)
if (!require("DT")) install.packages("DT")
library(DT)


study <- getStudyData()

# compute the results
myResults <- lapply(levels(study$Lab), function(Lab) {
    calAndPredict(study[study$Lab == Lab, ])
})
myResults <- do.call(rbind, myResults)

# get lists of labs, materials, and targets to populate dropdown UI selectors
labDropdownOptions <- levels(myResults$Lab)
materialDropdownOptions <- fct_drop(factor(myResults$SamName[which(!myResults$isCalib & !myResults$isCtrl)]))
targetDropdownOptions <- levels(fct_drop(factor(myResults[which(myResults$Lab == labDropdownOptions[1]), ]$Target)))


ui <- fluidPage(
    theme = shinytheme("united"),
    titlePanel("CSWG RNA Harmonization Study -- preliminary results", windowTitle = "CSWG RNA Harmonization Study Analysis"),
    fluidRow(
        column(5, h4("Calibration Curve & Lab Data")),
        column(7, h4("Material Results"))
    ),
    fluidRow(
        column(
            5,
            mainPanel(
                fluidRow(
                    column(
                        6,
                        selectInput("selectedLab", label = h5("Lab:"), choices = as.character(labDropdownOptions))
                    ),
                    column(
                        6,
                        selectInput("selectedTarget", label = h5("Target:"), choices = as.character(targetDropdownOptions))
                    ),
                ),
                tabsetPanel(
                    tabPanel(
                        "Calibration Curve",
                        plotOutput("calPlot")
                    ),
                    tabPanel(
                        "Calibration Fit Stats",
                        uiOutput("calSummary")
                    ),
                    tabPanel(
                        "Calibration Results",
                        DT::dataTableOutput("labTargTable")
                    )
                )
            )
        ),
        column(
            7,
            # Show a tables and plots of the materials results
            mainPanel(
                tabsetPanel(
                    tabPanel(
                        "Material Tabs",
                        selectInput("selectedMaterial", label = h5("Material:"), choices = as.character(materialDropdownOptions), selected = 1),
                        tabsetPanel(
                            tabPanel("Material x Rank", plotOutput("rankPlot")),
                            tabPanel("Material x Lab", plotOutput("labPlot")),
                            tabPanel("Deviations x Lab", plotOutput("devPlot")),
                            tabPanel("Material x Lab Results", DT::dataTableOutput("matTable")),
                            tabPanel("Material Value Results", DT::dataTableOutput("matSummary")),
                            tabPanel("Results v. Nominal", plotOutput("compPlot")),
                            tabPanel(
                                "Full Dataset",
                                tags$div(
                                    style = "padding-top:20px; padding-bottom:20px",
                                    #downloadButton("downloadFullData", "Download Full Dataset as .csv", style = "background-color: #377b57;")
                                    downloadButton("downloadFullData", "Download Full Dataset as .csv")
                                ), DT::dataTableOutput("studyComplete")
                            )
                        )
                    ),
                    tabPanel(
                        "PoC Results",
                        br(),
                        h4("Results from the Proof of Concept Study"),
                        tags$div(
                            style = "padding-top:20px; padding-bottom:20px",
                            #downloadButton("downloadPoCData", "Download PoC Results as .csv", style = "background-color: #2861bd;")
                            downloadButton("downloadPoCData", "Download PoC Results as .csv")
                        ),
                        tabsetPanel(
                            tabPanel("PoC Shared Y Axis", plotOutput("POC_results_shared")),
                            tabPanel("PoC Free Y Axis", plotOutput("POC_results_free")),
                        )
                    )
                )
            )
        )
    )
)
