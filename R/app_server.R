# Server function for the shiny app

library(jtools)
if (!require("DT")) install.packages("DT")
library(DT)

server <- function(input, output) {
    observe({
        # update the options for the target dropdown given the selected lab
        updateSelectInput(
            session = getDefaultReactiveDomain(), inputId = "selectedTarget",
            choices = as.character(levels(fct_drop(factor(myResults[which(myResults$Lab == input$selectedLab), ]$Target))))
        )
        output$calPlot <- renderPlot(
            width = 500, height = "auto",
            # generate cal curve based on input$selectedLab and input$selectedTarget from ui.R
            plotCal(myResults[which((myResults$Lab == input$selectedLab) & (myResults$Target == input$selectedTarget)), ])
        )

        statistics_to_display <- c(N = "nobs", R2 = "r.squared", Adj_R2 = "adj.r.squared", Residual_Std_Err = "sigma", F_statistic = "statistic")
        output$calSummary <- renderUI(
            HTML(huxtable::to_html(export_summs(summ(fitCal(myResults[which((myResults$Lab == input$selectedLab) & (myResults$Target == input$selectedTarget)), ]), digits = 2), statistics = statistics_to_display)))
        )

        output$labTargTable <- DT::renderDataTable({
            DT::datatable(
                labTargTable(
                    myResults[which((myResults$Lab == input$selectedLab) &
                        (myResults$Target == input$selectedTarget)), ]
                )[, 1:7],
                caption = "Calibration Results",
                rownames = FALSE,
                options = list(pageLength = 25)
            ) %>% DT::formatSignif(c("Signal", "log10_IU_per_mL_pred", "Wts"), digits = 4, interval = 0)
        })

        output$rankPlot <- renderPlot(
            width = 700, height = 700,
            sampleByRank(myResults[which(myResults$SamName == input$selectedMaterial), ])
        )

        output$labPlot <- renderPlot(
            width = 700, height = 700,
            sampleByLabTarg(myResults[which(myResults$SamName == input$selectedMaterial), ])
        )

        output$devPlot <- renderPlot(
            width = 700, height = 700,
            plotDeviations(myResults)
        )

        output$compPlot <- renderPlot(
            width = 500, height = 500,
            comparePlot(materialSummary(myResults))
        )

        output$matTable <- DT::renderDataTable({
            DT::datatable(materialTable(myResults),
                caption = "Median log10 IU/mL Values",
                rownames = FALSE,
                options = list(
                    pageLength = 23,
                    lengthChange = FALSE,
                    paging = FALSE
                )
            )
        })

        output$matSummary <- DT::renderDataTable({
            DT::datatable(materialSummary(myResults),
                colnames = c("Material", "Median log10 IU/mL", "95% CI"),
                caption = "Median log10 IU/mL Values & 95% CI",
                rownames = FALSE,
                options = list(
                    pageLength = 8,
                    lengthChange = FALSE,
                    paging = FALSE
                )
            )
        })
        output$studyComplete <- DT::renderDataTable({
            DT::datatable(myResults, rownames = FALSE) %>% DT::formatSignif(c(10, 11, 13:16), digits = 3, interval = 0)
        })
        output$downloadFullData <- downloadHandler(
            filename = "harmonization_study_full_dataset.csv",
            content = function(file) {
                write.csv(myResults, file)
            }
        )
        output$downloadPoCData <- downloadHandler(
            filename = "proof_of_concept_dataset.csv",
            content = function(file) {
                file.copy(
                    from = "..\\data\\20230223 PoC.csv",
                    to = file
                )
            }
        )
        output$POC_results_shared <- renderPlot(
            width = 1000, height = 1000,
            plot_poc_results(shared_axis = TRUE)
        )
        output$POC_results_free <- renderPlot(
            width = 1000, height = 1000,
            plot_poc_results(shared_axis = FALSE)
        )
    })
}
