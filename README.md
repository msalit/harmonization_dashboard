# RNA Harmonization Analysis App

Welcome to the RNA Harmonization Study Analysis repository!

## Links

- The live dashboard is available !TODO.
- The paper for this study is available !TODO.

## Organization

- All data files from the study are under `/data`
  - 20210603 MassCPR.csv -- data from the tech dev platform at MassCPR
  - 20210603 NIST-NML-NIB.csv -- data from National Labs, different experiment design
  - 20210810_normal_labs.csv -- data from other labs, formatted consistently
  - 20230223 PoC.csv -- Proof of Concept experiment data
  
- Source code is in `/R`, which is organized into folders according to purpose.
  
- The app is composed from a [shiny](https://shiny.posit.co/) server and ui, specified in [app_server.R](R/app_server.R) and [app_ui.R](R/app_ui.R) respectively.

- Calibration functions are specified in [calibration_utils.R](R/calibration_utils.R):
  - `calAndPredict()`
  - `fitCal()`
  - `fitGroup()`
  - `calcResults()`
- Data i/o functions are specified in [data_intake.R](R/data_intake.R)
  - `getPoC()`: reads in the Proof of Concept data from a .csv.
  - `getStudyData()`: combines data read in from the rest of the .csvs associated with the study.
- Functions which plot visualizations from the study data are specified in [plotting_utils.R](R/plotting_utils.R)
  - `plotCal()`
  - `sampleByRank()`
  - `sampleByLabTarg()`
  - `comparePlot()`
- Functions which create tables from the study data are specified in [table_utils.R](R/table_utils.R)
  - `materialTable()`
  - `materialSummary()`
  - `labTargTable()`
- And finally, functions which process the Proof of Concept data reside in [poc_utils.R](R/poc_utils.R)
  - `plot_poc_results()`
  - `calibPoC()`
  - `pubPlots()`
