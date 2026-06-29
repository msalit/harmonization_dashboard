# RNA Harmonization Analysis App

Welcome to the RNA Harmonization Study Analysis repository!

## Links

- The live dashboard is available at https://msalit.shinyapps.io/HarmonizationDashboard/.
- The paper for this study, *A playbook for harmonization of standards for emergent viral pathogens*, is available as a medRxiv preprint at https://doi.org/10.1101/2024.10.05.24314949 (v1: https://www.medrxiv.org/content/10.1101/2024.10.05.24314949v1).

## Running the app locally

The app entry point is [app.R](app.R) at the project root. From an R session in the project root:

```r
shiny::runApp()
```

Or open the project in RStudio and click **Run App** with `app.R` open.

### Required packages

All from CRAN. Install once:

```r
install.packages(c(
  "shiny", "shinythemes", "DT",
  "tidyverse",
  "scales", "jtools", "huxtable"
))
```

`tidyverse` covers `ggplot2`, `dplyr`, `forcats`, `readr`, `tidyr`, and `stringr`, all of which the app uses.

### About `renv.lock.archive`

`renv.lock.archive` is the exact package environment used at manuscript submission (R 4.4.0 + pinned package versions). It is **not** activated by this project — installing the packages above against current R is the supported path. The archive is kept only so a curious reader (or reviewer) can reconstruct the submission environment with `renv::restore(lockfile = "renv.lock.archive")` if they want to.

## Organization

- All data files from the study are under [`data/`](data/)
  - [20210603_MassCPR.csv](data/20210603_MassCPR.csv) -- data from the tech dev platform at MassCPR
  - [20210603_NIST-NML-NIB.csv](data/20210603_NIST-NML-NIB.csv) -- data from National Labs, different experiment design
  - [20210810_normal_labs.csv](data/20210810_normal_labs.csv) -- data from other labs, formatted consistently
  - [20230223_PoC.csv](data/20230223_PoC.csv) -- Proof of Concept experiment data

- Source code is in [`R/`](R/), one file per concern:

- The app is composed from a [shiny](https://shiny.posit.co/) server and ui, specified in [app_server.R](R/app_server.R) and [app_ui.R](R/app_ui.R) respectively.

- Calibration functions are specified in [calibration_utils.R](R/calibration_utils.R):
  - `calAndPredict()`
  - `fitCal()`
  - `fitGroup()`
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
- Functions which process the Proof of Concept data reside in [poc_utils.R](R/poc_utils.R)
  - `plot_poc_results()`
  - `computePoC()`
- And finally, functions to create and save publication-ready plots are in [makePubPlots.R](R/makePubPlots.R)
  - `pubBoxPlots1b()`
  - `pubResults2a()`
  - `pubPoC2b()`

## Adapting this dashboard for a new study

This app is **exemplary, not general-purpose**: it is wired to the specific materials, calibrants, replication structure, and reference values of the SARS-CoV-2 RNA Harmonization Study. Reusing the codebase for a new harmonization study requires a deliberate sweep through the study-specific points below.

1. **Replace the data files in [`data/`](data/)** with your study's CSVs. The schemas the code expects are:
   - **Wide-format lab data** (e.g. `20210810_normal_labs.csv`, `20210603_NIST-NML-NIB.csv`): columns `Lab, Target, isCtrl, Batch, dPCR, Sample code, SamName, Dilution, R1, R2, R3, R4`. The replicate columns are pivoted longer in [data_intake.R](R/data_intake.R) by index (`cols = 9:12`).
   - **Long-format lab data** (e.g. `20210603_MassCPR.csv`): columns `Lab, Target, isCtrl, Batch, dPCR, Sample code, SamName, Dilution, Rep, Sig`.
   - **PoC data** (`20230223_PoC.csv`): columns `Lab, Target, SamName, isCtrl, isCalib, Replicate, dPCR, Raw.Value, Experiment`.

2. **Update calibrant identifiers and concentrations** in [data_intake.R](R/data_intake.R) lines 65–73. The current code hardcodes the WHO International Standard dilution series (`WHO-IS_dilution_0` … `WHO-IS_dilution_6`) and assigns concentrations of `10^7.7, 10^6.7, …, 10^1.7 IU/mL`. Change both the names (used by `str_detect` on line 76 to flag calibrants) and the values for your study.

3. **Audit lab-specific special cases.** [data_intake.R](R/data_intake.R) lines 78–80 halve the calibrant concentration for the `MUSC` lab because they rehydrated the WHO IS in 1.0 mL instead of 0.5 mL. Anything similar in your dataset needs its own conditional.

4. **Update PoC calibrant IU values** in [poc_utils.R](R/poc_utils.R) lines 63–67. The current code hardcodes `AcroMetrix` (`log10IU = 4.033`) and `LGC` (`log10IU = 3.911`) with their CIs. The `Asuragen` filter on line 38 (`subset(PoC, SamName != "Asuragen")`) drops one sample because it was uncalibrated; revisit that filter if your PoC has a different uncalibrated set.

5. **Update the reference (nominal) values for materials.** Both [plotting_utils.R](R/plotting_utils.R) (line 225, in `comparePlot()`) and [makePubPlots.R](R/makePubPlots.R) (in `pubResults2a()`) define `myNominals` as a length-8 numeric literal **in alphabetical order of `SamName`**. These are the manuscript's fixed truth values for the 8 study materials, traceable to the WHO candidate-material preparation described in [WHO/BS/2020.2403](#references). For a new study, update both values and length, and verify the alphabetical ordering matches `materialSummary()`'s output.

6. **Update the replicate column count.** [data_intake.R](R/data_intake.R) lines 25–28 use `pivot_longer(cols = c(9:12))` to fold four replicate columns (`R1`–`R4`) into long form. If your new wide-format CSVs have a different number of replicate columns, change the column indices.

7. **Update the factor relevel index list** at [data_intake.R](R/data_intake.R) line 106:
   ```r
   allData_l$SamName <- fct_drop(factor(..., levels = levels(allData_l$SamName)[c(9:15, 1:8, 16)]))
   ```
   The hardcoded `c(9:15, 1:8, 16)` puts the WHO-IS dilution series first (positions 9–15 in the alphabetical factor), then the eight study materials (1–8), then the leftover (16). For a different sample roster, the index list must be recomputed by hand.

8. **Update study-specific UI text** in [app_ui.R](R/app_ui.R) — `titlePanel("CSWG RNA Harmonization Study -- preliminary results", …)` (line 25) and any other prose that names this study.

9. **Replace the publication-plot font.** [makePubPlots.R](R/makePubPlots.R) uses `family = "Trade Gothic LT Std"` in several `theme(text = element_text(...))` calls. Substitute a font available on your system, or comment those theme lines out to use ggplot2 defaults.

10. **Replace `myNominals`'s nominal-vs-IU regression line.** [plotting_utils.R](R/plotting_utils.R) line 244 draws `geom_abline(intercept = 0, slope = 7.7/8, color = "grey")`, which encodes the reference relationship between log10 genome-copies/mL and log10 IU/mL specific to the WHO IS. The slope `7.7/8` is `7.70` (the assigned potency, in log10 IU/mL, of the WHO IS after reconstitution in 0.5 mL) divided by `8` (the log10 of the `1×10⁸` genome-copies/mL concentration the bulk candidate material was prepared to contain) — both values from [WHO/BS/2020.2403](#references). Update the slope (and possibly intercept) for your study's reference standard.

## References

The nominal (reference) values and the genome-copies/mL ↔ IU/mL relationship used in `comparePlot()` ([plotting_utils.R](R/plotting_utils.R)) and `pubResults2a()` ([makePubPlots.R](R/makePubPlots.R)) derive from the WHO collaborative study that established the International Standard. The material-preparation details quoted in [data_intake.R](R/data_intake.R) (the `1×10⁸` genome-copies/mL bulk concentration and the `7.70` log10 IU/mL assigned potency) are from page 5 of that report:

- Bentley E, Mee ET, Routley S, Mate R, Fritzsche M, Hurley M, Le Duff Y, Anderson R, Hockley J, Rigsby P, Page M, Rose N, Mattiuzzo G, and the Collaborative Study Group. *Collaborative Study for the Establishment of a WHO International Standard for SARS-CoV-2 RNA.* WHO Expert Committee on Biological Standardization, Geneva, 9–10 December 2020. **WHO/BS/2020.2403**.
