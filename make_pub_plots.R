# make_pub_plots.R
# ----------------------------------------------------------------------------
# Standalone driver to regenerate the publication figures from RStudio.
#
# Open RNAStudy_for_Publication.Rproj first (so the working directory is the
# project root), then run:   source("make_pub_plots.R")
# ...or step through it interactively after the source() calls below.
#
# Outputs (written to the project root; 180 x 180 mm, 600 dpi TIFF):
#   <Material>.tiff        -- Fig 1b: one rank-ordered box plot per material
#   fig2a.tiff             -- Fig 2a: study IU vs nominal copies (free axes)
#   fig2a_equalaxes.tiff   -- Fig 2a: identical x/y scale range (reviewer ask)
#   fig2b.tiff             -- Fig 2b: Proof-of-Concept clinical-sample results
#
# NOTE: the figures use the "Trade Gothic LT Std" font. If it isn't installed
# on this machine, edit the `family = ...` lines in R/makePubPlots.R (see README
# adaptation item #9) before running, or the TIFFs will warn / fall back.
# ----------------------------------------------------------------------------

# fail early if we're not at the project root (data/ and R/ must be visible)
if (!dir.exists("data") || !dir.exists("R")) {
  stop("Working directory is not the project root. Open ",
       "RNAStudy_for_Publication.Rproj (or setwd() to the repo root) first.")
}

# required packages -- only what the publication plots actually use
library(tidyverse)   # ggplot2, dplyr, forcats, readr, tidyr, stringr
library(scales)      # label_log(), breaks_log(), trans_breaks(), math_format()

# source the utility code, in dependency order, into the global environment
source("R/data_intake.R")        # getStudyData(), getPoC()
source("R/calibration_utils.R")  # calAndPredict(), fitCal(), fitGroup()
source("R/table_utils.R")        # materialSummary()  (used by pubResults2a)
source("R/plotting_utils.R")     # shared plotting helpers
source("R/poc_utils.R")          # computePoC()       (used by pubPoC2b)
source("R/makePubPlots.R")       # pubBoxPlots1b(), pubResults2a(), pubPoC2b()

# --- generate the figures (comment out any you don't need to rebuild) --------
pubBoxPlots1b()   # Fig 1b: one <Material>.tiff per study material
pubResults2a()    # Fig 2a: fig2a.tiff (free axes) + fig2a_equalaxes.tiff
pubPoC2b()        # Fig 2b: fig2b.tiff

message("Done. TIFFs written to: ", normalizePath(getwd()))
