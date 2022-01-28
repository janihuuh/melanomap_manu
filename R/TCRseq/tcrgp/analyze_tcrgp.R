
## An analysis file for all TCRGP-related work
source("src/jani/R/tcrb/main.R")

source("src/jani/R/tcrb/tcrgp/fun_tcrgp.R")
source("src/jani/R/tcrb/tcrgp/fun_plotTcrgp.R")



## First, do summarisation from the TCRGP-output. These are done on two batches, which we need to merge
source("src/jani/R/tcrb/tcrgp/run_tcrgpSummary.R")
source("src/jani/R/tcrb/tcrgp/run_tcrgpMerge.R")


## Analyze TCRGP from baseline; evolution and expanded
source("src/jani/R/tcrb/tcrgp/plot_tcrgpBaseline.R")
source("src/jani/R/tcrb/tcrgp/plot_tcrgpEvolution.R")
source("src/jani/R/tcrb/tcrgp/plot_tcrgpEvolution2.R")
