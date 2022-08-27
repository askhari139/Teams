### Load the package : devtools::install_github("askhari139/Teams)
### Or, download the repo and load the individual script files as follows:

sourceLocation <- "D:/Teams/R"
setwd(sourceLocation)
source("setupScript.R")
source("influenceAndGs.R")
source("topoGenAndSim.R")
source("coherence.R")
source("expressionCorrelation.R")
source("stateLabeller.R")
source("utils.R")
source("compileData.R")




library(ggrepel)
library(grid)
library(gridExtra)

### Setup data folder structure
mainFolder <- "D:/Teams/Simulations" # folder where all the data is stored.
SetupFunc(mainFolder = mainFolder,
          topoFolder = paste0(mainFolder, "/TopoFiles"), # folder containing topofiles to analyze
          numThreads = 7, # number of threads to be used for simulations and coherence calculations
          julia = F) # Option to specify whether to install the Bmodel Julia package, for simulations

source(paste0(cwd, "/plotUtils.R"))
source(paste0(cwd, "/wtPlotUtils.R"))
