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
 


### Setup data folder structure
mainFolder <- "D:/Teams/Simulations" # folder where all the data is stored.
SetupFunc(mainFolder = mainFolder,
          topoFolder = paste0(mainFolder, "/TopoFiles"), # folder containing topofiles to analyze
          numThreads = 7, # number of threads to be used for simuations and coherece calculations
          julia = F) # Option to specify whether to install the Bmodel Julia package, for simulations


### Random network Analysis ----

# Generating random networks

RandomNetworks(numRand = 500) # Generate random networks for all files in topoFileFolder


# Simulating random networks
SimulateBoolean(rand = T) # Simulate all random networks

## Generating team strength, correlations, coherence, state strength and emt score data

setwd(randRaw)
sapply(netList, function(net) {
  setwd(net)
  GroupStrengthAll(net) # influence matrix and group strength
  topoFiles <- list.files(".", ".topo$")
  sapply(topoFiles, correlationMatBool) # correlation matrices for all random networks
  CoherenceSingleNode() # Single node coherence for all steady states in all WT and random networks
  ScoreNStrengthAll() # State phenotypes, EMTscores, state strengths and other misc. calculations.
  setwd("..")
})
setwd(randRaw)
sapply(EMPNets, function(net) {
  setwd(net)
  ScoreNStrengthAll()
  AllDataFile(net) # Compile all metric summary, only for pure steady states
  AllDataFileNoFlag(net) # compile summary for all states
  setwd("..")
})
# Multinode Perturbation
setwd(randRaw)
sapply(netList, function(net) {
  setwd(net)
  topoFiles <- list.files(".", ".topo$")
  topoDone <- ""
  if (dir.exists("PhenotypicTransition"))
  {
    setwd("PhenotypicTransition")
    topoDone <- list.files(".", "_allNodeCoherence") %>%
      str_remove("_allNode.*") %>% paste0(".topo")
    setwd("..")
  }
  t <- which(topoFiles %in% topoDone)
  if (length(t) != 0)
    topoFiles <- topoFiles[-which(topoFiles %in% topoDone)]
  plan(multisession, workers = numThreads)
  future_sapply(topoFiles, CoherenceAllNode)
  future:::ClusterRegistry("stop")
  setwd("..")
})

## Signal Node perturbation ---
sapply(netList, coherenceSignal)


### Causation network Analysis ----

# Generating causation networks

CausationNetworks(nMax = 20)

# Simulating causation networks

SimulateBoolean(cause = T)

# Generating team strength, coherence, state strength and emt score data

setwd(gsCausation)
sapply(netList, function(net) {
  setwd(net)
  LogFileGen()
  GroupStrengthAll(net)
  topoFiles <- list.files(".", ".topo$")
  sapply(topoFiles, correlationMatBool, logDf = logDf)
  CoherenceSingleNode()
  ScoreNStrengthAll()
  AllDataFile(net)
  AllDataFileNoFlag(net)
  setwd("..")
})

# # Multinode Perturbation
setwd(gsCausation)
sapply(netList, function(net) {
  setwd(net)
  topoFiles <- list.files(".", ".topo$")
  logDf <- read.csv("LogFile.csv")
  sapply(topoFiles, CoherenceAllNode, logDf = logDf)
  write.csv(logDf, "LogFile.csv", row.names = F)
  setwd("..")
})


