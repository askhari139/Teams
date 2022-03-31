### Load the package : devtools::install_github("askhari139/Teams)

### Random network Analysis ----

# Generating random networks

RandomNetworks(numRand = 10)


# Simulating random networks

SimulateBoolean(rand = T)

# Generating team strength, correlations, coherence, state strength and emt score data

setwd(randRaw)
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

# Multinode Perturbation
setwd(randRaw)
sapply(netList, function(net) {
    setwd(net)
    topoFiles <- list.files(".", ".topo$")
    logDf <- read.csv("LogFile.csv")
    sapply(topoFiles, CoherenceAllNode, logDf = logDf)
    write.csv(logDf, "LogFile.csv", row.names = F)
    setwd("..")
})


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

# Multinode Perturbation
setwd(gsCausation)
sapply(netList, function(net) {
    setwd(net)
    topoFiles <- list.files(".", ".topo$")
    logDf <- read.csv("LogFile.csv")
    sapply(topoFiles, CoherenceAllNode, logDf = logDf)
    write.csv(logDf, "LogFile.csv", row.names = F)
    setwd("..")
})


