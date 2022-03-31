EMTScoreCalc <- function(states, groupLabels, nodes) {
    states <- states %>% str_remove_all("'") %>%
        str_split("") %>% sapply(as.integer) %>% t %>% data.frame %>%
        set_names(nodes)
    Es <- states %>% select(groupLabels[[1]]) %>% rowSums
    Ms <- states %>% select(groupLabels[[2]]) %>% rowSums
    Es/length(groupLabels[[1]]) - (Ms/length(groupLabels[[2]]))
}

LabellerPhen <- function(x) {
    Phenotype <- rep("H", length(x))
    Phenotype[x == 1] <- "E"
    Phenotype[x == -1] <- "M"
    return(Phenotype)
}

strengthCalc <- function(state, inflMat, groupLabels) {
    nodes <- colnames(inflMat)
    Estate <- ifelse(nodes %in% groupLabels[[1]], state, 0)
    Mstate <- ifelse(nodes %in% groupLabels[[2]], state, 0)
    Strength <- inflMat[state == 1, state == 1] %>% sum
    Epithelial <- inflMat[Estate == 1, Estate == 1] %>% sum
    Mesenchymal <- inflMat[Mstate == 1, Mstate == 1] %>% sum
    return(c(Strength, Epithelial, Mesenchymal))
}

ScoreNStrength <- function(topoFile, logDf = NULL) {
    writeLog <- F
    if (is.null(logDf)) {
        writeLog <- T
        logDf <- read_csv("LogFile.csv", col_types = cols(), lazy = F)
    }

    freqDf <- bmodelFreqFileReader(topoFile)
    states <- freqDf$states
    nodes <- readLines(str_replace(topoFile, ".topo", "_nodes.txt"))
    if (!file.exists(str_replace(topoFile, ".topo", ".team")))
        GroupStrength(topoFile, getTeams = T, logDf = logDf)
    groupLabels <- readLines(str_replace(topoFile, ".topo", ".teams")) %>% str_split(",")
    inflMat <- read.csv(paste0("Influence/", str_replace(topoFile,".topo", "_fullInfl.csv")),
                        row.names = 1) %>% as.matrix
    inflMat <- inflMat[nodes, nodes]
    statesDf <- states %>% str_remove_all("'") %>%
        str_split("") %>% sapply(as.integer) %>% t %>% data.frame %>%
        set_names(nodes)
    scores <- apply(statesDf, 1, function(x){
        strengthCalc(x, inflMat, groupLabels)
    }) %>% t %>% data.frame %>% set_names(c("Strength", "Epithelial", "Mesenchymal"))

    freqDf <- freqDf %>%
        mutate(EMTScore = EMTScoreCalc(states, groupLabels, nodes)) %>%
        mutate(Phenotype = LabellerPhen(EMTScore)) %>%
        mutate(Strength = scores$Strength, Epithelial = scores$Epithelial,
               Mesenchymal = scores$Mesenchymal)
    write_csv(freqDf, str_replace(topoFile, ".topo", "_finFlagFreq.csv"), quote = "none", na = "")
    logDf <- logDf %>% mutate(Scores = ifelse(Files == topoFile, "Yes", Scores))
    if(writeLog) {
        write_csv(logDf, "LogFile.csv", quote = "none")
    }
    else
        logDf <<- logDf
}

ScoreNStrengthAll <- function() {
    logDf <- read.csv("LogFile.csv")
    topoFiles <- list.files(".", ".topo$")
    sapply(topoFiles, ScoreNStrength, logDf = logDf)
    write_csv(logDf, "LogFile.csv", quote = "none")
}
