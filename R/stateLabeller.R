#' Title
#'
#' @param states
#' @param groupLabels
#' @param nodes
#'
#' @return
#' @export
#'
#' @examples
EMTScoreCalc <- function(states, groupLabels, nodes) {
    states <- states %>% str_remove_all("'") %>%
        str_split("") %>% sapply(as.integer) %>% t %>% data.frame %>%
        set_names(nodes)
    Es <- states %>% select(groupLabels[[1]]) %>% rowSums
    Ms <- states %>% select(groupLabels[[2]]) %>% rowSums
    Es / length(groupLabels[[1]]) - (Ms / length(groupLabels[[2]]))
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
LabellerPhen <- function(x) {
    Phenotype <- rep("H", length(x))
    Phenotype[x == 1] <- "E"
    Phenotype[x == -1] <- "M"
    return(Phenotype)
}

#' Title
#'
#' @param state
#' @param inflMat
#' @param groupLabels
#'
#' @return
#' @export
#'
#' @examples
strengthCalc <- function(state, inflMat, groupLabels) {
    nodes <- colnames(inflMat)
    Estate <- ifelse(nodes %in% groupLabels[[1]], state, 0)
    Mstate <- ifelse(nodes %in% groupLabels[[2]], state, 0)
    Strength <- inflMat[state == 1, state == 1] %>% sum
    Epithelial <- inflMat[Estate == 1, Estate == 1] %>% sum
    Mesenchymal <- inflMat[Mstate == 1, Mstate == 1] %>% sum
    frust <- sapply(1:length(state), function(i) {
        sapply(1:length(state), function(j) {
            inflMat[i, j] * state[i] * state[j] < 0
        }) %>% sum
    }) %>% sum
    return(c(Strength, Epithelial, Mesenchymal, frust / (length(nodes) ^
                                                             2)))
}

#' Title
#'
#' @param topoFile
#' @param logDf
#'
#' @return
#' @export
#'
#' @examples
ScoreNStrength <- function(topoFile) {

    freqDf <- bmodelFreqFileReader(topoFile, flagFilter = F)
    if (nrow(freqDf) == 0)
        return()
    states <- freqDf$states
    nodes <- readLines(str_replace(topoFile, ".topo", "_nodes.txt"))
    if (!file.exists(str_replace(topoFile, ".topo", ".team")))
        GroupStrength(topoFile, getTeams = T, logDf = logDf)
    groupLabels <-
        readLines(str_replace(topoFile, ".topo", ".teams")) %>% str_split(",")
    inflMat <-
        read.csv(paste0(
            "Influence/",
            str_replace(topoFile, ".topo", "_fullInfl.csv")
        ),
        row.names = 1) %>% as.matrix
    inflMat <- inflMat[nodes, nodes]
    statesDf <- states %>% str_remove_all("'") %>%
        str_split("") %>% sapply(as.integer) %>% t %>% data.frame %>%
        set_names(nodes)
    scores <- apply(statesDf, 1, function(x) {
        strengthCalc(x, inflMat, groupLabels)
    }) %>% t %>% data.frame %>% set_names(c(
        "Strength",
        "Epithelial",
        "Mesenchymal",
        "InflFrustration"
    ))

    freqDf <- freqDf %>%
        mutate(EMTScore = EMTScoreCalc(states, groupLabels, nodes)) %>%
        mutate(Phenotype = LabellerPhen(EMTScore)) %>%
        mutate(
            Strength = scores$Strength,
            Epithelial = scores$Epithelial,
            Mesenchymal = scores$Mesenchymal,
            InflFrust = scores$InflFrustration
        )
    write_csv(
        freqDf,
        str_replace(topoFile, ".topo", "_finFlagFreq.csv"),
        quote = "none",
        na = ""
    )
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
ScoreNStrengthAll <- function() {
    topoFiles <- list.files(".", ".topo$")
    sapply(topoFiles, ScoreNStrength)
}
