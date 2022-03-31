hamming <- function(x)
{
    x <- str_split(x, "")
    sum(x[[1]]!=x[[2]])/(length(x[[1]])-2)
}
hamming <- cmpfun(hamming)

dec2bin <- function(x, l=32)
{
    paste(rev(rev(as.integer(rev(intToBits(x))))[1:l]), collapse = "")
} %>% cmpfun


state_gen <- function(n){
    s_list <- 1:(2^n)
    sapply(s_list, dec2bin, n)
} %>% cmpfun

inwards <- function(x,df){
    colnames(df) <- c("State", "Step")
    df %>% filter(Step == x) %>% select(State) %>% unlist
} %>% cmpfun

asyncInit <- function(sInit, update_matrix, maxT = 1000, init = NULL)
{#browser()
    f <- 0 #flag
    N_nodes <- length(sInit)
    up_ind_list <- sample(1:N_nodes, maxT, replace = T) # generate the update indices
    s <- sInit
    for (k in 1:maxT){
        s_dummy <- s%*%update_matrix %>% sign # update the state
        if (all(s_dummy == s))
            f <- 1 # flag 1 means it is a proper steady state

        up_ind <- up_ind_list[k]
        s[up_ind] <- s_dummy[up_ind]

        if (f == 1)
            break
    }
    if (is.null(init))
        init <- paste0("'",
                       paste0(
                           replace(sInit, which(s==-1),0),
                           collapse = ""),
                       "'")
    fin <- paste0("'",
                  paste0(
                      replace(s, which(s==-1),0),
                      collapse = ""),
                  "'")
    return(c(init, fin, f))
}

asyncInit <- cmpfun(asyncInit)

coherenceIter <- function(sInit, update_matrix, maxT = 1000, nNode = 1,
                          randChoice = F, nState = 100, nIter = 10)
{#browser()
    init <- sInit
    sInit <- str_split(sInit, "") %>% unlist
    sInit <- as.integer(sInit[2:(length(sInit)-1)])
    sInit <- ifelse(sInit == 0, -1, 1)
    if (nNode>1)
        randChoice <- T
    nNodes <- length(sInit)
    pertChoice <- sapply(1:nState, function(x){
        sample(1:nNodes, nNode)
    })
    if (nNode > 1)
        pertChoice <- t(pertChoice)
    else
        pertChoice <- matrix(pertChoice, ncol = 1)

    stateList <- apply(pertChoice, 1, function(x){
        s <- sInit
        s[x] <- -1*s[x]
        s
    })

    if(!randChoice)
        stateList <- sapply(1:nNodes, function(x){
            s <- sInit
            s[x] <- -1*s[x]
            s
        })
    # browser()
    df <- apply(stateList,2, function(x){
        sapply(1:nIter, function(i){
            # browser()
            asyncInit(x, update_matrix, maxT, init)
        }) %>% t %>% data.frame %>% set_names(c("init", "fin", "flag"))
    }, simplify = F) %>% reduce(rbind.data.frame) %>%
        group_by(init, fin, flag) %>%
        summarise(Freq = n(), .groups = 'drop') %>% ungroup
    df$Freq <- round(df$Freq/sum(df$Freq),4)
    df

}
coherenceIter <- cmpfun(coherenceIter)


coherence <- function(topoFile, maxT = 1000, nNode = 1,
                      randChoice = F, nState = 100, nIter = 10, write=T)
{
    net <- topoFile %>% str_remove(".topo")
    nodeOrder <- readLines(paste0(net, "_nodes.txt"))
    ls <- TopoToIntMat(topoFile)
    update_matrix <- ls[[1]]
    nodes <- ls[[2]]
    colnames(update_matrix) <- rownames(update_matrix) <- nodes
    update_matrix <- update_matrix[nodeOrder, nodeOrder]
    update_matrix <- 2*update_matrix + diag(length(nodes))
    groupLabels <- readLines(str_replace(topoFile, ".topo", ".teams")) %>% str_split(",")
    statesDf <- read_csv(paste0(net, "_finFlagFreq.csv"), col_types = cols(), lazy = F) %>%
        filter(!is.na(Avg0), flag == 1)
    if (nrow(statesDf) == 0)
        return(NA)
    states <- statesDf$states
    # browser()
    df <- lapply(states, function(x){
        coherenceIter(x, update_matrix, maxT,nNode,
                      randChoice, nState, nIter)
    }) %>% reduce(rbind.data.frame)
    df$flag <- as.integer(df$flag)
    df$nNode <- nNode
    df$Hamming <- df %>% select(init, fin) %>% apply(1, hamming)
    df$initPhen <- df$init %>% EMTScoreCalc(groupLabels, nodes) %>% LabellerPhen
    df$finPhen <- df$fin %>% EMTScoreCalc(groupLabels, nodes) %>% LabellerPhen
    if (write) {
        directoryNav("CoherencesData")
        write_csv(df, paste0(net, "_", nNode, "R_coherence.csv"),
                  quote = "none")
    }

    else
        return(df)
}

coherence <- cmpfun(coherence)

CoherenceSingleNode <- function()
{
    logDf <- read.csv("LogFile.csv")
    topoFiles <- logDf %>% filter(Coherence == "No") %>% select(Files) %>% unlist
    nets <- topoFiles %>% str_remove(".topo")
    plan(multisession, workers = numThreads)
    dummy <- future_lapply(nets, function(x){
        # browser()
        df <- coherence(paste0(x, ".topo"), randChoice = F, nIter = 50, write = F)
        df <- df[df$init == df$fin,] %>% mutate(states = init, coherence0 = Freq) %>%
            select(states, coherence0)
        coherenceVec<- df$coherence0
        names(coherenceVec) <- df$states
        df <- read_csv(paste0(x, "_finFlagFreq.csv"), show_col_types = F, lazy=F) %>%
            mutate(coherence0 = coherenceVec[states])
        write_csv(df, paste0(x, "_finFlagFreq.csv"), quote = "none", na = "")
        logDf <<- logDf %>% mutate(Coherence = ifelse(Files == paste0(x, ".topo"), "Yes", Coherence))
    })
    future:::ClusterRegistry("stop")
    write_csv(logDf, "LogFile.csv", quote = "none")
}
CoherenceSingleNode <- cmpfun(CoherenceSingleNode)



CoherenceAllNode <- function(topoFile, logDf = NULL, maxT = 1000, randChoice = T,
                             nState = 100, nIter = 10, reps = 1)
{
    writeLog <- F
    if (is.null(logDf)) {
        writeLog <- T
        logDf <- read_csv("LogFile.csv", col_types = cols(), lazy = F)
    }

    net <- topoFile %>% str_remove(".topo")
    l <- readLines(paste0(net, "_nodes.txt")) %>% length
    # browser()
    dfList <- lapply(1:reps, function(rep){
        plan(multisession, workers = numThreads)
        df <- lapply(1:l, function(x){
            coherence(topoFile, maxT, x, randChoice, nState, nIter, write = F)
        }) %>% reduce(rbind.data.frame) %>% mutate(Rep = rep)
        future:::ClusterRegistry("stop")
        return(df)
    })
    if (reps == 1)
        df <- dfList[[1]]
    else
        df <- dfList %>% reduce(rbind.data.frame)
    DirectoryNav("PhenotypicTransition")
    write_csv(df, paste0(net, "_allNodeCoherence_nPert",
                         nState,"_nIter",nIter,"_reps", reps, ".csv"), quote = "none")
    logDf <- logDf %>%
        mutate(MultiNodeCoherence = ifelse(Files == paste0(x, ".topo"), "Yes", MultiNodeCoherence))
    setwd("..")
    if(writeLog) {
        write_csv(logDf, "LogFile.csv", quote = "none")
    }
    else
        logDf <<- logDf
}
CoherenceAllNode <- cmpfun(CoherenceAllNode)
