#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
hamming <- function(x) {
    x <- str_split(x, "")
    sum(x[[1]] != x[[2]]) / (length(x[[1]]) - 2)
}
hamming <- cmpfun(hamming)

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
diffPos <- function(x) {
    x <- str_split(x, "")
    which(x[[1]] != x[[2]])
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
stateCompare <- function(x) {
    x <- str_split(x, "")
    p <- x[[1]] == x[[2]]
    p[2:(length(x[[1]]) - 1)]
}
stateCompare <- cmpfun(stateCompare)

#' Title
#'
#' @param x
#' @param l
#'
#' @return
#' @export
#'
#' @examples
dec2bin <- function(x, l = 32) {
    paste(rev(rev(as.integer(
        rev(intToBits(x))
    ))[1:l]), collapse = "")
} %>% cmpfun


#' Title
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
state_gen <- function(n) {
    s_list <- 1:(2 ^ n)
    sapply(s_list, dec2bin, n)
} %>% cmpfun

#' Title
#'
#' @param x
#' @param df
#'
#' @return
#' @export
#'
#' @examples
inwards <- function(x, df) {
    colnames(df) <- c("State", "Step")
    df %>% filter(Step == x) %>% select(State) %>% unlist
} %>% cmpfun

#' Title
#'
#' @param sInit
#' @param update_matrix
#' @param maxT
#' @param init
#'
#' @return
#' @export
#'
#' @examples
asyncInit <-
    function(sInit,
             update_matrix,
             maxT = 1000,
             init = NULL)
    {
        f <- 0 #flag
        N_nodes <- length(sInit)
        up_ind_list <-
            sample(1:N_nodes, maxT, replace = T) # generate the update indices
        s <- sInit
        for (k in 1:maxT) {
            s_dummy <- s %*% update_matrix %>% sign # update the state
            if (all(s_dummy == s))
                f <- 1 # flag 1 means it is a proper steady state

            up_ind <- up_ind_list[k]
            s[up_ind] <- s_dummy[up_ind]

            if (f == 1)
                break
        }
        if (is.null(init))
            init <- paste0("'",
                           paste0(replace(sInit, which(sInit == -1), 0),
                                  collapse = ""),
                           "'")
        fin <- paste0("'",
                      paste0(replace(s, which(s == -1), 0),
                             collapse = ""),
                      "'")
        return(c(init, fin, f))
    }

asyncInit <- cmpfun(asyncInit)

#' Title
#'
#' @param sInit
#' @param update_matrix
#' @param maxT
#' @param nNode
#' @param randChoice
#' @param nState
#' @param nIter
#' @param stateList1
#'
#' @return
#' @export
#'
#' @examples
coherenceIter <-
    function(sInit,
             update_matrix,
             maxT = 1000,
             nNode = 1,
             randChoice = F,
             nState = 100,
             nIter = 10,
             stateList1 = NULL)
    {
        init <- sInit
        sInit <- str_split(sInit, "") %>% unlist
        sInit <- as.integer(sInit[2:(length(sInit) - 1)])
        sInit <- ifelse(sInit == 0,-1, 1)
        if (nNode > 1)
            randChoice <- T
        nNodes <- length(sInit)
        if (is.null(stateList1)) {
            pertChoice <- sapply(1:nState, function(x) {
                sample(1:nNodes, nNode)
            })
            if (nNode > 1)
                pertChoice <- t(pertChoice)
            else
                pertChoice <- matrix(pertChoice, ncol = 1)

            stateList <- apply(pertChoice, 1, function(x) {
                s <- sInit
                s[x] <- -1 * s[x]
                s
            })

            if (!randChoice)
                stateList <- sapply(1:nNodes, function(x) {
                    s <- sInit
                    s[x] <- -1 * s[x]
                    s
                })
        }

        else
            stateList <- stateList1

        df <- apply(stateList, 2, function(x) {
            sapply(1:nIter, function(i) {
                if (is.null(stateList1))
                    asyncInit(x, update_matrix, maxT, init)
                else
                    asyncInit(x, update_matrix, maxT)
            }) %>% t %>% data.frame %>% set_names(c("init", "fin", "flag"))
        }, simplify = F) %>% reduce(rbind.data.frame) %>%
            group_by(init, fin, flag) %>%
            summarise(Freq = n(), .groups = 'drop')
        dfSum <- df %>% group_by(init) %>% summarise(Sum = sum(Freq))
        df <- merge(df, dfSum, by = "init", all = T)
        df$Freq <- round(df$Freq / df$Sum, 4)
        df %>% select(-Sum)

    }
coherenceIter <- cmpfun(coherenceIter)


#' Title
#'
#' @param topoFile
#' @param maxT
#' @param nNode
#' @param randChoice
#' @param nState
#' @param nIter
#' @param write
#'
#' @return
#' @export
#'
#' @examples
coherence <- function(topoFile,
                      maxT = 1000,
                      nNode = 1,
                      randChoice = F,
                      nState = NULL,
                      nIter = 10,
                      write = T)
{
    print(nNode)
    net <- topoFile %>% str_remove(".topo")
    nodeOrder <- readLines(paste0(net, "_nodes.txt"))
    ls <- TopoToIntMat(topoFile)
    update_matrix <- ls[[1]]
    nodes <- ls[[2]]
    colnames(update_matrix) <- rownames(update_matrix) <- nodes
    update_matrix <- update_matrix[nodeOrder, nodeOrder]
    update_matrix <- 2 * update_matrix + diag(length(nodes))
    groupLabels <-
        readLines(str_replace(topoFile, ".topo", ".teams")) %>% str_split(",")
    statesDf <-
        read_csv(paste0(net, "_finFlagFreq.csv"),
                 col_types = cols(),
                 lazy = F) %>%
        filter(!is.na(Avg0), flag == 1)
    if (nrow(statesDf) == 0)
        return(NA)
    if (!is.null(nState) && nState < nrow(statesDf)) {
        statesDf <- statesDf %>% arrange(Avg0)
        nS <- round(nState / 2)
        print(nS)
        print(nrow(statesDf))
        print(nrow(statesDf) - nS + 1)
        statesDf <-
            statesDf[c(1:nS, (nrow(statesDf) - nS + 1):nrow(statesDf)), ]
    }
    else {
        nState <- nrow(statesDf)
    }
    states <- statesDf$states
    df <- lapply(states, function(x) {
        coherenceIter(x, update_matrix, maxT, nNode,
                      randChoice, nState, nIter)
    }) %>% reduce(rbind.data.frame)
    df$flag <- as.integer(df$flag)
    df$nNode <- nNode
    df$Hamming <- df %>% select(init, fin) %>% apply(1, hamming)
    df$initPhen <-
        df$init %>% EMTScoreCalc(groupLabels, nodeOrder) %>% LabellerPhen
    df$finPhen <-
        df$fin %>% EMTScoreCalc(groupLabels, nodeOrder) %>% LabellerPhen
    if (write) {
        DirectoryNav("CoherencesData")
        write_csv(df, paste0(net, "_", nNode, "R_coherence.csv"),
                  quote = "none")
    }

    else
        return(df)
}

coherence <- cmpfun(coherence)

#' Title
#'
#' @param nState
#'
#' @return
#' @export
#'
#' @examples
CoherenceSingleNode <- function(nState = NULL)
{
    topoFiles <- list.files(".", ".topo")
    nets <- topoFiles %>% str_remove(".topo")
    plan(multisession, workers = numThreads)
    dummy <- future_lapply(nets, function(x) {
        freqDf <-
            read_csv(paste0(x, "_finFlagFreq.csv"),
                     lazy = F,
                     col_types = cols())
        if ("coherence0" %in% colnames(freqDf)) {
            return()
        }
        print(x)
        df <-
            coherence(
                paste0(x, ".topo"),
                randChoice = F,
                nIter = 50,
                write = F,
                nState = nState
            )
        if (is.na(df)) {
            return()
        }
        df <-
            df[df$init == df$fin, ] %>% mutate(states = init, coherence0 = Freq) %>%
            select(states, coherence0)
        coherenceVec <- df$coherence0
        names(coherenceVec) <- df$states
        df <-
            read_csv(paste0(x, "_finFlagFreq.csv"),
                     show_col_types = F,
                     lazy = F) %>%
            mutate(coherence0 = coherenceVec[states])
        write_csv(df,
                  paste0(x, "_finFlagFreq.csv"),
                  quote = "none",
                  na = "")
    })
    future:::ClusterRegistry("stop")

}
CoherenceSingleNode <- cmpfun(CoherenceSingleNode)



#' Title
#'
#' @param topoFile
#' @param maxT
#' @param randChoice
#' @param nState
#' @param nIter
#' @param reps
#'
#' @return
#' @export
#'
#' @examples
CoherenceAllNode <-
    function(topoFile,
             maxT = 1000,
             randChoice = T,
             nState = 100,
             nIter = 10,
             reps = 1)
    {
        print(topoFile)
        net <- topoFile %>% str_remove(".topo")
        l <- readLines(paste0(net, "_nodes.txt")) %>% length
        dfList <- lapply(1:reps, function(rep) {
            # plan(multisession, workers = numThreads)
            df <- lapply(1:l, function(x) {
                coherence(topoFile, maxT, x, randChoice, nState, nIter, write = F)
            }) %>% reduce(rbind.data.frame) %>% mutate(Rep = rep)
            # future:::ClusterRegistry("stop")
            return(df)
        })
        if (reps == 1)
            df <- dfList[[1]]
        else
            df <- dfList %>% reduce(rbind.data.frame)
        DirectoryNav("PhenotypicTransition")
        write_csv(
            df,
            paste0(
                net,
                "_allNodeCoherence_nPert",
                nState,
                "_nIter",
                nIter,
                "_reps",
                reps,
                ".csv"
            ),
            quote = "none"
        )
        setwd("..")
        if (writeLog) {
            write_csv(logDf, "LogFile.csv", quote = "none")
        }
    }
CoherenceAllNode <- cmpfun(CoherenceAllNode)

#' Title
#'
#' @param net
#' @param topoFile
#'
#' @return
#' @export
#'
#' @examples
coherenceSignal <- function(net, topoFile = "wild.topo") {
    setwd(paste0(randRaw, "/", net))
    topoDf <- read.delim(topoFile, sep = "")
    nodes <- readLines(str_replace(topoFile, ".topo", "_nodes.txt"))

    ls <- TopoToIntMat(topoFile)
    update_matrix <- ls[[1]]
    rownames(update_matrix) <- colnames(update_matrix) <- ls[[2]]
    update_matrix <- update_matrix[nodes, nodes]
    update_matrix <- 2 * update_matrix + diag(length(nodes))
    groupLabels <-
        readLines(str_replace(topoFile, ".topo", ".teams")) %>% str_split(",")
    eNodes <- groupLabels[[1]]
    mNodes <- groupLabels[[2]]

    signals <- nodes[!(nodes %in% topoDf$Target)]
    outputs <- nodes[!(nodes %in% topoDf$Source)]
    freqDf <-
        read_csv(
            str_replace(topoFile, ".topo", "_finFlagFreq.csv"),
            lazy = F,
            col_types = cols()
        ) %>%
        filter(flag == 1) %>% filter(!is.na(Avg0))
    states <- freqDf %>%
        separate(
            col = states,
            into = c("dummy1", "dummy2", nodes, "dummy3"),
            sep = ""
        ) %>%
        select(-contains("dummy")) %>% select(all_of(nodes), Phenotype)

    df <- apply(states, 1, function(s) {
        p <- s[length(s)]
        # print(p)
        s1 <- s[-length(s)] %>% as.numeric
        s <- ifelse(s1 == 0,-1, 1)
        names(s) <- nodes
        stateList <- sapply(signals, function(i) {
            s1 <- s
            s1[i] <- -1 * s1[i]
            s1
        })
        d <-
            coherenceIter(s, update_matrix, stateList1 = stateList) %>%
            mutate(
                initPhen = p,
                finPhen = fin %>% EMTScoreCalc(groupLabels, nodes) %>% LabellerPhen,
                initState = ifelse(s == -1, 0, 1) %>%
                    paste0(collapse = "") %>%
                    paste0("'", ., "'")
            )

        d$signal <-
            nodes[d %>% select(initState, init) %>% apply(1, diffPos) - 1]
        d <- d %>% select(initState, fin) %>% apply(1, function(x) {
            n <- nodes[stateCompare(x)]
            epi <- sum(eNodes %in% n) / length(eNodes)
            mes <- sum(mNodes %in% n) / length(mNodes)
            c(epi, mes)
        }) %>% t %>% data.frame %>% set_names(c("Epithelial", "Mesenchymal")) %>%
            cbind.data.frame(d, .) %>%
            select(initState,
                   initPhen,
                   signal,
                   Freq,
                   Epithelial,
                   Mesenchymal)
        d
    }) %>% reduce(rbind.data.frame)
    DirectoryNav("signalCoherence")
    write_csv(df,
              str_replace(topoFile, ".topo", "_signalCoherence.csv"),
              quote = "none")
}
