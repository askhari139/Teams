#' Title
#'
#' @param x
#' @param p
#'
#' @return
#' @export
#'
#' @examples
Swap <- function(x, p) {
    d <- x[p[1]]
    x[p[1]] <- x[p[2]]
    x[p[2]] <- d
    return(x)
}
Swap <- cmpfun(Swap)

#' Title
#'
#' @param x
#' @param max
#' @param nSwap
#'
#' @return
#' @export
#'
#' @examples
UniquePerms2 <- function(x, max = 10, nSwap = 10) {
    N <- length(x)
    n1 <- sum(x == 1)
    n2 <- sum(x == 2)
    nVec <- c(n1, 1)
    if (n1 > n2)
        nVec <- c(n2, 2)
    notOne <- ifelse(nVec[2] == 1, 2, 1)
    maxNets <- factorial(N) / (factorial(n1) * factorial(n2))
    if (is.nan(maxNets))
        maxNets <- 1000
    if (max == 0 || max > maxNets) {
        max <- maxNets
    }
    if (max > 500) {
        max <- 500
    }
    perms <- rep(0, max + 1)
    perms[1] <- sum(2 ^ (which(x == nVec[2])))
    lists <- matrix(rep(notOne, max * N), ncol = N)
    count <- 1
    while (count < max) {
        nSamp <- x
        sapply(1:nSwap, function(x) {
            p <- sample(1:N, 2)
            nSamp <<- Swap(nSamp, p)
        })
        nNew <- nSamp
        nSamp <- which(nNew == nVec[2])
        perm <- sum(2 ^ nSamp)
        if (!any(perms == perm)) {
            perms[count + 1] <- perm
            lists[count,] <- nNew
            count <- count + 1
            print(perm)
        }
    }
    return(lists)
}
UniquePerms2 <- cmpfun(UniquePerms2)

#' Title
#'
#' @param numRand
#' @param nSwap
#'
#' @return
#' @export
#'
#' @examples
RandomNetworks <- function(numRand = 500, nSwap = 10) {
    wd <- getwd()
    setwd(topoFileFolder)
    topoFiles <- list.files(".", pattern = ".topo")
    topoDf <-
        lapply(topoFiles, read.delim, sep = " ") %>% set_names(topoFiles)

    for (topoFile in topoFiles) {
        net <- str_remove(topoFile, ".topo")
        df <- topoDf[[topoFile]]
        wt <- df$Type
        n2 <- sum(wt == 2)
        perMax <- choose(length(wt), n2)
        if (perMax < numRand) {
            numRand <- perMax
        }
        setwd(randRaw)
        if (!dir.exists(net)) {
            dir.create(net)
        }
        setwd(net)
        write_delim(df, "wild.topo", delim = " ", quote = "none")
        onetwo <- df[, 1:2]
        rand_orders <- UniquePerms2(wt, max = numRand, nSwap = nSwap)
        dummy <- sapply(1:nrow(rand_orders), function(x) {
            y <- rand_orders[x,]
            df1 <- cbind.data.frame(onetwo, y)
            colnames(df1) <- c("Source", "Target", "Type")
            write_delim(df1, paste0(net, "_", x, ".topo"), delim = " ")
        })
        setwd("..")
    }
    setwd(wd)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
EdgeDeletion <- function() {
    wd <- getwd()
    setwd(topoFileFolder)
    topoFiles <- list.files(".", pattern = ".topo")
    topoDf <-
        lapply(topoFiles, read.delim, sep = " ") %>% set_names(topoFiles)

    for (topoFile in topoFiles) {
        net <- str_remove(topoFile, ".topo")
        df <- topoDf[[topoFile]]
        setwd(edgeDel)
        if (!dir.exists(net)) {
            dir.create(net)
        }
        setwd(net)
        sapply(1:nrow(df), function(x) {
            write_delim(df,
                        "wild.topo",
                        delim = " ",
                        quote = "none")
            dfNew <- df[-x,]
            nam <-
                paste0(net, "_", df$Source[x], "_", df$Target[x], "_del.topo")
            write_delim(dfNew, nam, delim = " ", quote = "none")
            dfNew <- df
            dfNew$Type[x] <- ifelse(dfNew$Type[x] == 1, 2, 1)
            nam <-
                paste0(net,
                       "_",
                       df$Source[x],
                       "_",
                       df$Target[x],
                       "_change.topo")
            write_delim(dfNew, nam, delim = " ", quote = "none")
        })
    }
}


#' Title
#'
#' @param topoFile
#'
#' @return
#' @export
#'
#' @examples
findMax <- function(topoFile) {
    topoOriginal <- list.files(".", ".topo")
    net <- topoFile %>% str_remove(".topo")
    topoDf <- read.delim(topoFile, sep = "")
    GWT <- GroupStrength(topoFile)[5]
    if (is.null(GWT) ||
        GWT == 0)
        return("Stop")
    GsChange <- sapply(1:nrow(topoDf), function(x) {
        df <- topoDf[-x,]
        topo <- paste0(net, "_", x, ".topo")
        write_delim(df, topo, delim = " ", quote = "none")
        G <- GroupStrength(topo)[5]
        GWT - G
    })
    n <- which.max(GsChange)
    topoFile <- paste0(net, "_", n, ".topo")
    topoOriginal <- c(topoOriginal, topoFile)
    topoFiles <- list.files(".", ".topo")
    topoFiles <- topoFiles[-which(topoFiles %in% topoOriginal)]
    sapply(topoFiles, file.remove)
    return(topoFile)
}


#' Title
#'
#' @param nMax
#'
#' @return
#' @export
#'
#' @examples
CausationNetworks <- function(nMax = 20) {
    wd <- getwd()
    setwd(topoFileFolder)
    topoFiles <- list.files(".", pattern = ".topo")
    topoDf <-
        lapply(topoFiles, read.delim, sep = " ") %>% set_names(topoFiles)
    sapply(topoFiles, function(topoFile) {
        setwd(gsCausation)
        net <- topoFile %>% str_remove(".topo")
        DirectoryNav(net)
        RemoveAllFiles()
        write_delim(topoDf[[topoFile]],
                    "wild.topo",
                    delim = " ",
                    quote = "none")
        topo <- "wild.topo"
        LogFileGen()
        sapply(1:nMax, function(i) {
            # browser()
            if (topo == "Stop")
                return()
            topo <<- findMax(topo)
        })
        setwd("Influence")
        RemoveAllFiles()
        setwd("..")
    })

}

LogFileGen <- function() {
    topoFiles <- list.files(".", ".topo")
    df <- data.frame(
        Files = topoFiles,
        Influence = "No",
        Gs = "No",
        Coherence = "No",
        MultiNodeCoherence = "No",
        InfluencePlot = "No",
        Correlation = "No",
        CorrelationPlot = "No",
        Simulated = "No",
        TeamComposition = "No",
        Strength = "No",
        Scores = "No",
        InteractionPlot = "No"
    )
    write.csv(df, "LogFile.csv", row.names = F)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
simulation <- function() {
    os <- .Platform$OS.type
    script <- "script.jl"
    if (os == "windows")
        script <- "scriptWindows.jl"
    file.copy(paste0(simPackage, "/", script), ".", overwrite = T)
    if (os != "windows") {
        Sys.setenv(JULIA_NUM_THREADS = as.character(numThreads))
        command <- "julia script.jl"
        system(command)
    }
    else {
        topoFiles <- list.files(".", ".topo$")
        size <- floor(length(topoFiles) / numThreads)
        topoList <- lapply(1:numThreads, function(x) {
            k <- (x - 1) * size
            id <- 1:size + k
            if (size + k > length(topoFiles))
                id <- id[1]:(length(topoFiles))
            topoFiles[id] %>% paste0(collapse = " ")
        })
        Sys.setenv(JULIA_NUM_THREADS = as.character(1))
        plan(multisession, workers = numThreads)
        simulater <- future_lapply(topoList, function(x) {
            command <- paste0("julia ", script, " ", x)
            system(command)
        })
        future:::ClusterRegistry("stop")
    }
    logDf <- read.csv("LogFile.csv")
    simulated <-
        list.files(".", "_finFlagFreq.csv") %>% str_replace("_finFlagFreq.csv", ".topo")
    logDf$Simulated[logDf$Files %in% simulated] <- "Yes"
    write.csv(logDf, "LogFile.csv", row.names = F)

}

#' Title
#'
#' @param rand
#' @param edge
#' @param cause
#'
#' @return
#' @export
#'
#' @examples
SimulateBoolean <- function(rand = F,
                            edge = F,
                            cause = F) {
    if (rand) {
        setwd(randRaw)
        sapply(netList, function(net) {
            setwd(net)
            simulation()
            setwd("..")
        })
    }

    if (edge) {
        setwd(edgeDel)
        sapply(netList, function(net) {
            setwd(net)
            simulation()
            setwd("..")
        })
    }

    if (cause) {
        setwd(gsCausation)
        sapply(netList, function(net) {
            setwd(net)
            simulation()
            setwd("..")
        })
    }

    if (!edge && !rand && !cause)
    {
        simulation()
    }
}
