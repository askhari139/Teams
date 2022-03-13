Swap <- function(x, p)
{
    d <- x[p[1]]
    x[p[1]] <- x[p[2]]
    x[p[2]] <- d
    return(x)
}
Swap <- cmpfun(Swap)

UniquePerms2 <- function(x, max = 10, nSwap = 10){
    # browser()
    N <- length(x)
    n1 <- sum(x == 1)
    n2 <- sum(x == 2)
    nVec <- c(n1,1)
    if (n1>n2) nVec <- c(n2,2)
    notOne <- ifelse(nVec[2] == 1, 2, 1)
    if (max == 0)
        max <- factorial(N)/(factorial(n1)*factorial(n2))
    if (max > 500)
        max <- 500
    perms <- rep(0,max+1)
    perms[1] <- sum(2^(which(x == nVec[2])))
    lists <- matrix(rep(notOne, max*N), ncol = N)
    count <- 1
    while(count<max){
        nSamp <- x
        sapply(1:nSwap, function(x){
            p <- sample(1:N, 2)
            nSamp <<- Swap(nSamp, p)
        })
        nNew <- nSamp
        nSamp <- which(nNew == nVec[2])
        perm <- sum(2^nSamp)
        if (!any(perms == perm))
        {
            perms[count+1] <- perm
            lists[count,] <- nNew
            count <- count + 1
            print(perm)
        }
    }
    return(lists)
}
UniquePerms2 <- cmpfun(UniquePerms2)

RandomNetworks <- function(numRand = 500, nSwap = 10) {
    wd <- getwd()
    setwd(topoFileFolder)
    topoFiles <- list.files(".", pattern = ".topo")
    topoDf <- lapply(topoFiles, read.delim, sep = " ") %>% set_names(topoFiles)
    
    for (topoFile in topoFiles){
        net <- str_remove(topoFile, ".topo")
        df <- topoDf[[topoFile]]
        wt <- df$Type
        n2 <- sum(wt == 2)
        perMax <- choose(length(wt), n2)
        if (perMax < numRand)
        {
            numRand <- perMax
        }
        setwd(randRaw)
        if(!dir.exists(net))
            dir.create(net)
        setwd(net)
        write_delim(df, "wild.topo", delim = " ", quote = "none")
        onetwo <- df[, 1:2]
        rand_orders <- UniquePerms2(wt, max = numRand, nSwap = nSwap)
        dummy <- sapply(1:nrow(rand_orders), function(x){
            y <- rand_orders[x,]
            df1 <- cbind.data.frame(onetwo, y)
            colnames(df1) <- c("Source", "Target", "Type")
            write_delim(df1, paste0(name, "_", x, ".topo"), delim = " ")
        })
        setwd("..")
    }
    setwd(wd)
}

EdgeDeletion <- function() {
    wd <- getwd()
    setwd(topoFileFolder)
    topoFiles <- list.files(".", pattern = ".topo")
    topoDf <- lapply(topoFiles, read.delim, sep = " ") %>% set_names(topoFiles)
    
    for (topoFile in topoFiles) {
        net <- str_remove(topoFile, ".topo")
        df <- topoDf[[topoFile]]
        setwd(edgeDel)
        if(!dir.exists(net))
            dir.create(net)
        setwd(net)
        sapply(1:nrow(df), function(x) {
            write_delim(df, "wild.topo", delim = " ", quote = "none")
            dfNew <- df[-x, ]
            nam <- paste0(net, df$Source[x], "_", df$Target[x], "_del.topo")
            write_delim(dfNew, nam, delim = " ", quote = "none")
            dfNew <- df
            dfNew$Type[x] <- ifelse(dfNew$Type[x] == 1, 2, 1)
            nam <- paste0(net, df$Source[x], "_", df$Target[x], "_change.topo")
            write_delim(dfNew, nam, delim = " ", quote = "none")
        })
    }
}


SimulateBoolean <- function(rand = T, edge = T) {
    if(rand) {
        setwd(randRaw)
        sapply(netList, function(net) {
            setwd(net)
            file.copy(paste0(simPackage, "/script.jl"), ".")
            julia_source("script.jl")
            setwd("..")
        })
    }
    
    if(edge) {
        setwd(edgeDel)
        sapply(netList, function(net) {
            setwd(net)
            file.copy(paste0(simPackage, "/script.jl"), ".")
            julia_source("script.jl")
            setwd("..")
        })
    }
}