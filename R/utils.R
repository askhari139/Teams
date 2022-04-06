RemoveAllFiles <- function() {
    filz <- list.files(".")
    sapply(filz, file.remove)
}

DirectoryNav <- function(d) {
    if(!dir.exists(d))
        dir.create(d)
    setwd(d)
}

bmodelFreqFileReader <- function(topoFile) {
    df <- read_csv(str_replace(topoFile, ".topo", "_finFlagFreq.csv"), lazy = F) %>%
        filter(flag == 1)
    return(df)
}

skewness<-function(x, finite=TRUE){
    n=length(x)
    S=(1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
    if(finite==FALSE){
        S=S
    }else{
        S=S*(sqrt(n*(n-1)))/(n-2)
    }
    return(S)
}

kurtosis<-function(x, finite){
    n=length(x)
    K=(1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2) - 3
    if(finite==FALSE){
        K=K
    }
    else{
        K=((n-1)*((n+1)*K - 3*(n-1))/((n-2)*(n-3))) +3
    }
    return(K)
}

bimodality_coefficient<-function(x, finite=TRUE,...){
    if(finite==TRUE){
        G=skewness(x,finite)
        sample.excess.kurtosis=kurtosis(x,finite)
        K=sample.excess.kurtosis
        n=length(x)
        B=((G^2)+1)/(K+ ((3*((n-1)^2))/((n-2)*(n-3))))
    }
    else{
        G=skewness(x,FALSE)
        K=kurtosis(x,FALSE)
        B=((G^2)+1)/(K)
    }
    return(B)
}

bimodality_coefficient <- function(vec)
{
    vec <- na.omit(vec)
    s <- sd(vec)
    m <- mean(vec)
    n <- length(vec)
    sk <- sum((vec-mean(vec))^3)/((length(vec)-1)*sd(vec)^3)
    ku <- sum((vec-mean(vec))^4)/((length(vec)-1)*sd(vec)^4) -3
    (sk^2 + 1)/(ku + 3*((n-1)^2)/((n-2)*(n-3)))
}

cfgAnalysis <- function(net)
{
    cfgFile <- paste0(net, ".cfg")
    cfgDat <- readLines(cfgFile)
    numEdges <- cfgDat[30] %>% str_extract("\\d+") %>% as.integer
    numNodes <- cfgDat[31] %>% str_extract("\\d+") %>% as.integer
    nodes <- cfgDat[32:(31+numNodes)] %>% str_remove("^\\d+?\t") %>%
        str_trim %>% str_replace_all(regex("\\W+"), "")
    topoDat <- cfgDat[(32+numNodes):length(cfgDat)] %>% str_split("\t") %>%
        reduce(rbind.data.frame) %>% set_names(c("ID", "S", "Tar", "Type")) %>%
        mutate(Type = ifelse(Type == "2", -1, 1), S = as.integer(S),
               Tar = as.integer(Tar))
    if (!file.exists(paste0(net, ".topo")))
    {
        topoDf <- topoDat %>% select(-ID) %>% mutate(S = nodes[S], Tar = nodes[Tar]) %>%
            set_names(c("Source", "Target", "Type")) %>%
            mutate(Type = ifelse(Type == -1, 2, 1))
        write_delim(topoDf, paste0(net, ".topo"), delim = " ", quote = "none")
    }

    return(list(nodes, topoDat))
}

frustCalcRAC <- function(state, topoDat)
{
    sum(state[topoDat$S]*state[topoDat$Tar] == topoDat$Type)/nrow(topoDat)
}
frustCalcRAC <- frustCalcRAC %>% cmpfun


discretize <- function(net)
{
    setwd(paste0(RACIPE_WT, "/", net))
    ls <- cfgAnalysis(net)
    nodes <- ls[[1]]
    topoDat <- ls[[2]]
    solutionDf <- paste0(net, "_solution.dat") %>% read_delim(delim = "\t") %>%
        set_names(c("ParIndex", "nStates", "Count", nodes)) %>%
         mutate(Count = Count/max(Count)) %>%
        select(all_of(nodes), Count)
    dots <- lapply(nodes, as.symbol)
    states <- solutionDf %>% select(-Count) %>% sapply(function(x){
        y <- (x-mean(x))/sd(x)
        ifelse(y>0, 1, -1)
    }) %>% data.frame %>% set_names(nodes) %>% mutate(Count = solutionDf$Count) %>%
        group_by(across(nodes)) %>% summarise(Frequency = sum(Count))
    states$Frequency <- states$Frequency/sum(states$Frequency)
    states$Frustration <- states %>% select(all_of(nodes)) %>%
                   apply(1, function(x){frustCalcRAC(x, topoDat)})
    statesDf <- states %>%
        unite(State, all_of(nodes), sep = "") %>%
        mutate(State = paste0("'",State %>% str_replace_all("-1", "0"),"'"))
    write_csv(statesDf, paste0(net, "_discreteStates.csv"), quote = "none")
    print(net)
}

SecondarySignals <- function(topoDf, sig) {
    targets <- topoDf %>% filter(Source == sig) %>%
        select(Target) %>% unlist
    secSigs <- c()

}

getEMSONodes <- function(topoFile)
{
    wd <- getwd()
    ls <- TopoToIntMat(topoFile)
    intMat <- ls[[1]]
    nodes <- ls[[2]]
    colnames(intMat) <- rownames(intMat) <- nodes
    signal <- which(apply(intMat, 2, function(x){all(x==0)}))
    output <- which(apply(intMat, 1, function(x){all(x==0)}))
    secondary_signal <- which(apply(intMat, 2, function(x){
        if (length(signal) !=0)
            all(x[-signal] == 0)
        else
            F
    }))
    secondary_output <- which(apply(intMat, 1, function(x){
        if (length(output) != 0)
            all(x[-output] == 0)
        else
            F
    }))
    ls <- readLines(str_replace(topoFile, ".topo", ".teams")) %>% str_split(",")
    names(ls[[1]]) <- paste0("E", 1:length(ls[[1]]))
    names(ls[[2]]) <- paste0("M", 1:length(ls[[2]]))
    groups <- ls %>% unlist
    sigs <- unique(nodes[c(signal, secondary_signal)])
    names(sigs) <- paste0("S", 1:length(sigs))
    outs <- unique(nodes[c(output, secondary_output)])
    if (length(outs) > 0)
        names(outs) <- paste0("O", 1:length(outs))
    nodes <- unique(c(groups, sigs, outs))
    nodesInGroups <- which(nodes %in% groups)
    dummy <- sapply(1:length(nodes), function(x){
        n <- nodes[x]
        if(n %in% groups)
            names(nodes)[x] <<- names(groups[groups == nodes[x]])
        if(n %in% sigs)
            names(nodes)[x] <<- names(sigs[sigs == nodes[x]])
        if(n %in% outs)
            names(nodes)[x] <<- names(outs[outs == nodes[x]])
    })

    setwd(wd)
    return(nodes)
}

correlGrob <- function(df, x, y, xPos = NULL, yPos = NULL, method = "pearson")
{
    corr <- cor.test(df[[x]], df[[y]], method = method)
    pVal <- ifelse(corr$p.value < 0.05, "*", "")
    xPos <- ifelse(!is.null(xPos), xPos, 0.5)
    yPos <- ifelse(!is.null(yPos), yPos, 0.9)
    grob <- grobTree(textGrob(paste0("\u03c1 : ", round(corr$estimate, 2), pVal),
                              x=xPos,  y=yPos, hjust=0,
                              gp=gpar(col="black", fontsize=18, fontface="bold")))
    return(grob)
}



multiFactorCorrelation <- function(df,x, y,z, label = T, method = "pearson")
{
    facts <- unique(df[[x]])
    s <- sapply(facts, function(f){
        d <- df[df[[x]] == f, ] %>% select(z,y)
        d <- cor.test(d[[y]], d[[z]], method = method)
        if(label)
            paste0("\u03c1: ",round(d$estimate,2), ifelse(d$p.value < 0.05, "*", ""))
        else
            c(d$estimate, d$p.value)
    })
    if (label)
    {
        names(s) <- facts
    }
    else
    {
        s <- s %>% t %>% data.frame %>% set_names(c("Correlation", "pValue")) %>%
            mutate(Factors = facts)
    }
    s

}

multiFactorCorrelationAnova <- function(df,x, y,z, label = T, method = "pearson")
{
    facts <- unique(df[[x]])
    s <- sapply(facts, function(f){
        d <- df[df[[x]] == f, ] %>% select(z,y)
        d <- cor(d[[y]], d[[z]], method = method)
        d1 <- df %>% select(all_of(c(x,y,z))) %>% set_names(c("X", "Y", "Z"))
        p <-
        if(label)
            paste0("\u03c1: ",round(d$estimate,2), ifelse(d$p.value < 0.05, "*", ""))
        else
            c(d$estimate, d$p.value)
    })
    if (label)
    {
        names(s) <- facts
    }
    else
    {
        s <- s %>% t %>% data.frame %>% set_names(c("Correlation", "pValue")) %>%
            mutate(Factors = facts)
    }
    s

}


