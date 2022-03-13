topo_to_int_mat <- function(topo_file) {
    # print(topo_file)
    df <- read.delim(topo_file, sep = " ", stringsAsFactors = F)
    if (ncol(df) != 3) {
        df <- read.delim(topo_file, stringsAsFactors = F)
    }
    # browser()
    colnames(df) <- c("Source", "Target", "Type")
    df <- df %>% 
        mutate(Source = str_remove_all(Source, "\\s")) %>%
        mutate(Target = str_remove_all(Target, "\\s")) %>%
        mutate(Type = ifelse(Type == 2, -1, 1))
    
    nodes <- unique(c(df$Source, df$Target)) %>% 
        sort(decreasing = T)
    n_nodes <- length(nodes)
    intmat <- rep(0, n_nodes * n_nodes) %>% 
        matrix(ncol = n_nodes)
    df1 <- df %>% 
        mutate(Source = sapply(Source, function(x) {which(nodes == x)})) %>% 
        mutate(Target = sapply(Target, function(x) {which(nodes == x)}))
    # browser()
    dummy <- apply(df1, 1, function(x) {
        # browser()
        i <- x[1]
        j <- x[2]
        k <- x[3]
        intmat[i,j] <<- k
    })
    return(list(intmat, nodes))
}

compute_power_matrix <- function(mat, power) {
    res <- mat
    if (power == 1)
    {
        return(res)
    }
    for (i in 2:power) {
        res <- res %*% mat
    }
    return(res)
}
compute_power_matrix <- cmpfun(compute_power_matrix)

influence_matrix <- function(intmat, nodes, lmax = 10, write = T) {
    intmax <- intmat
    intmax[which(intmax == -1)] <- 1
    res <- 0
    for (l in 1:lmax) {
        intM <- compute_power_matrix(intmat, l)
        maxM <- compute_power_matrix(intmax, l)
        r1 <- intM / maxM 
        r1[is.nan(r1)] <- intM[is.nan(r1)]
        res <- res + r1
    }
    res <- res / lmax
    
    nodes <- nodes %>% str_replace_all(regex("\\W+"), "")
    
    influence_mat <- res
    colnames(influence_mat) <- rownames(influence_mat) <- nodes
    signal <- which(apply(intmat, 2, function(x){all(x==0)}))
    output <- which(apply(intmat, 1, function(x){all(x==0)}))
    secondary_signal <- which(apply(intmat, 2, function(x){
        if (length(signal) !=0)
            all(x[-signal] == 0)
        else
            F
    }))
    secondary_output <- which(apply(intmat, 1, function(x){
        if (length(output) != 0)
            all(x[-output] == 0)
        else
            F
    }))
    nonEssentials <- c(signal, output, secondary_signal, secondary_output)
    if(length(nonEssentials))
    {
        influence_reduced <- influence_mat[-nonEssentials, 
                                           -nonEssentials]
        nodes_reduced <- nodes[-nonEssentials] %>% c  %>% str_replace_all(regex("\\W+"), "")
    }
    else
    {
        influence_reduced <- influence_mat
        nodes_reduced <- nodes %>% str_replace_all(regex("\\W+"), "")
    }
    if (length(nodes_reduced) < 2)
        return()
    rownames(influence_reduced) <- colnames(influence_reduced) <- nodes_reduced
    if (write)
    {
        
    }
    influence_reduced
}
influence_matrix <- cmpfun(influence_matrix)


gsPathLength <- function(topoFile, pathLength = 10)
{
    ls <- topo_to_int_mat(topoFile)
    intmat <- ls[[1]]
    nodes <- ls[[2]]
    
    
    inflMat <- influence_matrix(intmat, nodes, pathLength)
    
    nodes <- rownames(inflMat)
    df <- inflMat
    df1 <- apply(df, 2, function(x){
        ifelse(x > 0, 1, -1)
    })
    df1 <- cbind(df1, t(df1))
    hc <- hclust(dist(df1))
    clust <- cutree(hc, 2)
    g1 <- nodes[clust == 1] %>% sort
    g2 <- nodes[clust == 2] %>% sort
    if(g1[length(g1)]> g2[1])
    {
        g0 <- g1
        g1 <- g2
        g2 <- g0
    }
    nOrder <- c(g1,g2)
    df2 <- data.frame(df) %>% mutate(nodes1 = nodes) %>%
        gather(key = "Nodes", value = "Influence", -nodes1) %>%
        mutate(nodes1 = factor(nodes1, levels = nOrder), Nodes = factor(Nodes, levels = nOrder))
    
    g11 <- df[g1,g1] %>% as.vector %>% mean
    g22 <- df[g2,g2] %>% as.vector %>% mean
    g12 <- df[g1,g2] %>% as.vector %>% mean
    g21 <- df[g2,g1] %>% as.vector %>% mean
    c(g11, g22, g12, g21) %>% abs %>% mean
    
}
