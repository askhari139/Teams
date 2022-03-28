TopoToIntMat <- function(topoFile, plotOut = F, nOrder = NULL) {
  df <- read.delim(topoFile, sep = " ", stringsAsFactors = F)
  df <- df %>%
    mutate(Type = ifelse(Type == 2, -1, 1))

  nodes <- unique(c(df$Source, df$Target)) %>%
    sort(decreasing = T)
  n_nodes <- length(nodes)
  intmat <- rep(0, n_nodes * n_nodes) %>%
    matrix(ncol = n_nodes)
  df1 <- df %>%
    mutate(Source = sapply(Source, function(x) {
      which(nodes == x)
    })) %>%
    mutate(Target = sapply(Target, function(x) {
      which(nodes == x)
    }))
  # browser()
  dummy <- apply(df1, 1, function(x) {
    # browser()
    i <- x[1]
    j <- x[2]
    k <- x[3]
    intmat[i, j] <<- k
  })
  if (plotOut) {
      if(is.null(nOrder)) {
          nOrder <- getEMSONodes(topoFile)
      }
      df <- data.frame(intmat) %>%
          set_names(nodes) %>%
          mutate(nodes1 = nodes) %>%
          gather(key = "Nodes", value = "Influence", -nodes1) %>%
          mutate(nodes1 = factor(nodes1, levels = nOrder), Nodes = factor(Nodes, levels = nOrder))
      ggplot(df, aes(x = Nodes, y = nodes1, fill = Influence)) + geom_tile() +
          theme_Publication() + scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) +
          labs(fill = "Edge") +
          theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                axis.text.x = element_text(angle = 90), legend.position = "right",
                legend.direction = "vertical", legend.key.height = unit(0.5, "in"))
      directoryNav("MatrixPlots")
      ggsave(paste0(str_replace(topoFile, ".topo", "_interaction.png")),
             width = 7, height = 6)
      setwd("..")

  }
  return(list(intmat, nodes))
}

ComputePowerMatrix <- function(mat, power) {
  res <- mat
  if (power == 1) {
    return(res)
  }
  for (i in 2:power) {
    res <- res %*% mat
  }
  return(res)
}
ComputePowerMatrix <- cmpfun(ComputePowerMatrix)

InfluenceMatrix <- function(net, intmat, nodes, lmax = 10, write = T) {
  intmax <- intmat
  intmax[which(intmax == -1)] <- 1
  res <- 0
  for (l in 1:lmax) {
    intM <- ComputePowerMatrix(intmat, l)
    maxM <- ComputePowerMatrix(intmax, l)
    r1 <- intM / maxM
    r1[is.nan(r1)] <- intM[is.nan(r1)]
    res <- res + r1
  }
  res <- res / lmax

  nodes <- nodes

  influence_mat <- res
  colnames(influence_mat) <- rownames(influence_mat) <- nodes
  signal <- which(apply(intmat, 2, function(x) {
    all(x == 0)
  }))
  output <- which(apply(intmat, 1, function(x) {
    all(x == 0)
  }))
  secondary_signal <- which(apply(intmat, 2, function(x) {
    if (length(signal) != 0) {
      all(x[-signal] == 0)
    } else {
      F
    }
  }))
  secondary_output <- which(apply(intmat, 1, function(x) {
    if (length(output) != 0) {
      all(x[-output] == 0)
    } else {
      F
    }
  }))
  nonEssentials <- c(signal, output, secondary_signal, secondary_output)
  if (length(nonEssentials)) {
    influence_reduced <- influence_mat[
      -nonEssentials,
      -nonEssentials
    ]
    nodes_reduced <- nodes[-nonEssentials] %>%
      c() %>%
      str_replace_all(regex("\\W+"), "")
  } else {
    influence_reduced <- influence_mat
    nodes_reduced <- nodes %>% str_replace_all(regex("\\W+"), "")
  }
  if (length(nodes_reduced) < 2) {
    return()
  }
  rownames(influence_reduced) <- colnames(influence_reduced) <- nodes_reduced
  if (write) {
    directoryNav("Influence")
    write.csv(influence_mat, paste0(net, "_fullInfl.csv"))
    write.csv(influence_reduced, paste0(net, "_reducedInfl.csv"))
    setwd("..")
  }
  influence_reduced
}
InfluenceMatrix <- cmpfun(InfluenceMatrix)


GroupStrength <- function(topoFile, pathLength = 10, plotOut = F, getTeams = F) {
  ls <- topo_to_int_mat(topoFile)
  intmat <- ls[[1]]
  nodes <- ls[[2]]


  inflMat <- InfluenceMatrix(intmat, nodes, pathLength)

  nodes <- rownames(inflMat)
  df <- inflMat
  df1 <- apply(df, 2, function(x) {
    ifelse(x > 0, 1, -1)
  })
  df1 <- cbind(df1, t(df1))
  hc <- hclust(dist(df1))
  clust <- cutree(hc, 2)
  g1 <- nodes[clust == 1] %>% sort()
  g2 <- nodes[clust == 2] %>% sort()
  if (g1[length(g1)] > g2[1]) {
    g0 <- g1
    g1 <- g2
    g2 <- g0
  }
  nOrder <- c(g1, g2)
  if (getTeams) {
    l <- list(g1, g2)
    mirdetect <- c(sum(str_detect(g1 %>% str_to_upper, "MIR")),
                   sum(str_detect(g2 %>% str_to_upper, "MIR")))
    egroup <- which.max(mirdetect)
    mgroup <- ifelse(egroup == 1, 2, 1)
    names(l)[c(egroup, mgroup)] <- c("E", "M")
    return(l)
  }
  df2 <- data.frame(df) %>%
    mutate(nodes1 = nodes) %>%
    gather(key = "Nodes", value = "Influence", -nodes1) %>%
    mutate(nodes1 = factor(nodes1, levels = nOrder), Nodes = factor(Nodes, levels = nOrder))
  if (plotOut) {
      nOrder <- c(g1,g2)
      ggplot(df2, aes(x = Nodes, y = nodes1, fill = Influence)) + geom_tile() +
          theme_Publication() + scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) +
          theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                axis.text.x = element_text(angle = 90), legend.position = "right",
                legend.direction = "vertical", legend.key.height = unit(0.5, "in"))

      directoryNav("MatrixPlots")
      ggsave(str_replace(topoFile, ".topo", "_group.png"), width = 7, height = 6)
      setwd("..")
      return()
  }

  g11 <- df[g1, g1] %>%
    as.vector() %>%
    mean()
  g22 <- df[g2, g2] %>%
    as.vector() %>%
    mean()
  g12 <- df[g1, g2] %>%
    as.vector() %>%
    mean()
  g21 <- df[g2, g1] %>%
    as.vector() %>%
    mean()
  c(g11, g22, g12, g21) %>%
    abs() %>%
    mean() %>%
      c(g11, g22, g12, g21, .)
}


GroupStrengthAll <- function(net) {
    wd <- getwd()
    topoFiles <- list.files(".", ".topo")
    df <- sapply(topoFiles, function(topoFile) {
        groupStrength(topoFile)
    }) %>% t %>%
        data.frame %>%
        set_names(c("G11", "G22", "G12", "G21", "Gs")) %>%
        mutate(Network = topoFiles %>% str_remove(".topo$"))
    directoryNav("CompiledData")
    write_csv(df, paste0(net, "_GroupStrengths.csv"), quote = "none")
}

